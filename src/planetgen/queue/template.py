# planetgen/queue/template.py

"""
A warm template and the jobs forked from it (OPS.226, D1 of `docs/design/operations/job-engine.md`, section 3.3).

A **template** is a process that imports the code, builds one warm set (the bright-star sampling tables of `units.warm_star_tables`, say, which every worker builds for itself on main and
which takes about a minute and 1.5 GB), holds no connection and runs one thread, so it is safe to fork. **Every job is a fresh fork of it**: nothing waits holding tables, and the
tables are shared by every job the template has forked. The template is the parent of every job, so it is told how each one ended (its exit code, or the signal that killed it) and reports
it at once; the supervisor holds a kernel handle (`pidfd`) on each job as well, so the end of a job is known even if the template is the one that died.

This module is the template process and the job process. `jobhost` is the supervisor's side of it.

* **The channel.** The supervisor and the template talk over one `SOCK_SEQPACKET` socket pair, a JSON object a message. The supervisor sends `{"op": "fork", "id", "runner", "job"}` and
  `{"op": "stop"}`; the template sends `{"event": "ready"}`, `{"event": "forked", "id", "pid"}` (with the supervisor's end of the job's own pipe and a `pidfd` of the job attached, `socket.send_fds`; the template opens the `pidfd` itself, while the job is still its child to reap),
  `{"event": "exit", "id", "pid", "code", "signal"}` and `{"event": "error", "id", "error"}`.
* **The pipe.** Each job has a pipe of its own to the supervisor (`Link`): the job says what it is doing (`link.send`) and hears what it is asked (`link.requests`, `link.stop_requested`). When the
  supervisor goes, the pipe closes: the job notices at its next safe point (`link.check`), finishes the commit it is in and exits (`EXIT_SUPERVISOR_GONE`), so no job runs on for a
  supervisor nobody is watching.
* **What a job may run.** `runner` is a dotted function name that `allowlist.check_function` allows, as the function of a queued job is; it is called `runner(job, link)` in the job process and
  its return is the exit code (`None` is 0). Nothing the supervisor sends names anything else to run.

Built and tested, not yet used: the bulk lane moves to it in OPS.228.
"""

import contextlib
import gc
import importlib
import json
import os
import select
import signal
import socket
import sys
import threading
import time

from planetgen.queue import allowlist

EXIT_SUPERVISOR_GONE = 75
"""int: The exit code of a job whose supervisor went away (`SupervisorGone`): it finished what it was committing and stopped (BSD's `EX_TEMPFAIL`: try again)."""

EXIT_FAILED = 1
"""int: The exit code of a job whose runner raised."""

USE_PIDFD = os.environ.get("PLANETGEN_TEMPLATE_PIDFD", "1") != "0" and hasattr(os, "pidfd_open")
"""bool: Whether the template hands the supervisor a `pidfd` of each job (Linux). Off where the kernel has none (macOS), and when `PLANETGEN_TEMPLATE_PIDFD=0` says so, which is how the path without them is tried
on Linux: the template's report of each exit is then the way, and a template that has died leaves the supervisor to look for each job's pid (`jobhost.JobHost.poll`)."""

MESSAGE_BYTES = 1 << 20
"""int: The most one message may hold."""


class SupervisorGone(Exception):
    """The pipe to the supervisor closed: whoever watches this job is gone, so it finishes what it is committing and exits."""


def _function(name):
    """The function the dotted `name` names, which a job may name (`allowlist.check_function`)."""
    allowlist.check_function(name, "the job's runner")
    module, _, function = name.rpartition(".")
    return getattr(importlib.import_module(module), function)


def encode(message):
    return json.dumps(message, separators=(",", ":")).encode("utf-8")


def decode(raw):
    return json.loads(raw.decode("utf-8"))


class Link:
    """A job's side of its pipe to the supervisor."""

    def __init__(self, sock):
        self.sock = sock
        self.sock.setblocking(False)
        self.gone = False
        self.stop_requested = False
        self.pending = []

    def send(self, message):
        """Says `message` (a JSON-able dict) to the supervisor. Raises `SupervisorGone` when the pipe is closed."""
        if self.gone:
            raise SupervisorGone("the supervisor closed its end of the pipe")
        try:
            self.sock.send(encode(message))
        except BlockingIOError:        # a supervisor that is not reading: what a job says is progress, and the next message says it again
            return
        except (BrokenPipeError, ConnectionResetError, OSError) as exc:
            self.gone = True
            raise SupervisorGone(f"the pipe to the supervisor is closed ({exc})") from exc

    def _drain(self):
        while not self.gone:
            try:
                raw = self.sock.recv(MESSAGE_BYTES)
            except BlockingIOError:
                return
            except OSError:
                self.gone = True
                return
            if not raw:
                self.gone = True
                return
            message = decode(raw)
            if message.get("op") == "stop":
                self.stop_requested = True
            self.pending.append(message)

    def requests(self):
        """The requests the supervisor has sent since the last call (dicts, oldest first); `{"op": "stop"}` also sets `stop_requested`."""
        self._drain()
        taken, self.pending = self.pending, []
        return taken

    def check(self):
        """Raises `SupervisorGone` if the supervisor has closed its end; a safe point of a job calls this between two pieces (what was sent stays for `requests`)."""
        self._drain()
        if self.gone:
            raise SupervisorGone("the supervisor closed its end of the pipe")


def run_job(runner, job, sock):
    """What a forked job process does: runs `runner(job, link)` and returns the exit code (`SupervisorGone` is `EXIT_SUPERVISOR_GONE`, any other exception `EXIT_FAILED`, said on the pipe if it still is one)."""
    link = Link(sock)
    try:
        code = _function(runner)(job, link)
        return 0 if code is None else int(code)
    except SupervisorGone:
        return EXIT_SUPERVISOR_GONE
    except SystemExit as exc:
        return exc.code if isinstance(exc.code, int) else (0 if exc.code is None else EXIT_FAILED)
    except BaseException as exc:  # noqa: BLE001 -- a job's failure is its exit code and one line on the pipe
        with contextlib.suppress(Exception):
            link.send({"event": "error", "error": f"{type(exc).__name__}: {exc}"[:500]})
        return EXIT_FAILED


class Template:
    """The template process: warm set built, nothing open, one thread, a loop that forks a job on request and reports how each one ended."""

    def __init__(self, control, warm=(), clock=time.time):
        self.control = control
        self.warm = list(warm)
        self.clock = clock
        self.children = {}          # pid -> the id the supervisor gave the job
        self.wake_read = self.wake_write = None
        self.stopping = False

    def build(self):
        """Builds the warm set (`warm`, dotted names of functions called with no argument), freezes what it made so a fork never copies it by collecting, and checks the process is safe to fork."""
        began = self.clock()
        for name in self.warm:
            _function(name)()
        gc.collect()
        gc.freeze()
        threads = threading.active_count()
        if threads != 1:
            raise RuntimeError(f"a template must run one thread to be safe to fork; it runs {threads}")
        return self.clock() - began

    def _watch_children(self):
        """Routes a child's end to the loop at once: SIGCHLD writes a byte to a pipe the loop selects on."""
        self.wake_read, self.wake_write = os.pipe()
        os.set_blocking(self.wake_read, False)
        os.set_blocking(self.wake_write, False)
        signal.set_wakeup_fd(self.wake_write, warn_on_full_buffer=False)
        signal.signal(signal.SIGCHLD, lambda number, frame: None)

    def say(self, message, fds=()):
        socket.send_fds(self.control, [encode(message)], list(fds)) if fds else self.control.send(encode(message))

    def reap(self):
        """Reports every child that has ended (`os.waitpid`, no wait): `{"event": "exit", "id", "pid", "code", "signal"}`. Returns how many."""
        count = 0
        while self.children:
            try:
                pid, status = os.waitpid(-1, os.WNOHANG)
            except ChildProcessError:
                break
            if pid == 0:
                break
            job_id = self.children.pop(pid, None)
            if job_id is None:
                continue
            count += 1
            code = os.waitstatus_to_exitcode(status)
            self.say({"event": "exit", "id": job_id, "pid": pid, "code": code if code >= 0 else None, "signal": -code if code < 0 else None})
        return count

    def fork(self, job_id, runner, job):
        """Forks one job: its pipe is made, the child runs `runner(job, link)` and ends with its exit code, and the supervisor is sent the child's pid with its end of the pipe."""
        supervisor_end, job_end = socket.socketpair(socket.AF_UNIX, socket.SOCK_SEQPACKET)
        pid = os.fork()
        if pid == 0:
            code = EXIT_FAILED
            try:
                signal.set_wakeup_fd(-1)
                signal.signal(signal.SIGCHLD, signal.SIG_DFL)
                for descriptor in (self.control.fileno(), self.wake_read, self.wake_write, supervisor_end.fileno()):
                    with contextlib.suppress(OSError):
                        os.close(descriptor)
                code = run_job(runner, job, job_end)
            finally:
                sys.stdout.flush()
                sys.stderr.flush()
                os._exit(code)
        job_end.close()
        self.children[pid] = job_id
        pidfd = os.pidfd_open(pid) if USE_PIDFD else None       # opened here, while the child is still this process's to reap: the pid cannot be gone or another process's
        try:
            self.say({"event": "forked", "id": job_id, "pid": pid, "pidfd": pidfd is not None}, fds=[supervisor_end.fileno()] + ([pidfd] if pidfd is not None else []))
        finally:
            supervisor_end.close()
            if pidfd is not None:
                os.close(pidfd)

    def handle(self, message):
        operation = message.get("op")
        if operation == "stop":
            self.stopping = True
        elif operation == "fork":
            try:
                allowlist.check_function(message["runner"], "the job's runner")
            except allowlist.ForbiddenFunction as exc:
                self.say({"event": "error", "id": message.get("id"), "error": str(exc)})
                return
            self.fork(message["id"], message["runner"], message.get("job") or {})

    def serve(self):
        """The template's loop: builds the warm set, says it is ready, and forks jobs and reports their ends until it is told to stop or the supervisor closes the channel; then it waits for the jobs it has left and ends."""
        seconds = self.build()
        self._watch_children()
        self.say({"event": "ready", "seconds": round(seconds, 3), "pid": os.getpid()})
        open_ = True
        while open_ and not self.stopping:
            readable, _, _ = select.select([self.control, self.wake_read], [], [], 1.0)
            if self.wake_read in readable:
                with contextlib.suppress(BlockingIOError):
                    os.read(self.wake_read, 4096)
            self.reap()
            if self.control in readable:
                raw = self.control.recv(MESSAGE_BYTES)
                if not raw:
                    open_ = False
                else:
                    self.handle(decode(raw))
        while self.children:         # the jobs are not stopped by the end of the template: they end on their own, at their next safe point when the pipe is closed
            select.select([self.wake_read], [], [], 0.5)
            with contextlib.suppress(BlockingIOError):
                os.read(self.wake_read, 4096)
            try:
                self.reap()
            except OSError:
                break


def main(argv=None):
    """`python -m planetgen.queue.template FD [WARM ...]`: the template process, whose channel to the supervisor is the open descriptor FD."""
    argv = list(sys.argv[1:] if argv is None else argv)
    control = socket.socket(fileno=int(argv[0]))
    Template(control, argv[1:]).serve()
    return 0


if __name__ == "__main__":
    sys.exit(main())
