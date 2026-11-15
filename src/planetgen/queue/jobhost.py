# planetgen/queue/jobhost.py

"""
The supervisor's side of a warm template (OPS.226, section 3.3 of `docs/design/operations/job-engine.md`): starts the template process, asks it to fork a job, holds the job's pipe and a
kernel handle on it, and tells the supervisor how each job ended. `template` is the other side.

* **`JobHost.start`** runs `python -m planetgen.queue.template` with one end of a `SOCK_SEQPACKET` pair and waits for its `ready` (the warm set is built by then).
* **`JobHost.fork(runner, job)`** asks for one fork and returns the `Job`: its id, its pid, the supervisor's end of its pipe (`Job.send`, `Job.messages`) and a `pidfd` on it where the kernel has
  them (Linux; elsewhere the template's report is the whole way, with `psutil` as the check).
* **`JobHost.poll(timeout)`** returns the `Exit` of every job that ended since the last call (its pipe stays open for what it said last until `Job.close`). The template's report carries the exit code or the signal. Where the template has died the
  kernel handle says the job ended and the `Exit` has neither (`source` is `pidfd`): the job's status was the template's to read. A job is reported once.
* **A job whose pipe is closed finishes and exits** (`template.EXIT_SUPERVISOR_GONE`): closing the `Job` (or ending this process) is how a supervisor lets go of it.

Built and tested, not yet used: the bulk lane moves to it in OPS.228.
"""

import dataclasses
import itertools
import os
import select
import signal
import socket
import subprocess
import sys
import time

from planetgen.queue import template
from planetgen.util.format import format_duration_seconds

READY_SECONDS = 600.0
"""float: How long `start` waits for the template to build its warm set (the bright-star tables take about a minute on a cold cache)."""

FORK_SECONDS = 30.0
"""float: How long `fork` waits for the template to answer."""


class HostError(Exception):
    """The template did not do what it was asked: it never got ready, it died, or it refused a job."""


@dataclasses.dataclass
class Exit:
    """How a job ended. `code` is its exit code, `signal` the signal that killed it; both are `None` when only the kernel's handle saw it end (`source` is `pidfd`)."""

    id: int
    pid: int
    code: int = None
    signal: int = None
    source: str = "template"


def _pid_alive(pid):
    """Whether a process with this pid still runs (a zombie does not count): the way to know a job ended where the kernel has no handle to give."""
    try:
        import psutil
        return psutil.pid_exists(pid) and psutil.Process(pid).status() != psutil.STATUS_ZOMBIE
    except ImportError:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True
    except psutil.NoSuchProcess:
        return False


class Job:
    """A running job as the supervisor holds it: its id, its pid, its pipe and a kernel handle on it."""

    def __init__(self, job_id, pid, sock, pidfd=None):
        self.id, self.pid, self.sock, self.pidfd = job_id, pid, sock, pidfd
        self.sock.setblocking(False)
        self.ended = False

    def send(self, message):
        """Asks the job something (`{"op": "stop"}` is the request to stop between two pieces)."""
        self.sock.send(template.encode(message))

    def messages(self):
        """What the job has said since the last call (dicts, oldest first)."""
        said = []
        while True:
            try:
                raw = self.sock.recv(template.MESSAGE_BYTES)
            except (BlockingIOError, OSError):
                return said
            if not raw:
                return said
            said.append(template.decode(raw))

    def kill(self, number=signal.SIGKILL):
        """Sends the job a signal (the last resort when it does not honour a request; `ProcessLookupError` is a job that is already gone)."""
        os.kill(self.pid, number)

    def release_handle(self):
        """Lets go of the kernel handle (a job that has ended needs none); the pipe stays, so what the job said last can be read."""
        if self.pidfd is not None:
            try:
                os.close(self.pidfd)
            except OSError:
                pass
            self.pidfd = None

    def close(self):
        """Lets go of the job: its pipe closes, which tells it to finish and exit, and the kernel handle is released."""
        for name in ("sock", "pidfd"):
            held = getattr(self, name)
            if held is None:
                continue
            try:
                held.close() if name == "sock" else os.close(held)
            except OSError:
                pass
            setattr(self, name, None)


class JobHost:
    """One template process and the jobs forked from it."""

    def __init__(self, python=None, warm=()):
        self.python, self.warm = python or sys.executable, list(warm)
        self.control = self.process = None
        self.jobs = {}              # id -> Job, until its end is reported
        self.template_pid = None
        self.template_seconds = None
        self._ids = itertools.count(1)
        self._buffer = []           # (message, fds) read ahead while waiting for another
        self._dead = False

    def start(self, ready_seconds=READY_SECONDS):
        """Starts the template and waits until its warm set is built. Returns the seconds that took."""
        mine, theirs = socket.socketpair(socket.AF_UNIX, socket.SOCK_SEQPACKET)
        self.process = subprocess.Popen([self.python, "-m", "planetgen.queue.template", str(theirs.fileno()), *self.warm], pass_fds=(theirs.fileno(),))
        theirs.close()
        self.control = mine
        message, _fds = self._expect("ready", ready_seconds)
        self.template_pid, self.template_seconds = message["pid"], message["seconds"]
        return message["seconds"]

    def _read(self, timeout):
        """One message from the template with the descriptors that came with it, or `None` after `timeout` seconds; the channel closing makes the template dead (`HostError` only when something was waiting on it)."""
        if self._buffer:
            return self._buffer.pop(0)
        if self._dead:
            return None
        readable, _, _ = select.select([self.control], [], [], max(timeout, 0.0))
        if not readable:
            return None
        raw, fds, _flags, _address = socket.recv_fds(self.control, template.MESSAGE_BYTES, 4)
        if not raw:
            self._dead = True
            return None
        return template.decode(raw), fds

    def _expect(self, event, seconds, job_id=None):
        """The next `event` message (about `job_id` when given); the others read meanwhile are kept for `poll`."""
        deadline = time.monotonic() + seconds
        kept = []
        try:
            while True:
                found = self._read(deadline - time.monotonic())
                if found is None:
                    if self._dead:
                        raise HostError("the template ended before it answered")
                    raise HostError(f"the template did not send {event!r} in {format_duration_seconds(seconds)}")
                message, fds = found
                if message.get("event") == event and (job_id is None or message.get("id") == job_id):
                    return message, fds
                if message.get("event") == "error" and message.get("id") == job_id:
                    raise HostError(message.get("error", "the template refused the job"))
                kept.append(found)
        finally:
            self._buffer = kept + self._buffer

    def fork(self, runner, job=None, seconds=FORK_SECONDS):
        """Asks the template to fork one job that runs `runner(job, link)` (a dotted name the allow-list permits). Returns the `Job`."""
        job_id = next(self._ids)
        self.control.send(template.encode({"op": "fork", "id": job_id, "runner": runner, "job": job or {}}))
        message, fds = self._expect("forked", seconds, job_id)
        sock = socket.socket(fileno=fds[0])
        pidfd = fds[1] if len(fds) > 1 else None
        held = Job(job_id, message["pid"], sock, pidfd)
        self.jobs[job_id] = held
        return held

    def poll(self, timeout=0.0):
        """The `Exit` of every job that ended since the last call, waiting up to `timeout` seconds for the first. A job is reported once."""
        ended = []
        deadline = time.monotonic() + timeout
        while True:
            while True:
                found = self._read(0.0)
                if found is None:
                    break
                message, fds = found
                for descriptor in fds:
                    os.close(descriptor)
                if message.get("event") == "exit" and message["id"] in self.jobs:
                    held = self.jobs.pop(message["id"])
                    held.ended = True
                    ended.append(Exit(message["id"], message["pid"], message.get("code"), message.get("signal")))
                    held.release_handle()       # its pipe stays open: what the job said last is read, and the caller closes it
            if self._dead:
                for held in list(self.jobs.values()):       # only the kernel can say: the template that would have read the status is gone
                    if held.pidfd is not None:
                        gone, source = bool(select.select([held.pidfd], [], [], 0.0)[0]), "pidfd"
                    else:                                    # no kernel handle (macOS): the pid is looked for, as `psutil` does
                        gone, source = not _pid_alive(held.pid), "pid"
                    if gone:
                        self.jobs.pop(held.id)
                        ended.append(Exit(held.id, held.pid, source=source))
                        held.release_handle()
            if ended or time.monotonic() >= deadline:
                return ended
            watched = [] if self._dead else [self.control]
            watched += [held.pidfd for held in self.jobs.values() if held.pidfd is not None and self._dead]
            if not watched and not (self._dead and self.jobs):      # nothing can arrive: the template is gone and no job is left to look for
                return ended
            select.select(watched, [], [], min(0.05, max(deadline - time.monotonic(), 0.0)))

    def stop(self, wait=10.0):
        """Tells the template to stop forking and lets go of every job (their pipes close, so each finishes and exits); waits for the template to end."""
        try:
            self.control.send(template.encode({"op": "stop"}))
        except OSError:
            pass
        for held in list(self.jobs.values()):
            held.close()
        if self.process is not None:
            try:
                self.process.wait(timeout=wait)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
        if self.control is not None:
            self.control.close()


WARM_SETS = {
    "tables": ("planetgen.queue.units.warm_star_tables",),
    "plain": (),
}
"""dict: The warm set of each template of a `TemplateSet`: the lane whose units draw stars (`units.TABLE_KINDS`) has the sampling tables built once, every other lane has nothing to warm."""


def template_key(kind):
    """Which template forks a unit of command `kind`: `tables` for the commands that draw from the sampling tables (`units.TABLE_KINDS`, the rule `units.warm_unit` has today), else `plain`."""
    from planetgen.queue import units
    return "tables" if kind in units.TABLE_KINDS else "plain"


class TemplateSet:
    """The supervisor's templates: one per warm set, each started when the first job that wants it is asked for (a host with no star-drawing work never builds the tables) and stopped together."""

    def __init__(self, python=None, warm_sets=None, ready_seconds=READY_SECONDS):
        self.python, self.warm_sets, self.ready_seconds = python, dict(WARM_SETS if warm_sets is None else warm_sets), ready_seconds
        self.hosts = {}

    def host(self, key):
        """The running `JobHost` of template `key`, started now if it is not (the warm set is built before this returns)."""
        if key not in self.hosts or self.hosts[key]._dead:
            made = JobHost(python=self.python, warm=self.warm_sets[key])
            made.start(self.ready_seconds)
            self.hosts[key] = made
        return self.hosts[key]

    def fork(self, kind, runner, job=None):
        """Forks the unit of `kind` from its template. Returns `(key, Job)`."""
        key = template_key(kind)
        return key, self.host(key).fork(runner, job)

    def poll(self, timeout=0.0):
        """`(key, Exit)` of every job that ended since the last call, over every template."""
        found = []
        for key, host in self.hosts.items():
            found += [(key, exit_) for exit_ in host.poll(timeout if len(self.hosts) == 1 else 0.0)]
        return found

    def stop(self):
        for host in self.hosts.values():
            host.stop()
        self.hosts = {}
