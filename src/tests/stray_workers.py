# tests/stray_workers.py

"""
Queue workers a test run left behind (TEST.393).

A test that starts real workers (`workers.ensure_workers`, `redisqueue.start_detached`) starts them in a session of their
own so that they outlive the process that started them; one that is stopped part way (a test past its time limit, a
phase the runner stopped, a pytest worker that was killed) never reaches the code that tells them to stop, and each keeps
its memory (about 300 MB with its horse) until it has been idle for `workers.IDLE_SECONDS`, which it never is while
it holds a unit of a database that has been dropped. Fifty of them filled a 16 GB machine and MySQL was killed.

Every test worker's Redis is a database of its own, 1 to 15 (`private_redis`), and a real worker is started with
`--url` naming it, so a worker whose `--url` names one of those databases is a test's, never the machine's. The
conftest ends its own xdist worker's strays when its session ends, and `scripts/test.py` ends and reports all of them
after a run.
"""

import os
import signal
import subprocess
import time
from urllib.parse import urlsplit

WORKER_COMMAND = "planetgen.cli.worker"
"""str: What every planetgen worker process has in its command line (`planetgen.queue.workers.WORKER_COMMAND`)."""

TEST_DATABASES = range(1, 16)
"""range: The Redis databases the test workers use (`private_redis.url_for_worker`); the machine's own is 0."""

TERM_WAIT_SECONDS = 10
"""int: How long a worker has to end after SIGTERM before SIGKILL."""


def listing():
    """The process list `parse` reads (`pid pgid seconds-running args`), or "" where `ps` cannot say."""
    try:
        return subprocess.run(["ps", "-axww", "-o", "pid=,pgid=,etimes=,args="], capture_output=True, text=True,
                              timeout=30, check=False).stdout
    except (OSError, subprocess.SubprocessError):
        return ""


def database_of(args):
    """The Redis database number the `--url` in a command line names, or `None` when it has none."""
    words = args.split()
    for index, word in enumerate(words[:-1]):
        if word == "--url":
            path = urlsplit(words[index + 1]).path.strip("/")
            return int(path) if path.isdigit() else 0
    return None


def parse(text):
    """
    The worker processes in `text` (`listing()`'s output) that lead a session of their own, as dicts with `pid`, `seconds`
    (how long it has run) and `database` (the Redis database its `--url` names, `None` for none). A work horse is a fork
    of its worker in the worker's session, so ending the worker's group ends it; it is not listed again.
    """
    found = []
    for line in text.splitlines():
        parts = line.split(None, 3)
        if len(parts) < 4 or not all(part.isdigit() for part in parts[:3]):
            continue
        pid, group, seconds, args = int(parts[0]), int(parts[1]), int(parts[2]), parts[3]
        if pid == group and WORKER_COMMAND in args and "pytest" not in args:
            found.append({"pid": pid, "seconds": seconds, "database": database_of(args)})
    return found


def strays(text=None, databases=None, older_than=0):
    """
    The test workers in `text` (default: the process list now): those whose `--url` names one of `databases` (default
    `TEST_DATABASES`) and that have run `older_than` seconds or more.
    """
    databases = TEST_DATABASES if databases is None else databases
    return [worker for worker in parse(listing() if text is None else text)
            if worker["database"] in databases and worker["seconds"] >= older_than]


def reap(found, killpg=os.killpg, sleep=time.sleep, exists=None, wait=TERM_WAIT_SECONDS):
    """
    Ends the workers in `found` and their horses: SIGTERM to the session, SIGKILL to what is left `wait` seconds later.

    Returns:
        int: How many workers were told to end.
    """
    exists = exists or _exists
    sent = 0
    for worker in found:
        try:
            killpg(worker["pid"], signal.SIGTERM)
            sent += 1
        except ProcessLookupError:
            pass
        except PermissionError:
            pass
    if sent:
        sleep(wait)
        for worker in found:
            if exists(worker["pid"]):
                try:
                    killpg(worker["pid"], signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    pass
    return sent


def _exists(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def reap_database(url):
    """Ends the workers whose `--url` names the database of the Redis `url`, and returns how many there were."""
    path = urlsplit(url).path.strip("/")
    if not path.isdigit() or int(path) not in TEST_DATABASES:
        return 0
    return reap(strays(databases=[int(path)]))


def reap_all():
    """Ends every test worker still running and returns the workers it found (each a dict, see `parse`)."""
    found = strays()
    reap(found)
    return found
