"""The functions the template tests run in a template and its jobs (`tests/test_template.py`); named to the allow-list by the test."""

import os
import signal
import threading
import time

from planetgen.queue import template

WARMED = []


def warm():
    """A warm set: built once in the template, held by every job forked from it."""
    WARMED.append(os.getpid())


def warm_with_a_thread():
    threading.Thread(target=time.sleep, args=(30,), daemon=True).start()


def quick(job, link):
    return job.get("code", 0)


def uses_the_warm_set(job, link):
    return 0 if WARMED and WARMED[0] != os.getpid() else 3       # built before the fork, so by the template, not by this job


def kills_itself(job, link):
    os.kill(os.getpid(), signal.SIGKILL)


def explodes(job, link):
    raise ValueError("the piece broke")


def chats(job, link):
    link.send({"hello": os.getpid()})
    for _ in range(400):
        if link.stop_requested or any(request.get("op") == "stop" for request in link.requests()):
            link.send({"stopped": True})
            return 0
        time.sleep(0.02)
    return 9


def pieces(job, link):
    """Counts out pieces, checking the pipe between two of them; each piece is a line in the job's file, a commit that is finished before the job looks at the pipe again."""
    with open(job["file"], "a", encoding="ascii") as out:
        for number in range(job["count"]):
            out.write(f"start {number}\n")
            out.flush()
            time.sleep(job.get("piece_seconds", 0.1))
            out.write(f"done {number}\n")
            out.flush()
            link.check()
    return 0


def lingers(job, link):
    """Lives until killed, checking nothing: a job that outlasts its template."""
    time.sleep(60)


def unit_in_the_tests_queue(job, link):
    """`forkedjob.run` for a test: the key prefix and queue names the test's own queue uses (`tests/inline_queue.py`) are the process's own settings in production, so a test hands them over."""
    from planetgen.queue import forkedjob, memorygate, redisqueue, runrecord
    runrecord.PREFIX, memorygate.PREFIX = job["prefix"], job["gate_prefix"]
    redisqueue.BULK_QUEUE, redisqueue.MAINTENANCE_QUEUE = job["bulk_queue"], job["maintenance_queue"]
    return forkedjob.run(job, link)
