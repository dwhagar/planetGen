"""The jobs of tests/test_requires.py: each writes when it started and ended to a Redis list, so a test can read the order."""

import time

from planetgen.queue import redisqueue


def mark(label, seconds=0.0, fail=False):
    connection = redisqueue.connect()
    connection.rpush("test:requires:order", f"start {label}")
    if seconds:
        time.sleep(seconds)
    if fail:
        connection.rpush("test:requires:order", f"fail {label}")
        raise RuntimeError(f"{label} failed on purpose")
    connection.rpush("test:requires:order", f"end {label}")
    return label
