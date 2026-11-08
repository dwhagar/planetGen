# planetgen/web/lib/estimatejob.py

"""
The Generate page's size-and-time estimate, worked out in the work queue (ADM.141).

`planetgen sectors ... --estimate-only` reads the whole plan and can take minutes on a big galaxy. Run inside the web
request it outlived Apache's `WSGIDaemonProcess request-timeout` (60 s), which restarts the daemon and answers
"Truncated or oversized response headers"; and the page's live size line asked again on every change, each ask
starting the command once more. Now the request only checks the cache (`estimatecache.estimate_in_background`),
queues `work` and answers "working"; the page asks again every few seconds, and `work` stores the answer (or why
there is none) under the key the request reads.
"""

import json
import os

from planetgen.queue import api_jobs
from planetgen.web.lib import estimatecache

TIMEOUT_SECONDS = 600
"""int: How long the estimate command may run (`estimatecache.CLAIM_SECONDS` outlives it)."""

ESTIMATE_PREFIX = "ESTIMATE "
"""str: `generation.run_common.ESTIMATE_PREFIX`: the line the estimate is read from."""

TIMED_OUT_MESSAGE = "Working out the size and time took too long. Nothing was generated."


class EstimateError(Exception):
    """The estimate could not be made; the message says why (it is shown to the admin)."""


def read(done):
    """
    The estimate in the finished run of the estimate command.

    Args:
        done (dict): `api_jobs.run_command`'s answer.

    Returns:
        dict: `stats.Estimate.as_dict()` plus `what`.

    Raises:
        EstimateError: The command timed out, did not start, or printed no estimate (its own last lines are the
            message, e.g. a shell too large to generate without a limit).
    """
    if done["timed_out"]:
        raise EstimateError(TIMED_OUT_MESSAGE)
    if done["error"]:
        raise EstimateError(f"The size and time couldn't be worked out: {done['error']}")
    for line in reversed(done["stdout"].splitlines()):
        if line.startswith(ESTIMATE_PREFIX):
            try:
                return json.loads(line[len(ESTIMATE_PREFIX):])
            except ValueError:
                break
    tail = [line.strip() for line in (done["stderr"] + "\n" + done["stdout"]).splitlines() if line.strip()]
    message = " ".join(tail[-3:]) if tail else f"the generator exited with {done['returncode']}"
    raise EstimateError(f"Nothing was generated: {message}")


def work(key, versions, argv, env, cwd, timeout=TIMEOUT_SECONDS):
    """
    The queued body of an estimate (a module-level function so a worker can import it by name): runs the command,
    stores its estimate under `key` (`estimatecache.store`, which gives up the claims), or why there is none.

    Args:
        key (str): `estimatecache.cache_key`.
        versions: What the estimate is made under (`estimatecache.estimate_in_background`'s `current`).
        argv (list): The job step's command line; `--estimate-only` is added here.
        env (dict): The step's MySQL environment (`jobs.mysql_env`), laid over this process's.
        cwd (str): The repository.
        timeout (int): Seconds the command may run.

    Returns:
        dict: `{"stored": "estimate"}` or `{"stored": "error"}`.
    """
    try:
        done = api_jobs.run_command(argv + ["--estimate-only"], {**os.environ, **env}, cwd, timeout)
        estimate = read(done)
    except EstimateError as exc:
        estimatecache.store(key, versions, error=exc)
        return {"stored": "error"}
    except BaseException as exc:
        estimatecache.store(key, versions, error=f"The size and time couldn't be worked out: {exc}")
        raise
    estimatecache.store(key, versions, estimate=estimate)
    return {"stored": "estimate"}


def start(key, versions, argv, env, cwd):
    """
    Queues `work` (once: the caller holds the claim of `key`).

    Returns:
        bool: Whether it was queued; `False` with no work queue to take it (no Redis server, an update running).
    """
    try:
        api_jobs.submit(work, key, versions, list(argv), dict(env), cwd)
    except api_jobs.NoQueue:
        return False
    return True
