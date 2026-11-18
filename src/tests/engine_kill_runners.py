"""The unit the real-process kill tests of the job engine run (`tests/test_engine_kill_processes.py`, OPS.236): a job forked from a template that takes its attempt on the control database as a unit does, does
its work in numbered pieces written to a file (a line when a piece starts and one when it is done, so a test can see every piece finished once and none twice), looks at its pipe to the supervisor between
two pieces, and ends its attempt and its row as `forkedjob.run` would. It is named to the allow-list by the test (`allowlist.EXTRA_ENV`)."""

import os
import signal
import time

from planetgen.db import store
from planetgen.queue import ledger


def _config(named):
    return store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"], password=named["password"], database=named["database"])


def _done_pieces(path):
    if not os.path.exists(path):
        return set()
    with open(path, encoding="ascii") as handle:
        return {int(line.split()[1]) for line in handle if line.startswith("done")}


def unit(job, link):
    """Runs the unit `job` names. `job["args"]`: `database` (the control database's parameters), `marks` (the file), `pieces`, `piece_seconds`, and the faults: `die_before_attempt` (the first process to
    run it kills itself before it takes its attempt), `kill_at` (`[attempt, piece]`: the process of that attempt kills itself at the start of that piece, as the kernel would)."""
    args = job["args"]
    marks = args["marks"]
    if args.get("die_before_attempt") and not os.path.exists(marks + ".started"):
        open(marks + ".started", "w").close()
        os.kill(os.getpid(), signal.SIGKILL)
    control = store.get_control_connection(_config(args["database"]))
    attempt = ledger.begin_attempt(control, job["row"], pid=os.getpid())
    finished = _done_pieces(marks)             # a unit that resumes starts after the pieces that were committed
    for piece in range(args["pieces"]):
        if piece in finished:
            continue
        if args.get("kill_at") == [attempt, piece]:
            os.kill(os.getpid(), signal.SIGKILL)
        with open(marks, "a", encoding="ascii") as out:
            out.write(f"start {piece} {attempt}\n")
        time.sleep(args.get("piece_seconds", 0.1))
        with open(marks, "a", encoding="ascii") as out:
            out.write(f"done {piece} {attempt}\n")
        link.check()
    ledger.end_attempt(control, job["row"], attempt, "done")
    ledger.set_job_state(control, job["row"], "done")
    control.commit()
    return 0
