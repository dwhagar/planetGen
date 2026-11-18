"""A supervisor's launcher in a process of its own, for the kill tests of the job engine (`tests/test_engine_kill_processes.py`, OPS.236): it forks the units of the ledger from a real template, turn after turn,
until it is killed. The test names the runner and the control database on the command line as JSON."""

import json
import sys
import time

from planetgen.db import store
from planetgen.queue import jobhost, launcher, memorygate


def main():
    spec = json.loads(sys.argv[1])
    named = spec["database"]
    config = store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"], password=named["password"], database=named["database"])
    control = store.get_control_connection(config)
    forker = jobhost.TemplateSet(warm_sets={"plain": []}, ready_seconds=60)
    one = launcher.Launcher(control, forker, runner=spec["runner"], read_host=lambda: memorygate.Host(16000.0, 12000.0, 8, 0.2))
    while True:
        one.turn()
        time.sleep(0.1)


if __name__ == "__main__":
    main()
