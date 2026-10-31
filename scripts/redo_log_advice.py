#!/usr/bin/env python3
"""
Says how big MySQL's InnoDB redo log should be for this machine (PERF.222), for scripts/deploy-common.sh
(`show_redo_log_advice`, run by install.sh and update.sh) and by hand.

    redo_log_advice.py [--seconds N]

Reads how much redo the server has written (`Innodb_os_log_written`) and the redo capacity it has now, and prints the
capacity that holds `HOLD_SECONDS` of that write rate, with the line for the MySQL config and the `SET GLOBAL` that
changes it at run time on MySQL 8.0.30 and later. The write rate is the server's average since it started; with
`--seconds N` it is also measured over N seconds (run it that way while a galaxy is generating: the average of a quiet
server says little). Counts the checkpointer warnings (MY-014084) in the server's error log when it can read the log.
Never changes a setting or a file, and always exits 0: this only advises. Needs the planetGen virtual environment.
"""
import argparse
import math
import os
import sys
import time

GB = 1024 ** 3
HOLD_SECONDS = 1800
"""int: How much redo the log should hold at the measured rate: 30 minutes, a rule of thumb. The checkpointer has to
flush pages older than the log's tail before the log wraps; a log that holds only a few minutes of writes makes it
lag (the `log_checkpointer lagging` warning, MY-014084) whenever a run writes in bursts."""
MIN_BYTES = 1 * GB
MAX_BYTES = 16 * GB
"""int: The advice is held between these: a small log lags under a generation run; a very large one lengthens crash
recovery and takes disk the server may not have."""
WARNING_CODE = "MY-014084"


def recommend(bytes_per_second):
    """The redo log capacity in bytes that holds `HOLD_SECONDS` of `bytes_per_second`: a whole number of GB, between
    `MIN_BYTES` and `MAX_BYTES`."""
    wanted = max(0.0, bytes_per_second) * HOLD_SECONDS
    return int(min(MAX_BYTES, max(MIN_BYTES, math.ceil(wanted / GB) * GB)))


def human(size):
    """`size` bytes as MB or GB."""
    return f"{size / GB:.1f} GB" if size >= GB else f"{size / 2 ** 20:.0f} MB"


def capacity_of(variables):
    """`(capacity in bytes, settable)` from a dict of `SHOW GLOBAL VARIABLES` values: `innodb_redo_log_capacity` where
    the server has it (8.0.30 and later, which `SET GLOBAL` changes without a restart), else the older
    `innodb_log_file_size` times `innodb_log_files_in_group` (changed in the config file and a restart)."""
    if "innodb_redo_log_capacity" in variables:
        return int(variables["innodb_redo_log_capacity"]), True
    return int(variables.get("innodb_log_file_size", 0)) * int(variables.get("innodb_log_files_in_group", 1)), False


def advice(rate_bytes_per_second, capacity, settable, warnings=None, measured_over=""):
    """The lines to print: the measured rate, the capacity now, and what to set (nothing when it already holds
    `HOLD_SECONDS`). `warnings` is the count of checkpointer warnings in the error log, or `None` when unread."""
    wanted = recommend(rate_bytes_per_second)
    lines = [f"MySQL redo log: {human(rate_bytes_per_second)}/s written {measured_over}; capacity now {human(capacity)}."]
    if warnings:
        lines.append(f"  The error log holds {warnings} checkpointer warning(s) ({WARNING_CODE}, log_checkpointer lagging): "
                     "the log is too small for the writes.")
    if capacity >= wanted and not warnings:
        lines.append(f"  That holds {HOLD_SECONDS // 60} minutes of writes at this rate; nothing to change.")
        return lines
    target = wanted
    if warnings and target <= capacity:
        target = min(MAX_BYTES, capacity * 2)          # it lags at a size the rate says is enough: double it
    if settable:
        lines.append(f"  Recommended: innodb_redo_log_capacity = {target // GB}G, in the [mysqld] section of the MySQL config "
                     "(Ubuntu: /etc/mysql/mysql.conf.d/mysqld.cnf), or now without a restart:")
        lines.append(f"    mysql -e 'SET GLOBAL innodb_redo_log_capacity = {target};'")
    else:
        lines.append(f"  Recommended: innodb_log_file_size = {max(1, target // GB // 2)}G with innodb_log_files_in_group = 2 "
                     "in the [mysqld] section of the MySQL config, then restart MySQL (this version cannot change it live).")
    lines.append(f"  (Sized to hold {HOLD_SECONDS // 60} minutes of the measured writes, between {human(MIN_BYTES)} and "
                 f"{human(MAX_BYTES)}. Run `scripts/redo_log_advice.py --seconds 60` while a galaxy generates for a real rate.)")
    return lines


def count_warnings(path):
    """How many lines of the error log at `path` carry `WARNING_CODE`, or `None` if it cannot be read."""
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            return sum(1 for line in handle if WARNING_CODE in line)
    except OSError:
        return None


def _pairs(conn, sql):
    """`{name: value}` from a `SHOW ... WHERE` statement."""
    found = {}
    for row in conn.execute(sql).fetchall():
        values = list(row.values()) if isinstance(row, dict) else list(row)
        found[values[0]] = values[1]
    return found


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--seconds", type=float, default=0.0, help="also measure the write rate over this many seconds")
    args = parser.parse_args(argv)
    try:
        from planetgen.db import store
        conn = store.get_connection(store.MySQLConfig(), ensure_schema=False)
        try:
            status = _pairs(conn, "SHOW GLOBAL STATUS WHERE Variable_name IN ('Innodb_os_log_written', 'Uptime')")
            variables = _pairs(conn, "SHOW GLOBAL VARIABLES WHERE Variable_name IN ('innodb_redo_log_capacity', "
                                     "'innodb_log_file_size', 'innodb_log_files_in_group', 'log_error')")
            written = int(status["Innodb_os_log_written"])
            rate, over = written / max(1, int(status["Uptime"])), "on average since the server started"
            if args.seconds > 0:
                time.sleep(args.seconds)
                after = int(_pairs(conn, "SHOW GLOBAL STATUS WHERE Variable_name = 'Innodb_os_log_written'")
                            ["Innodb_os_log_written"])
                rate, over = (after - written) / args.seconds, f"over the last {args.seconds:g} seconds"
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- advice only: a server that cannot be asked is skipped
        print(f"MySQL redo log: not checked ({type(exc).__name__}: {exc})", file=sys.stderr)
        return 0
    capacity, settable = capacity_of(variables)
    log_path = variables.get("log_error")
    warnings = count_warnings(log_path) if log_path and os.path.isabs(log_path) else None
    print("\n".join(advice(rate, capacity, settable, warnings, over)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
