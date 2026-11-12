# planetgen/db/connection_loss.py

"""
Why MySQL dropped a connection (DB.115)
=======================================

A long job stage once ended with `InterfaceError: (0, '')`, what pymysql raises on a connection that is already closed. The rollback that hid the real error and the sky unit's retry
on a new connection are in (`store.Connection.rollback`, `run_galaxy._sky_unit`); this module names the *cause* in the log when it happens, so the next loss answers the question instead of
leaving it to guesses. Nothing here changes a result.

A unit takes `facts(conn)` when it opens its connection. When the connection is lost it takes them again on the new one and `verdict` compares the two with how long no batch had finished:

* **MySQL restarted** (or was killed and brought back, for example by the kernel's out-of-memory killer): the server's `Uptime` is lower than it was, or has grown by less than the time that passed.
* **An idle timeout**: no batch had finished for `wait_timeout` seconds or more, so the server closed a connection it saw idle.
* **The server aborted the client** (`Aborted_clients` went up): a read or write timeout (`net_read_timeout`, `net_write_timeout`) or a packet over `max_allowed_packet`.
* **Something else closed the socket**: the same server, a connection idle for less than its timeout and no abort counted (a firewall, a proxy or a `KILL`).

The status is read the way `stats.server_info` reads it: an account that cannot read it gives `None` and the verdict says it cannot tell.
"""

from planetgen.util import log

_STATUS = ("Uptime", "Aborted_clients", "Aborted_connects")
"""tuple: The `SHOW GLOBAL STATUS` counters the verdict reads."""

RESTART_SLACK_S = 5
"""int: How much shorter than the time that passed an uptime may fall before the server counts as restarted (clocks and the time the status takes to read)."""


def facts(conn):
    """The server's `Uptime` (s), `Aborted_clients`, `Aborted_connects` and `wait_timeout`, `net_read_timeout`, `net_write_timeout` (s) and `max_allowed_packet` (bytes) as a dict; a value the account cannot read is `None`."""
    out = dict.fromkeys(("uptime", "aborted_clients", "aborted_connects", "wait_timeout", "net_read_timeout", "net_write_timeout", "max_allowed_packet"))
    try:
        rows = conn.execute("SHOW GLOBAL STATUS WHERE Variable_name IN ('Uptime', 'Aborted_clients', 'Aborted_connects')").fetchall()
        status = {row["Variable_name"]: int(row["Value"]) for row in rows}
        out.update(uptime=status.get("Uptime"), aborted_clients=status.get("Aborted_clients"), aborted_connects=status.get("Aborted_connects"))
    except Exception:
        pass
    try:
        row = conn.execute("SELECT @@global.wait_timeout AS wait_timeout, @@global.net_read_timeout AS net_read_timeout,"
                           " @@global.net_write_timeout AS net_write_timeout, @@global.max_allowed_packet AS max_allowed_packet").fetchone()
        out.update({key: int(row[key]) for key in ("wait_timeout", "net_read_timeout", "net_write_timeout", "max_allowed_packet")})
    except Exception:
        pass
    return out


def verdict(before, after, idle_s, passed_s):
    """
    One sentence for the log: what the two `facts` say dropped the connection. `idle_s` is the time since the unit last finished a batch, `passed_s` the time since `before` was taken.
    """
    if before.get("uptime") is None or after.get("uptime") is None:
        return "cannot tell why: this MySQL account cannot read the server status"
    if after["uptime"] < before["uptime"] or after["uptime"] + RESTART_SLACK_S < before["uptime"] + passed_s:
        return (f"MySQL was restarted: its uptime is {after['uptime']:,} s now and was {before['uptime']:,} s {passed_s:,.0f} s ago "
                "(look in `journalctl -u mysql` and `dmesg | grep -i -E 'oom|killed process'` for the time)")
    timeout = after.get("wait_timeout")
    if timeout is not None and idle_s >= timeout:
        return f"idle timeout: no batch had finished for {idle_s:,.0f} s and the server's wait_timeout is {timeout:,} s"
    grown = (after.get("aborted_clients") or 0) - (before.get("aborted_clients") or 0)
    if grown > 0:
        return (f"the server aborted the client ({grown:,} aborted clients since the unit began): a read or write timeout (net_read_timeout {after.get('net_read_timeout')} s, "
                f"net_write_timeout {after.get('net_write_timeout')} s) or a packet over max_allowed_packet ({after.get('max_allowed_packet')} bytes)")
    return (f"the same server (uptime {after['uptime']:,} s) and a connection idle for {idle_s:,.0f} s under a wait_timeout of {timeout} s: "
            "something between the job and MySQL closed the socket (a KILL, a firewall or a proxy)")


def explain(conn, before, idle_s, passed_s):
    """Logs `verdict` for a connection that was lost, from a new connection `conn`; never raises."""
    try:
        log.warning(f"MySQL connection lost: {verdict(before, facts(conn), idle_s, passed_s)}.")
    except Exception as exc:
        log.warning(f"MySQL connection lost: the cause could not be read ({exc!r}).")


def describe(after):
    """What `facts` read after a loss, in one sentence, for a place that has no earlier `facts` to compare with (OPS.171): a small uptime says the server restarted."""
    if after.get("uptime") is None:
        return "cannot tell why: this MySQL account cannot read the server status"
    return (f"the server's uptime is {after['uptime']:,} s (a small figure means it restarted), its wait_timeout {after.get('wait_timeout')} s, "
            f"{after.get('aborted_clients')} clients aborted so far and max_allowed_packet {after.get('max_allowed_packet')} bytes")


def explain_lost(config, exc):
    """Logs, as a warning, that a statement lost its connection (`exc`) and what a new connection to `config` says about the server; never raises (OPS.171)."""
    try:
        from planetgen.db import store
        conn = store.get_connection(config)
        try:
            what = describe(facts(conn))
        finally:
            conn.close()
        log.warning(f"MySQL connection lost ({exc!r}): {what}.")
    except Exception as again:
        log.warning(f"MySQL connection lost ({exc!r}): the cause could not be read ({again!r}).")
