#!/usr/bin/env python3
"""
Says what protects MySQL and Redis from the kernel's out-of-memory killer on this host (PERF.223), for
scripts/deploy-common.sh (`show_oom_advice`, run by install.sh and update.sh) and by hand.

    oom_advice.py

On 2026-10-31 a bright-star scatter run took the whole 16 GB host, the kernel killed mysqld first (3:22 AM EDT, crash
recovery on restart), then the Python workers. Two things keep that from happening again:

* `OOMScoreAdjust=-900` on the MySQL and Redis units, so the kernel picks something else first (the workers can be
  started again, a killed database cannot be asked politely to resume). This prints the drop-in for each unit that is
  loaded and not at -900 or lower.
* Optional: a memory limit on `apache2.service`. The generation workers are children of the web app (PERF.95 starts
  them from it), so they sit in Apache's cgroup by design. `MemoryHigh` there makes the kernel reclaim from that group
  first and slow it, and `MemoryMax` makes an out-of-memory kill happen inside the group (it ends a worker, never
  mysqld). The scatter stages already hold their workers to half the host's memory (PERF.220), so the limits below sit
  above that. A separate workers slice is not offered: a worker started by the web app cannot move itself into one
  without root, and the limit on Apache's unit does the same work with one file.

Prints advice only: never changes a unit or a file, and always exits 0. Linux with systemd; elsewhere it prints nothing.
"""
import os
import shutil
import subprocess
import sys

WANTED_ADJUST = -900
"""int: The `OOMScoreAdjust` the databases get (the lowest the kernel allows for a service is -1000, which never dies)."""
DATABASE_UNITS = ("mysql", "mysqld", "redis-server", "redis")
"""tuple: The unit names the databases go by on Debian, Ubuntu and Red Hat; the ones that are not loaded are skipped."""
HIGH_SHARE = 0.55
MAX_SHARE = 0.65
"""float: `MemoryHigh` and `MemoryMax` for Apache's unit as shares of the host's memory, above the 50 percent a scatter
stage's workers are held to (`queue.workers.STAGE_MEMORY_SHARE`), leaving the rest to MySQL, Redis and the kernel."""


def parse_show(text):
    """`{property: value}` from `systemctl show` output."""
    found = {}
    for line in text.splitlines():
        key, sep, value = line.partition("=")
        if sep:
            found[key.strip()] = value.strip()
    return found


def needs_adjust(properties):
    """Whether a loaded unit's `OOMScoreAdjust` is above `WANTED_ADJUST` (a unit that is not loaded, or does not say,
    is left alone)."""
    if properties.get("LoadState") != "loaded":
        return False
    try:
        return int(properties.get("OOMScoreAdjust", "0")) > WANTED_ADJUST
    except ValueError:
        return False


def adjust_commands(unit):
    """The shell lines that give `unit` the drop-in."""
    return [f"sudo mkdir -p /etc/systemd/system/{unit}.service.d",
            f"printf '[Service]\\nOOMScoreAdjust={WANTED_ADJUST}\\n' | sudo tee /etc/systemd/system/{unit}.service.d/oom.conf >/dev/null"]


def apache_limits(memory_mb):
    """`(MemoryHigh, MemoryMax)` in whole GB for a host of `memory_mb`, at least 1 and 2."""
    return max(1, int(memory_mb * HIGH_SHARE / 1024)), max(2, int(memory_mb * MAX_SHARE / 1024))


def advice(unprotected, apache_loaded, apache_limited, memory_mb):
    """The lines to print: the drop-in for each unit in `unprotected`, and the optional Apache limit when `apache_loaded`
    and not `apache_limited`. Empty when there is nothing to say."""
    lines = []
    if unprotected:
        lines.append(f"warning: {', '.join(unprotected)} can be chosen first by the kernel's out-of-memory killer "
                     f"(OOMScoreAdjust above {WANTED_ADJUST}). Protect them:")
        for unit in unprotected:
            lines += [f"  {line}" for line in adjust_commands(unit)]
        lines.append("  sudo systemctl daemon-reload && sudo systemctl restart " + " ".join(unprotected))
    if apache_loaded and not apache_limited and memory_mb:
        high, maximum = apache_limits(memory_mb)
        lines += ["optional: the generation workers run inside Apache's service, so a limit there keeps a memory spike "
                  "away from MySQL (an out-of-memory kill then ends a worker):",
                  "  sudo mkdir -p /etc/systemd/system/apache2.service.d",
                  f"  printf '[Service]\\nMemoryHigh={high}G\\nMemoryMax={maximum}G\\n' | sudo tee "
                  "/etc/systemd/system/apache2.service.d/memory.conf >/dev/null",
                  "  sudo systemctl daemon-reload && sudo systemctl restart apache2"]
    return lines


def _show(unit, *properties):
    done = subprocess.run(["systemctl", "show", f"{unit}.service", *[f"-p{name}" for name in properties]],
                          capture_output=True, text=True, timeout=10, check=False)
    return parse_show(done.stdout)


def _memory_mb():
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // 2 ** 20
    except (ValueError, OSError, AttributeError):
        return 0


def main():
    if sys.platform != "linux" or shutil.which("systemctl") is None:
        return 0
    try:
        unprotected = [unit for unit in DATABASE_UNITS if needs_adjust(_show(unit, "LoadState", "OOMScoreAdjust"))]
        apache = _show("apache2", "LoadState", "MemoryMax")
        lines = advice(unprotected, apache.get("LoadState") == "loaded", apache.get("MemoryMax", "infinity") != "infinity",
                       _memory_mb())
    except (OSError, subprocess.SubprocessError) as exc:
        print(f"OOM protection: not checked ({type(exc).__name__}: {exc})", file=sys.stderr)
        return 0
    if lines:
        print("\n".join(lines), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
