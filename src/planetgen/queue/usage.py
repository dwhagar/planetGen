# planetgen/queue/usage.py

"""
How much memory and processor a unit or job used over its life, as a peak and an average (OPS.155; `jobprofiles` keeps the result).

A `Usage` is fed by whatever already wakes up while the work runs: the unit's beat thread (`unitbeat.UnitBeat` calls `tick` on every tick of its thread) or the gate's
heartbeat (`memorygate.refresh`). It looks at the memory (the `resident` function the caller gives, MB, shared pages counted once as `units.resident_mb` does) and at the
processor seconds of this process and the children it has waited for, at most once per `interval` seconds, so a unit that finishes in a second costs two looks and a long one
a look every few seconds. The work loop is never touched.

* **Memory.** The peak is the largest reading; the average is the mean of the readings.
* **Processor.** In cores: processor seconds over wall seconds. The peak is the busiest interval between two looks, the average the whole life (`cpu_seconds` over `seconds`).
"""

import os
import time

from planetgen.util import log

KINDS = ("mem_peak_mb", "mem_avg_mb", "cpu_peak", "cpu_avg")

INTERVAL_SECONDS = 5.0
"""float: How often at most a `Usage` looks."""


def cpu_seconds(pid=None):
    """Processor seconds (user and system, with the children it waited for) of process `pid` (this one by default); 0.0 where the system does not say."""
    try:
        if pid is None or pid == os.getpid():
            times = os.times()
            return float(times.user + times.system + times.children_user + times.children_system)
        import psutil
        times = psutil.Process(pid).cpu_times()
        return float(times.user + times.system + getattr(times, "children_user", 0.0) + getattr(times, "children_system", 0.0))
    except Exception as exc:        # noqa: BLE001 -- a process that is gone or a system without the figure
        log.swallowed("queue", f"reading the processor time of process {pid}", exc)
        return 0.0


class Usage:
    def __init__(self, resident, clock=time.monotonic, cpu=cpu_seconds, interval=INTERVAL_SECONDS):
        self.resident, self.clock, self.cpu, self.interval = resident, clock, cpu, interval
        self.started = self.clock()
        self._cpu0 = self._cpu_last = self.cpu()
        self._t_last = None
        self._memory = []
        self._peak_cpu = 0.0
        self._looked = None

    def tick(self, force=False):
        """Looks now if `interval` seconds have passed since the last look (or `force`)."""
        now = self.clock()
        if not force and self._looked is not None and now - self._looked < self.interval:
            return
        self._looked = now
        try:
            reading = float(self.resident())
        except Exception as exc:        # noqa: BLE001
            log.swallowed("queue", "reading the resident memory of a unit", exc)
            reading = 0.0
        if reading > 0.0:
            self._memory.append(reading)
        used = self.cpu()
        since = self._t_last if self._t_last is not None else self.started
        if now - since >= 0.25:
            self._peak_cpu = max(self._peak_cpu, max(used - self._cpu_last, 0.0) / (now - since))
            self._cpu_last, self._t_last = used, now

    def summary(self):
        """`{"mem_peak_mb", "mem_avg_mb", "cpu_peak", "cpu_avg", "cpu_seconds", "seconds"}`; a figure nothing was read for is `None`. Looks once more first."""
        self.tick(force=True)
        seconds = max(self.clock() - self.started, 0.0)
        used = max(self.cpu() - self._cpu0, 0.0)
        return {"mem_peak_mb": max(self._memory) if self._memory else None, "mem_avg_mb": sum(self._memory) / len(self._memory) if self._memory else None,
                "cpu_peak": self._peak_cpu if self._peak_cpu > 0.0 else None, "cpu_avg": used / seconds if seconds >= 0.25 and used > 0.0 else None,
                "cpu_seconds": used, "seconds": seconds}
