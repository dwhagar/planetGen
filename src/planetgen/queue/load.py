# planetgen/queue/load.py

"""
The server's load for the admin queue page (ADM.10), as the three
numbers "x / x / x" over the last 1, 5 and 15 minutes.

Linux and macOS have a load average (`os.getloadavg`): runnable
processes, decaying over each window. Windows has none, so there the
three numbers are the CPU in use, in percent, averaged over the same
windows the same way (Boss, 2026-10-01): a background thread in this
process samples the CPU every `SAMPLE_SECONDS` (`GetSystemTimes`) and
keeps one decaying average per window, so the numbers fill in over the
first minutes after the web server starts.
"""

import math
import os
import threading
import time

WINDOWS_SECONDS = (60, 300, 900)
"""tuple: The three averaging windows: 1, 5 and 15 minutes."""

SAMPLE_SECONDS = 5.0
"""float: How often the Windows sampler reads the CPU."""


def load_average():
    """
    The server's load now.

    Returns:
        dict: `values` (three floats, or `None` while the Windows sampler
            has no reading yet), `kind` ("load", or "cpu" for CPU percent)
            and `text` ("0.52 / 0.61 / 0.70", "12% / 9% / 8%", or
            "measuring" / "unknown").
    """
    if hasattr(os, "getloadavg"):
        try:
            values = [float(value) for value in os.getloadavg()]
        except OSError:
            values = None
        if values is not None:
            return {"values": values, "kind": "load", "text": " / ".join(f"{value:.2f}" for value in values)}
    if os.name != "nt":
        return {"values": None, "kind": "load", "text": "unknown"}
    values = _sampler().averages()
    text = "measuring" if values is None else " / ".join(f"{value:.0f}%" for value in values)
    return {"values": values, "kind": "cpu", "text": text}


def decay(average, sample, elapsed, window):
    """One step of a decaying average over `window` seconds, the way
    the Unix load average is kept: the older the reading, the less it
    counts."""
    keep = math.exp(-elapsed / window)
    return average * keep + sample * (1.0 - keep)


class CpuSampler:
    """
    Decaying averages of the CPU in use, in percent, one per window.

    Args:
        read (callable): Returns `(idle, total)` CPU times since boot.
        clock (callable): Seconds, monotonic.
    """

    def __init__(self, read, clock=time.monotonic):
        self._read = read
        self._clock = clock
        self._lock = threading.Lock()
        self._last = None
        self._averages = None
        self._thread = None

    def sample(self):
        """Takes one reading and folds it into the averages."""
        idle, total = self._read()
        now = self._clock()
        with self._lock:
            if self._last is not None:
                last_idle, last_total, last_time = self._last
                spent = total - last_total
                if spent > 0:
                    busy = 100.0 * max(0.0, min(1.0, 1.0 - (idle - last_idle) / spent))
                    elapsed = max(now - last_time, 1e-6)
                    if self._averages is None:
                        self._averages = [busy] * len(WINDOWS_SECONDS)
                    else:
                        self._averages = [decay(average, busy, elapsed, window)
                                          for average, window in zip(self._averages, WINDOWS_SECONDS)]
            self._last = (idle, total, now)

    def averages(self):
        """The three averages, or `None` before two readings; starts the
        sampling thread on first use."""
        self._start()
        with self._lock:
            return list(self._averages) if self._averages is not None else None

    def _start(self):
        with self._lock:
            if self._thread is not None:
                return
            self._thread = threading.Thread(target=self._run, name="cpu-sampler", daemon=True)
        self._thread.start()

    def _run(self):
        while True:
            try:
                self.sample()
            except Exception:  # noqa: BLE001 -- keep sampling; a missed reading just counts less
                pass
            time.sleep(SAMPLE_SECONDS)


def _windows_cpu_times():
    """`(idle, total)` CPU time in 100 ns ticks (`GetSystemTimes`;
    kernel time includes idle time)."""
    import ctypes
    from ctypes import wintypes

    idle, kernel, user = wintypes.FILETIME(), wintypes.FILETIME(), wintypes.FILETIME()
    if not ctypes.windll.kernel32.GetSystemTimes(ctypes.byref(idle), ctypes.byref(kernel), ctypes.byref(user)):
        raise OSError("GetSystemTimes failed")

    def ticks(filetime):
        return (filetime.dwHighDateTime << 32) | filetime.dwLowDateTime

    return ticks(idle), ticks(kernel) + ticks(user)


_SAMPLER = None
_SAMPLER_LOCK = threading.Lock()


def _sampler():
    global _SAMPLER
    with _SAMPLER_LOCK:
        if _SAMPLER is None:
            _SAMPLER = CpuSampler(_windows_cpu_times)
        return _SAMPLER
