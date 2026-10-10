# planetgen/queue/load.py

"""
The server's load for the admin queue page (ADM.10), as the three
numbers "x / x / x" over the last 1, 5 and 15 minutes: the load average
(`os.getloadavg`), runnable processes decaying over each window.
"""

import os


def load_average():
    """
    The server's load now.

    Returns:
        dict: `values` (three floats, or `None` where the system has no
            load average), `kind` ("load") and `text` ("0.52 / 0.61 /
            0.70", or "unknown").
    """
    if hasattr(os, "getloadavg"):
        try:
            values = [float(value) for value in os.getloadavg()]
        except OSError:
            values = None
        if values is not None:
            return {"values": values, "kind": "load", "text": " / ".join(f"{value:.2f}" for value in values)}
    return {"values": None, "kind": "load", "text": "unknown"}
