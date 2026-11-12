# planetgen/queue/splitprofile.py

"""
What the queue engine knows of a kind of job before it cuts it (OPS.161): seconds and memory per item, from the stored profile.

**This is the stub.** The profile table is OPS.155 (the `job_profiles` table of the control database, with its writers): until it is built, `lookup` answers `None` for every
kind and the engine uses each job's suggested split, as the queue does today. The reader of the OPS.155 rows (`jobprofiles`, PR 6283) replaces the body of `lookup` and nothing else: `Reading` and the name and arguments
of `lookup` are the contract. A test (or a lane before OPS.155 lands) puts readings in `STUB`.
"""

from dataclasses import dataclass

MIN_SAMPLES = 3
"""int: The fewest samples a reading stands on before the engine trusts it (the design note's fallback chain: a profile with at least 3 samples on this machine shape)."""


@dataclass(frozen=True)
class Reading:
    """
    What the profile holds for one kind of job (`planetgen/queue/splitprofile.py`, docs/design/operations/self-tuning-jobs.md, "One stored profile per job kind"):
    `seconds_per_item` the unit's own worked seconds for one item; `mb_base` the memory a unit holds before its first item (the worker and its tables) and `mb_per_item` the
    slope of the peak memory per item (0 until the kind has samples from two size buckets); `samples` how many finished units this stands on; `trust` 1.0 for this machine's shape
    and this epoch, 0.5 for a profile measured elsewhere (the half trust of the design note).
    """

    seconds_per_item: float
    mb_base: float = 0.0
    mb_per_item: float = 0.0
    samples: int = 0
    trust: float = 1.0

    def usable(self):
        """Whether the engine may choose a cut from this: enough samples, a positive time per item and some trust."""
        return self.samples >= MIN_SAMPLES and self.seconds_per_item > 0 and self.trust > 0


STUB = {}
"""dict: Readings by job kind, for a test; empty in production until OPS.155."""


def lookup(kind):
    """The `Reading` of the units of command `kind`, or `None` when the profile holds nothing for it. OPS.155's rows are keyed `unit:<command>` (`jobprofiles.record`): the table's reader asks for `unit:` + `kind`."""
    return STUB.get(kind)
