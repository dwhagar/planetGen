"""
OPS.266: what a sector's generation tells the unit that runs it. A sector is one piece of a unit (`generate-sectors`, `generate-unplaced-sector`, `generate-system`), and its stars, planets and
moons are drawn inside that piece, so without these calls the unit's count stood still for as long as the sector took (18 to 42 s for a sector of one system on a small galaxy, minutes on a
crowded one). The system the sector is building is the second-level count's label ("System 3 of 40: bodies made") and each body made adds one to it; a body's air and moons, which come
later (GEN.1282), start the quiet clock again. Outside a unit (a test, `planetgen system` run by hand) every call does nothing, as `queue/unitbeat.py`'s own do.
"""

from planetgen.queue import unitbeat


def system(label):
    """The sector begins building a system: the second-level count starts again under `label`, counting the bodies the system gets."""
    unitbeat.step(label)


def made():
    """One body (a planet, a moon, a belt) exists: the count of the system under way goes up by one."""
    unitbeat.advance()


def settled():
    """A body's air, class or moons are done: the unit has moved on though no body was added."""
    unitbeat.moving()


class SaveStep:
    """What `store.save_sector` calls as it saves (`update(completed=, total=, description=)`, the shape of a `generation.steps` step): the unit's second-level count, "Saving the sector"."""

    def __init__(self, total):
        self.total = total
        unitbeat.step("Saving the sector", total)

    def update(self, completed=None, total=None, description=None, **_ignored):
        if total is not None:
            self.total = total
        if completed is not None:
            unitbeat.step_to(completed, self.total)

    def advance(self, amount=1):
        unitbeat.advance(amount)
