# planetgen/physics/rogue_age.py

"""
How old a rogue planet is, and how long ago it was thrown out (GEN.1112, RP5)
=============================================================================

The rogue planets had an age drawn uniformly from 0.1 to 10 Gy (`ROGUE_PLANET_AGE_RANGE_GY`), whatever made them. A planet that was ejected formed with its star, so it is
as old as that star, and it spent the first part of its life lit by it: the age is the parent star's age (from the galaxy's star-formation history), the **ejection delay** is
how long after formation the planet was thrown out, and **the time since ejection** is the rest (`rogue-planet-review.md` row 10, `rogue-planet-drift.md` section 5). The three
numbers feed the heat law (a giant's cooling runs from its formation, its stellar heating stopped at the ejection), the tide the star raised until the ejection, and the
drift (`rogue_drift`, which needs the time since ejection) and the history sentence of the rogue's page. This module draws the age and takes the delay from the record.

**The parent's age** is drawn from an `AgeLaw` (`galaxy.age_law`, the star-formation history's rate over the population's age window, GEN.546), so the rogues' ages follow
the galaxy's star-forming past and not a flat range. The window starts no younger than the shortest ejection delay, since the planet cannot be ejected before it exists.

**The ejection delay** is the time of the loss in the planet architecture: `Loss.time_yr`, the disk lifetime (`scatter.DISK_LIFETIME_YR`, 3 Myr), at which the architecture's
scatter step throws a body out. A donor's rogue reads it from the record (`from_record`); a top-up rogue, which has no donor and no record, takes the same time (`EJECTION_TIME_GY`).
The delay was a power-law draw before (index 1, an even share of each decade from 1 Myr to the parent's age, [R0747]:580 recalled) and that stand-in (RG2) is gone.

Everything is a pure function of the numbers it is given; the random number comes from the caller (one uniform number for the parent's age), so the draw is
the caller's own `draw.py` stream. The loss time is a stand-in of `scatter` (`docs/design/research/gap-filling-plan.md` section 6.127).
"""

import math
from dataclasses import dataclass

from planetgen.physics import scatter

EJECTION_TIME_GY = scatter.DISK_LIFETIME_YR / 1.0e9
"""float: When the planet architecture throws a body out, Gy after formation (`Loss.time_yr`): the disk lifetime."""


class RogueAgeError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class RogueAge:
    """
    A rogue planet's ages, Gy.

    Attributes:
        age_gy (float): Its age: the parent star's.
        delay_gy (float): How long after formation it was ejected; the time it was lit by its star.
        since_ejection_gy (float): The time since the ejection, `age_gy - delay_gy`; the time to drift (`rogue_drift`) and to cool in the dark.
    """

    age_gy: float
    delay_gy: float
    since_ejection_gy: float

    @property
    def lit_share(self):
        """The share of its life the planet spent bound to its star."""
        return self.delay_gy / self.age_gy


class _Fixed:
    """A source of random numbers that gives one number (what `AgeLaw.sample` asks for), so the draw is the caller's."""

    def __init__(self, u):
        self.u = u

    def random(self):
        return self.u

    def uniform(self, low, high):
        return low + self.u * (high - low)


def _unit(u, name):
    if not 0.0 <= u <= 1.0 or not math.isfinite(u):
        raise RogueAgeError(f"{name} must be a number from 0 to 1")


def draw_age(age_law, u_age, window=None, minimum_gy=EJECTION_TIME_GY):
    """
    The parent's age, Gy, drawn from the star-formation history's `AgeLaw` over its window (or part of it), not younger than `minimum_gy`.

    Args:
        age_law (AgeLaw): The population's age law (`galaxy.age_law.AgeLaw`, or `AgeLaw.uniform((low, high))`).
        u_age (float): A uniform random number, 0 to 1.
        window (tuple, optional): `(low, high)` Gy to draw within; default the law's whole window.
        minimum_gy (float): The youngest allowed.
    """
    _unit(u_age, "u_age")
    low, high = window if window is not None else age_law.window
    low = max(low, minimum_gy)
    if not high > low:
        raise RogueAgeError("the age window is empty once an ejection delay is allowed for")
    return age_law.sample(low, high, _Fixed(u_age))


def rogue_age(age_law, u_age, window=None, delay_gy=EJECTION_TIME_GY):
    """
    A rogue's age, ejection delay and time since ejection from one uniform random number.

    Args:
        age_law, window: As for `draw_age`.
        u_age (float): A uniform random number, 0 to 1: the parent's age.
        delay_gy (float): The ejection delay, the architecture's loss time by default; the youngest parent is no younger than it.

    Returns:
        RogueAge: The three ages.
    """
    age = draw_age(age_law, u_age, window, delay_gy)
    return from_record(age, delay_gy)


def from_record(parent_age_gy, ejection_delay_gy_):
    """The ages of a rogue whose ejection record (RP1) gives the parent's age and the time of the loss after formation."""
    if not parent_age_gy > 0.0 or not 0.0 <= ejection_delay_gy_ <= parent_age_gy:
        raise RogueAgeError("the parent's age must be positive and the ejection no later than it")
    return RogueAge(parent_age_gy, ejection_delay_gy_, parent_age_gy - ejection_delay_gy_)
