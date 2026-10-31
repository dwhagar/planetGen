# planetgen/physics/rogue_age.py

"""
How old a rogue planet is, and how long ago it was thrown out (GEN.1112, RP5)
=============================================================================

The rogue planets had an age drawn uniformly from 0.1 to 10 Gy (`ROGUE_PLANET_AGE_RANGE_GY`), whatever made them. A planet that was ejected formed with its star, so it is
as old as that star, and it spent the first part of its life lit by it: the age is the parent star's age (from the galaxy's star-formation history), the **ejection delay** is
how long after formation the planet was thrown out, and **the time since ejection** is the rest (`rogue-planet-review.md` row 10, `rogue-planet-drift.md` section 5). The three
numbers feed the heat law (a giant's cooling runs from its formation, its stellar heating stopped at the ejection), the tide the star raised until the ejection, and the
drift (`rogue_drift`, which needs the time since ejection) and the history sentence of the rogue's page. This module draws them; nothing calls it yet (GEN.1111 creates the
rogue rows).

**The parent's age** is drawn from an `AgeLaw` (`galaxy.age_law`, the star-formation history's rate over the population's age window, GEN.546), so the rogues' ages follow
the galaxy's star-forming past and not a flat range. The window starts no younger than the shortest ejection delay, since the planet cannot be ejected before it exists.

**The ejection delay** has the density `p(t)` proportional to `t^-DELAY_INDEX` between `EJECTION_DELAY_MIN_GY` and the parent's age: embryos and the planets of a crowded
disc are thrown out early (within the first few tens of Myr), and planet-planet scattering goes on at every later time with a roughly even share of each decade of time
(Chatterjee et al. 2008, ApJ 686:580, whose scattering times span 10^6 to 10^9 yr and beyond; recalled, not opened here). Index 1 is that logarithmic spread; 0 would be
uniform in time. When the ejection record (RP1, GEN.1119) gives the time of the loss itself the delay is read from it (`from_record`) and the index is not used.

Everything is a pure function of the numbers it is given; the random numbers come from the caller (one uniform number for the parent's age, one for the delay), so the draw is
the caller's own `draw.py` stream. Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.88.
"""

import math
from dataclasses import dataclass

EJECTION_DELAY_MIN_GY = 1.0e-3
"""float: The earliest an ejection is placed after formation, Gy (1 Myr: after the gas disc's first stage, the embryos of a crowded disc being thrown out). Source tier T8, tweak level L2."""
DELAY_INDEX = 1.0
"""float: `p(t)` proportional to `t^-DELAY_INDEX` for the ejection delay `t`; 1 is an even share of each decade of time, 0 uniform in time. Source tier T8, tweak level L3."""


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


def ejection_delay_gy(u, available_gy, minimum_gy=EJECTION_DELAY_MIN_GY, index=DELAY_INDEX):
    """
    The ejection delay, Gy, by the inverse cumulative of `p(t)` proportional to `t^-index` between `minimum_gy` and `available_gy`.

    Args:
        u (float): A uniform random number, 0 to 1.
        available_gy (float): The longest it can be: the parent's age.
        minimum_gy (float): The shortest.
        index (float): The power-law index.
    """
    _unit(u, "u")
    if not minimum_gy > 0.0 or not available_gy > 0.0:
        raise RogueAgeError("the delays must be positive")
    if available_gy <= minimum_gy:
        return available_gy
    if abs(index - 1.0) < 1e-12:
        return minimum_gy * (available_gy / minimum_gy) ** u
    power = 1.0 - index
    low, high = minimum_gy ** power, available_gy ** power
    return (low + u * (high - low)) ** (1.0 / power)


def draw_age(age_law, u_age, window=None, minimum_gy=EJECTION_DELAY_MIN_GY):
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


def rogue_age(age_law, u_age, u_delay, window=None, minimum_gy=EJECTION_DELAY_MIN_GY, index=DELAY_INDEX):
    """
    A rogue's age, ejection delay and time since ejection from two uniform random numbers.

    Args:
        age_law, window, minimum_gy: As for `draw_age`.
        u_age, u_delay (float): Uniform random numbers, 0 to 1: the parent's age and the delay.
        index (float): The delay's power-law index.

    Returns:
        RogueAge: The three ages.
    """
    age = draw_age(age_law, u_age, window, minimum_gy)
    delay = ejection_delay_gy(u_delay, age, minimum_gy, index)
    return RogueAge(age, delay, age - delay)


def from_record(parent_age_gy, ejection_delay_gy_):
    """The ages of a rogue whose ejection record (RP1) gives the parent's age and the time of the loss after formation."""
    if not parent_age_gy > 0.0 or not 0.0 <= ejection_delay_gy_ <= parent_age_gy:
        raise RogueAgeError("the parent's age must be positive and the ejection no later than it")
    return RogueAge(parent_age_gy, ejection_delay_gy_, parent_age_gy - ejection_delay_gy_)
