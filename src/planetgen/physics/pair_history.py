# planetgen/physics/pair_history.py

"""
What a pair of stars went through before it looked like this (GEN.1012; docs/design/systems/overlap-free-systems.md 4.3)
=====================================================================================================================

A star that is a white dwarf, a helium core or a neutron star's parent today was a giant in the past, and a giant of the past filled a lobe that the pair's
present orbit may not have had room for. `multistar.classify_pair` knew this for a white dwarf only (`mass_loss.engulfment_radius_au`, a lookup by the
progenitor's mass); this rule is the one that holds for any pair and any member.

For a pair at `axis_au` with eccentricity `e` now, each member has a present mass `mass`, the largest it ever had `mass_max`, and the largest radius it ever had
`radius_max_au` (the red-giant or AGB tip, or its present radius for a star that never shrank):

* **The orbit of the past.** The mass the pair lost on its giant branches left slowly and evenly, so the orbit kept its angular momentum and grew:
  `a (M_a + M_b)` is constant ([R1745, R2149, R4534]; the same law as `mass_loss.expansion_factor`), so the
  separation of the past was `a_then = a (M_a + M_b) / (M_a_max + M_b_max)`, smaller than the present one, with the same eccentricity.
* **The lobe of the past.** At the closest approach of that orbit, `a_then (1 - e)`, each member's largest radius is set against its Roche lobe at its largest
  mass ([R1222], `roche`). A member that reached its lobe, or whose radius reached the closest approach itself (it swallowed its partner), went through
  mass transfer or a common envelope, and the pair is allowed only as a post-common-envelope pair, with the axis the law for that case draws, not the general
  period law.
* **A pair of members that never shrank** has the present orbit for its past one, and its past filling is the present filling `classify_pair` already measures.

The adiabatic form needs the mass to leave slowly against the orbital period: it does for a giant branch and does not for a supernova ([R4534],
section 2.4), so a neutron star or a black hole in the pair is reported `sudden` and the separation of the past is only the orbit the loss would have had.
Stand-ins are the rows GM10 and GM11 of `docs/design/research/gap-filling-plan.md` section 6.60.
"""

import math
from dataclasses import dataclass

from planetgen.physics import roche

FILL_LIMIT = 1.0
"""float: A member whose largest radius is this multiple of its Roche lobe at the closest approach of the past orbit has filled the lobe. Source tier T1 [R1222], tweak
level L0."""
SHRINK_TOLERANCE = 1.0e-9
"""float: How much larger a member's largest radius or mass must be than its present one to count as having shrunk (float noise of a record that kept the same number twice). Source tier T6, tweak
level L0."""


@dataclass(frozen=True)
class Member:
    """One member of a pair as the rule reads it: masses in solar masses (any one unit, only ratios and sums matter), radii in AU (present, then the largest it had)."""
    mass: float
    mass_max: float
    radius_au: float
    radius_max_au: float
    sudden: bool = False
    """bool: The mass left in a supernova (a neutron star or a black hole), not slowly."""


@dataclass(frozen=True)
class PairHistory:
    """The result of the rule: the pair's past orbit and whether it interacted in it."""
    past_axis_au: float
    past_pericentre_au: float
    fill_a: float
    fill_b: float
    swallowed: bool
    regime: str

    @property
    def interacted(self):
        """Whether either member filled its lobe, or swallowed its partner, in the past."""
        return self.swallowed or max(self.fill_a, self.fill_b) >= FILL_LIMIT


def past_axis_au(axis_au, member_a, member_b):
    """The separation, AU, the pair had when its members had their largest masses: the present axis over the orbit's growth (`mass_loss.expansion_factor` with the pair's mass)."""
    return axis_au * (member_a.mass + member_b.mass) / (member_a.mass_max + member_b.mass_max)


def pair_history(axis_au, eccentricity, member_a, member_b):
    """
    The `PairHistory` of a pair at `axis_au` and `eccentricity`: its orbit of the past, the filling of each member's lobe by its largest radius at that orbit's closest
    approach, whether either member's radius reached the closest approach itself, and the regime of the mass loss (`adiabatic`, or `sudden` when a member lost its mass in
    a supernova).
    """
    then = past_axis_au(axis_au, member_a, member_b)
    peri = then * (1.0 - eccentricity)
    fill_a = member_a.radius_max_au / roche.roche_lobe_radius_au(peri, member_a.mass_max, member_b.mass_max)
    fill_b = member_b.radius_max_au / roche.roche_lobe_radius_au(peri, member_b.mass_max, member_a.mass_max)
    swallowed = max(member_a.radius_max_au, member_b.radius_max_au) >= peri
    regime = "sudden" if member_a.sudden or member_b.sudden else "adiabatic"
    return PairHistory(then, peri, fill_a, fill_b, swallowed, regime)


def _shrank(member):
    """Whether the member was larger or heavier in the past than it is."""
    return member.radius_max_au > member.radius_au * (1.0 + SHRINK_TOLERANCE) or member.mass_max > member.mass * (1.0 + SHRINK_TOLERANCE)


def must_be_post_common_envelope(axis_au, eccentricity, member_a, member_b):
    """Whether a pair at this present orbit can only be a post-common-envelope pair: its members interacted in the past (`PairHistory.interacted`). A pair whose
    members never shrank is never forced: its past orbit is the present one and the present filling decides."""
    if not (_shrank(member_a) or _shrank(member_b)):
        return False
    return pair_history(axis_au, eccentricity, member_a, member_b).interacted


def widest_forced_axis_au(eccentricity, member_a, member_b, upper_au=1.0e6):
    """
    The widest present axis, AU, at which the pair is forced into a post-common-envelope pair: the present orbit grows with the mass the pair lost, so a wider
    present orbit than this had a past orbit that cleared the largest radius. Found by bisection on `pair_history` (the past filling falls monotonically with the
    axis); 0 when the pair is not forced at any axis, and `upper_au` when even that axis is.
    """
    if not must_be_post_common_envelope(upper_au * 1.0e-12, eccentricity, member_a, member_b):
        return 0.0
    if must_be_post_common_envelope(upper_au, eccentricity, member_a, member_b):
        return upper_au
    low, high = upper_au * 1.0e-12, upper_au
    for _ in range(200):
        middle = math.sqrt(low * high)
        if must_be_post_common_envelope(middle, eccentricity, member_a, member_b):
            low = middle
        else:
            high = middle
        if high / low < 1.0 + 1.0e-12:
            break
    return low
