# planetgen/physics/hosting_gate.py

"""
The hosting gate: whether a body can have satellites (GEN.579, docs/design/unified-system-barycentre.md section 4.4)
====================================================================================================================

Every body other than the central object can be the central object of a system of its own, down to small bodies that get no
satellites. The recursion ends because this gate says no, not because of a depth counter. `can_host(host)` runs the three
parts of the design in order and names the first that fails:

1. **A satellite must fit.** The smallest satellite is `FLOOR_DIAMETER_M` (10 m) across and may not outweigh its host by more
   than `MASS_RATIO_LIMIT` (a tenth, the limit `physics.planets.moon_size_limits` already uses; equal-mass pairs are binary
   nodes, not satellites).
2. **A stable orbit must exist.** Outside the host's radius plus the smallest satellite's (and outside the fluid Roche limit
   (`relations.roche_limit`) when the host is fluid), inside the Domingos limit (`clearance_ledger.moon_hill_fraction`, 0.4895 for a circular orbit) of the host's Hill radius at its
   pericentre, and, when the smallest satellite is comparable in mass to its host (`COMPARABLE_MASS_RATIO`), passing the
   hierarchical stability screen of `physics.stability`.
3. **A chance that one exists at all.** Below `SIZE_CUT_M` (1 km across, Boss's card S5, a tunable stand-in) no body has satellites.
   Above it, the share of bodies of that population and size that have one comes from the gated table `satellite-share`
   (GEN.472 format), whose rows carry their source state. A moon has no row: it passes only through the tidal-survival gate
   (`tidal_survival_radius_m`, Kollmeier and Raymond 2019, equation 1) *and* a row that does not exist yet (the odds start at
   zero, Boss's card S6), so no moon draws a satellite until a sourced row is added.

A population is named for what the body is, never for a place in the Solar System: `inner-small-body` (rubble piles on orbits
that reach the inner system), `belt-small-body` (a belt between the snow line and the inner system), `trans-neptunian-cold`,
`trans-neptunian-hot`, `moon`. The caller classes the body; this module does not guess.

The rows are in `satellite-share`; their numbers and sources are in docs/design/composition/small-body-size-laws.md section 2. All are
partly verified or stand-ins. A change to a row changes draws, so it needs a generator epoch when a draw reads it. The moon draw reads this module
(GEN.761): `planets.generate_moons` asks `can_host`, and `Planet._draw_moons` takes a class D body's chance of moons from
`satellite_share` (`planets.hosting_population` names the body).
"""

import math
from dataclasses import dataclass

from planetgen.generation import possibility_table as pt
from planetgen.physics import clearance_ledger, constants, relations, stability, unit_registry as units

_yr_to_s = units.converter("yr", "s")

FLOOR_DIAMETER_M = 10.0
"""float: The smallest satellite's diameter (Boss's card S5: the 10 m floor)."""

MASS_RATIO_LIMIT = 0.1
"""float: A satellite's mass over its host's, at most (`physics.planets.moon_size_limits`); STAND-IN, tunable."""

SIZE_CUT_M = 1000.0
"""float: Below this diameter a body has no satellites (Boss's card S5). STAND-IN: GEN.580's note (section 2) shows
satellites of near-Earth rubble piles down to about 120 m; the cut stays until the owner of this item moves it."""

COMPARABLE_MASS_RATIO = 0.01
"""float: A smallest satellite this heavy compared with its host is checked with the hierarchical stability screen.
STAND-IN."""

SATELLITE_DENSITY_KG_M3 = 2000.0
"""float: The density assumed for the smallest satellite (the submoon density of Kollmeier and Raymond 2019); STAND-IN."""

# Tidal survival of a submoon (Kollmeier and Raymond 2019, MNRAS, arXiv:1810.03304, equation 1)
SUBMOON_HILL_FRACTION = clearance_ledger.moon_hill_fraction(0.0)[0]
"""float: A submoon's orbit is stable out to this fraction of its moon's Hill sphere (Domingos et al. 2006, via Kollmeier and
Raymond 2019); the circular-orbit value of `clearance_ledger.moon_hill_fraction`, 0.4895."""

SUBMOON_LOVE_NUMBER = 0.25
"""float: The moon's tidal Love number k2 (Moore and Schubert 2000, via Kollmeier and Raymond 2019, their default)."""

SUBMOON_TIDAL_Q = 100.0
"""float: The moon's tidal quality factor (Lainey 2016, via Kollmeier and Raymond 2019, their default)."""

SUBMOON_AGE_YR = 4.6e9
"""float: The time a submoon has to survive (their default: the age of the Solar System)."""

SUBMOON_DEFAULT_RADIUS_M = 10_000.0
"""float: The submoon size their figure 1 draws as the default."""

SATELLITE_SHARE_CLASSES = ("inner-small-body", "belt-small-body", "trans-neptunian-cold", "trans-neptunian-hot")
"""tuple: The populations the `satellite-share` table has rows for."""

PHOTOMETRIC_BINARY_MAX_DIAMETER_M = 13_000.0
"""float: The largest primary the lightcurve method finds a binary around, (939) Isberga, 13.4 +- 1.3 km (Margot et al.
2015); above it the share is that of a large body (adaptive optics and imaging)."""

LARGE_BODY_MIN_DIAMETER_M = 100_000.0
"""float: The size of the 'large' main-belt surveys (about 300 asteroids over 100 km, Margot et al. 2015)."""

TRANS_NEPTUNIAN_MIN_DIAMETER_M = 100_000.0
"""float: The smallest trans-Neptunian body the binary frequency was measured for (HST search at wide separations)."""


class HostingError(ValueError):
    """A host that cannot be judged: a missing or non-physical number."""


@dataclass(frozen=True)
class Host:
    """
    A body that might host satellites.

    Attributes:
        population (str): What it is: `inner-small-body`, `belt-small-body`, `trans-neptunian-cold`, `trans-neptunian-hot`
            or `moon` (the classes of `satellite-share`), or any other string for a body this table does not cover.
        diameter_m (float): Its diameter.
        mass_kg (float): Its mass.
        parent_mass_kg (float or None): The mass of what it orbits; `None` for a free body (no Hill limit).
        orbit_radius_m (float or None): Its semi-major axis about the parent.
        eccentricity (float): Its orbital eccentricity about the parent.
        fluid (bool): Whether it is fluid (the fluid Roche limit applies).
    """

    population: str
    diameter_m: float
    mass_kg: float
    parent_mass_kg: float | None = None
    orbit_radius_m: float | None = None
    eccentricity: float = 0.0
    fluid: bool = False

    def __post_init__(self):
        if not (self.diameter_m > 0 and self.mass_kg > 0):
            raise HostingError("a host needs a positive diameter and mass")
        if (self.parent_mass_kg is None) != (self.orbit_radius_m is None):
            raise HostingError("the parent's mass and the orbit's radius are given together or not at all")
        if self.parent_mass_kg is not None and not (self.parent_mass_kg > 0 and self.orbit_radius_m > 0):
            raise HostingError("the parent's mass and the orbit's radius must be positive")
        if not 0.0 <= self.eccentricity < 1.0:
            raise HostingError("a bound host has an eccentricity from 0 up to but not including 1")

    @property
    def density_kg_m3(self):
        return self.mass_kg / (math.pi / 6.0 * self.diameter_m ** 3)


@dataclass(frozen=True)
class Verdict:
    """
    The result of `can_host`.

    Attributes:
        allowed (bool): Whether the body can host satellites.
        part (int): The first part that failed (1, 2 or 3), or 0 when it passed.
        reason (str): Why, in words.
    """

    allowed: bool
    part: int
    reason: str


def smallest_satellite_mass_kg():
    """The mass of a floor-sized satellite at `SATELLITE_DENSITY_KG_M3`."""
    return SATELLITE_DENSITY_KG_M3 * math.pi / 6.0 * FLOOR_DIAMETER_M ** 3


def hill_radius_m(host):
    """The host's Hill radius at its pericentre, or `None` for a free body."""
    if host.parent_mass_kg is None:
        return None
    return relations.hill_radius(host.orbit_radius_m * (1.0 - host.eccentricity), host.mass_kg, host.parent_mass_kg)


def stable_zone_m(host):
    """
    The orbits a satellite can hold about `host`: `(inner_m, outer_m)` from the body's centre, `outer_m` `None` for a free
    body. The inner edge is the host's radius plus the smallest satellite's, or the fluid Roche limit when that is farther;
    the outer edge is the Domingos fraction of the Hill radius (`clearance_ledger.moon_hill_fraction`, 0.4895 for a circular orbit).
    """
    inner = (host.diameter_m + FLOOR_DIAMETER_M) / 2.0
    if host.fluid:
        inner = max(inner, relations.roche_limit(host.diameter_m / 2.0, host.density_kg_m3, SATELLITE_DENSITY_KG_M3))
    hill = hill_radius_m(host)
    return inner, None if hill is None else clearance_ledger.moon_hill_fraction(0.0)[0] * hill


def satellite_fits(host):
    """Part 1: a floor-sized satellite does not outweigh the host past `MASS_RATIO_LIMIT`, and the host is bigger than it."""
    mass = smallest_satellite_mass_kg()
    return mass <= MASS_RATIO_LIMIT * host.mass_kg and FLOOR_DIAMETER_M < host.diameter_m


def stable_orbit_exists(host):
    """Part 2: the stable zone is not empty and, for a comparable satellite, the triple passes the stability screen."""
    inner, outer = stable_zone_m(host)
    if outer is not None and not inner < outer:
        return False
    sat_mass = smallest_satellite_mass_kg()
    if host.parent_mass_kg is not None and sat_mass / host.mass_kg >= COMPARABLE_MASS_RATIO:
        return stability.is_triple_stable(host.mass_kg, sat_mass, host.parent_mass_kg, inner, host.orbit_radius_m,
                                          0.0, host.eccentricity, 0.0)
    return True


def tidal_survival_radius_m(moon_density_kg_m3, planet_mass_kg, moon_orbit_radius_m,
                            submoon_mass_kg=None, love_number=SUBMOON_LOVE_NUMBER, tidal_q=SUBMOON_TIDAL_Q,
                            age_yr=SUBMOON_AGE_YR):
    """
    The smallest moon radius that can hold a submoon of `submoon_mass_kg` for `age_yr`, metres (Kollmeier and Raymond 2019,
    equation 1, adapted from Barnes and O'Brien 2002):

        R_moon >= [ 39 M_sub k2 T sqrt(G) / (2 (4 pi rho_moon)^(8/3) Q) * (3 M_p / (f a_moon)^3)^(13/6) ]^(1/3)

    with `f = SUBMOON_HILL_FRACTION`. The default submoon is a 10 km radius ball at `SATELLITE_DENSITY_KG_M3`.
    """
    if submoon_mass_kg is None:
        submoon_mass_kg = SATELLITE_DENSITY_KG_M3 * 4.0 / 3.0 * math.pi * SUBMOON_DEFAULT_RADIUS_M ** 3
    if not (moon_density_kg_m3 > 0 and planet_mass_kg > 0 and moon_orbit_radius_m > 0 and submoon_mass_kg > 0
            and love_number > 0 and tidal_q > 0 and age_yr > 0):
        raise HostingError("the tidal gate needs positive numbers")
    seconds = _yr_to_s(age_yr)
    density_term = 3.0 * planet_mass_kg / (SUBMOON_HILL_FRACTION * moon_orbit_radius_m) ** 3
    bracket = (39.0 * submoon_mass_kg * love_number * seconds * math.sqrt(constants.G)
               / (2.0 * (4.0 * math.pi * moon_density_kg_m3) ** (8.0 / 3.0) * tidal_q) * density_term ** (13.0 / 6.0))
    return bracket ** (1.0 / 3.0)


def moon_survives_tides(moon_radius_m, moon_density_kg_m3, planet_mass_kg, moon_orbit_radius_m, **kwargs):
    """Whether a moon of this radius passes the tidal-survival gate for a submoon (the gate row of Boss's card S6)."""
    return moon_radius_m >= tidal_survival_radius_m(moon_density_kg_m3, planet_mass_kg, moon_orbit_radius_m, **kwargs)


def _no_satellites_weight(subject):
    others = [row for row in SATELLITE_SHARE.passing(subject) if row.name != "no satellites"]
    return max(0.0, 1.0 - sum(row.evaluate(subject) for row in others))


_no_satellites_weight.__qualname__ = "no_satellites_weight[satellite-share]"


SATELLITE_SHARE = pt.REGISTRY.register(pt.Table(
    "satellite-share", 1, [
        pt.Row("inner rubble pile with a satellite", ("inner-small-body",), 0.15,
               source="Pravec et al. 2006 (Icarus) via Margot et al. 2015 (Asteroids IV, arXiv:1504.00034): 15 +- 4 percent of "
                      "near-Earth asteroids over 0.3 km are binary",
               state=pt.PARTLY_VERIFIED,
               limit="second-hand (the 2006 paper was not read); measured above 0.3 km and radar down to 0.12 km, the gate's "
                     "cut is 1 km; 0.11 to 0.19 within the quoted error",
               conditions=(pt.Condition("diameter_m", SIZE_CUT_M, PHOTOMETRIC_BINARY_MAX_DIAMETER_M),)),
        pt.Row("belt rubble pile with a satellite", ("belt-small-body",), 0.16,
               source="Jacobson et al. 2016 (arXiv:1605.06157): the modelled steady-state binary fraction 16 +- 6 percent; "
                      "Margot et al. 2015: small main-belt binaries are discovered about as often as near-Earth ones",
               state=pt.STAND_IN, replacement="a measured small-body binary fraction",
               conditions=(pt.Condition("diameter_m", SIZE_CUT_M, PHOTOMETRIC_BINARY_MAX_DIAMETER_M),)),
        pt.Row("large small body with a satellite", ("inner-small-body", "belt-small-body"), 0.05,
               source="Margot et al. 2015: 16 satellites among over 300 large main-belt asteroids surveyed (about 0.05), "
                      "detectability-limited, 'substantially smaller than 16 percent'",
               state=pt.STAND_IN, replacement="a debiased satellite fraction of large asteroids",
               conditions=(pt.Condition("diameter_m", PHOTOMETRIC_BINARY_MAX_DIAMETER_M, None),)),
        pt.Row("cold trans-neptunian binary", ("trans-neptunian-cold",), 0.22,
               source="Stephens and Noll 2006 via Noll et al. 2008 (astro-ph/0703134): 22 +10 -5 percent of classical "
                      "trans-neptunian objects under 5 degrees inclination have detectable binaries",
               state=pt.PARTLY_VERIFIED,
               limit="separations over 60 mas and primaries over about 100 km only; the share below that size is not measured",
               conditions=(pt.Condition("diameter_m", TRANS_NEPTUNIAN_MIN_DIAMETER_M, None),)),
        pt.Row("hot trans-neptunian binary", ("trans-neptunian-hot",), 0.055,
               source="Stephens and Noll 2006 via Noll et al. 2008: 5.5 +4 -2 percent of the other dynamical classes",
               state=pt.PARTLY_VERIFIED,
               limit="same detection limits as the cold row",
               conditions=(pt.Condition("diameter_m", TRANS_NEPTUNIAN_MIN_DIAMETER_M, None),)),
        pt.Row("no satellites", (pt.ANY_CLASS,), _no_satellites_weight,
               source="the rest of the rows above (1 minus the chance of a satellite)", state=pt.STAND_IN,
               replacement="measured satellite fractions by size")],
    description="Whether a small body, by population and size, has a satellite; the size cut is `SIZE_CUT_M`."))


def satellite_share(host):
    """
    The chance that `host` has satellites, from `satellite-share`, or `None` when the host's population has no rows
    (a planet or star: its host-class table decides). A body below `SIZE_CUT_M` has none. A moon has none: it needs a row
    that does not exist yet.
    """
    if host.population == "moon" or host.diameter_m < SIZE_CUT_M:
        return 0.0
    if host.population not in SATELLITE_SHARE_CLASSES:
        return None
    subject = pt.Subject(host.population, {"diameter_m": host.diameter_m})
    return math.fsum(row.evaluate(subject) for row in SATELLITE_SHARE.passing(subject) if row.name != "no satellites")


def can_host(host):
    """
    Whether `host` can have satellites, running the three parts of the gate in order (module docstring).

    Returns:
        Verdict: `allowed`, the first failing `part` (0 when it passes) and the reason.
    """
    if not satellite_fits(host):
        return Verdict(False, 1, "a floor-sized satellite does not fit: it would outweigh the host past the mass-ratio "
                                 "limit or the host is no bigger than the floor")
    if not stable_orbit_exists(host):
        return Verdict(False, 2, "no stable orbit: the stable zone between the host's radius (or Roche limit) and the "
                                 "Hill limit is empty, or the triple fails the stability screen")
    if host.diameter_m < SIZE_CUT_M:
        return Verdict(False, 3, f"below the size cut of {SIZE_CUT_M:.0f} m: bodies this small have no satellites")
    if host.population == "moon":
        return Verdict(False, 3, "a moon hosts only through a sourced satellite-share row, and none exists yet")
    share = satellite_share(host)
    if share is None:
        return Verdict(True, 0, f"this gate has no rows for population {host.population!r}; its host-class table decides")
    if share <= 0.0:
        return Verdict(False, 3, "the satellite share of this population and size is zero")
    return Verdict(True, 0, f"fits, has a stable zone, and {share:.3f} of such bodies have satellites")
