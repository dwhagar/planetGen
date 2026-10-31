# planetgen/physics/ocean_tides.py

"""
The open-ocean equilibrium tide at a world's surface (GEN.1034)
================================================================

The height a liquid surface would reach if it could follow the tidal potential of each body that raises one, summed over the moons and the star into the
spring range (every raiser in line) and the neap range (the largest against the rest). It is the *equilibrium* value only (Pugh and Woodworth 2014, *Sea-Level
Science*, chapter 3, *Tides: A Scientific History*, recalled): the range on a coast depends on the shape and depth of the basin and on its
resonance, which this does not model, so it is an upper bound for an ocean that does not resonate and can be exceeded many times over in one that does.

* **The raiser.** A body of mass `m` at distance `d` raises on a body of mass `M` and radius `R` a bulge of `(m / M) (R / d)^3 R` above the mean surface on the
  line to it and a third of that below on the line across, so a high-to-low range of `(3/2) (m / M) (R / d)^3 R`: 0.535 m for the Moon on the Earth and 0.246 m
  for the Sun, against the published 0.54 and 0.25 m.
* **Eccentric orbits.** The force goes as `1 / d^3`, so the range at the closest approach is `1 / (1 - e)^3` of the range at the semi-major axis and at the
  farthest `1 / (1 + e)^3`.
* **A raiser the body is locked to** (the body turns at the raiser's orbital rate, either way round) holds its bulge fixed in the sky of the body: it raises no
  oscillating tide, and the water stands higher on one side for good. Only raisers the body is not locked to enter the ranges (`oscillates`). A tide that an eccentric
  orbit still varies on a locked body (the radial tide) is not counted.
* **The floor.** The sea floor rises with the tide, so the sea measures `1 + k2 - h2` of it (Love 1909; 0.69 for the Earth with `k2` 0.30 and `h2` 0.61), with `h2` from
  `k2` of the interior model (`H2_OVER_K2`) unless it is given.
* **The liquid.** The equilibrium height does not depend on the liquid's density: the liquid and the potential both scale with it. What the density sets here is the
  loading of the floor by the sea's own weight and its self-attraction, which for a thin layer multiplies the tide by `1 / (1 - 3 rho_liquid / (5 rho_mean))` (recalled): 1.12 for water on the Earth. The factor goes singular as the liquid approaches 5/3 of the body's mean density, so a liquid heavier than
  `LOADING_REFUSE` of that is refused (`OceanTideError`): a global layer of a heavy liquid is not thin, and needs the fluid Love numbers, not this.
* **An ice-covered ocean** has no free surface: the lid flexes, and the free-surface tide is not defined (zero, with the reason).

The depth of the basin, the dissipation of the tide and non-water liquids (GEN.1035) are not here. Stand-ins are OT1, OT2 and OT6 to OT8 of
`docs/design/research/gap-filling-plan.md` section 6.68.

**The stored columns** (`FIELDS`, on `planets` and `moons`; `generate`, called where `daylight.generate` is, after the moons and the spin are final):

* `tide_range_spring_m` and `tide_range_neap_m`, metres: the open-ocean equilibrium ranges above, for a body whose hydrosphere is `surface ocean` (a free
  sea); `0` for an `ice-covered ocean` (no free surface, with the reason in `ocean_tide`); `None` for every other body (dry, vapour, ice, a brine seep, a
  hycean world, a gas giant, a body with no spin on record).
* `ocean_liquid`: the liquid of that sea, `water` until the surface liquid inventory (GEN.1035) names others; `None` where the tide columns are.

The raisers of a planet are its star and each of its moons; those of a moon are its planet and the star (at the planet's distance). The other planets'
tides are not counted (a planet at one astronomical unit raises on the Earth 1 percent of the Moon's). The spin sense is read against each raiser's orbit
normal: a raiser that goes round the other way is a retrograde one for the body and can never lock it. `k2` is the interior model's
(`_love_k2`): the homogeneous rigid-sphere figure for rock (`tidal_heating.rigid_k2`), the ice sphere's for an icy body (`icy_tide.love_k2`, GEN.960), a gas giant's `tidal_heating.GIANT_K2`.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, icy_tide, spin, tidal_heating

FIELDS = ("tide_range_spring_m", "tide_range_neap_m", "ocean_liquid")
"""tuple: The stored columns on `planets` and `moons` (see the module text)."""

H2_OVER_K2 = 5.0 / 3.0
"""float: `h2 / k2` of a homogeneous incompressible elastic sphere (Love 1909), used when only `k2` of the interior model is known. The Earth's measured values (0.61 and 0.30) give 2.0, so
the floor's share of the tide is 0.80 here against 0.69 measured. Source tier T5, tweak level L1."""
LOADING_REFUSE = 0.9
"""float: The most `3 rho_liquid / (5 rho_mean)` the thin-layer loading factor is used at (the factor is 10 there); a heavier liquid is refused. Source tier T2, tweak level L0."""
LOCK_TOLERANCE = 1.0e-3
"""float: A body whose spin is within this share of a raiser's orbital rate is locked to it (the same share as `moon_history.SYNCHRONOUS_TOLERANCE`). Source tier T6, tweak level L0."""


class OceanTideError(ValueError):
    """A liquid or a body the thin-layer equilibrium tide cannot describe."""


@dataclass(frozen=True)
class Liquid:
    """A liquid an ocean can hold: its name and density, kg/m3."""
    name: str
    density_kg_m3: float


WATER = Liquid("water", 1000.0)
"""Liquid: Pure water at 4 degrees (NIST); a sea of salt water is 2.5 percent denser, which moves the loading factor by 0.3 percent. Source tier T1, tweak level L0."""


@dataclass(frozen=True)
class Raiser:
    """A body that raises a tide: its mass, its distance (semi-major axis), its orbital eccentricity, and its orbital rate about the body (rad/s, negative for a retrograde orbit)."""
    name: str
    mass_kg: float
    distance_m: float
    eccentricity: float
    orbital_rate_rad_s: float

    @classmethod
    def of(cls, name, mass_kg, distance_m, body_mass_kg, eccentricity=0.0, retrograde=False):
        """The raiser of mass `mass_kg` on an orbit of semi-major axis `distance_m` about a body of `body_mass_kg` (the two-body rate), going round the other way when `retrograde`."""
        rate = math.sqrt(constants.G * (body_mass_kg + mass_kg) / distance_m ** 3)
        return cls(name, mass_kg, distance_m, eccentricity, -rate if retrograde else rate)


@dataclass(frozen=True)
class RaiserTide:
    """One raiser's tide: the equilibrium range (m) before the floor and the loading, the range after both at the semi-major axis and at the closest and farthest approach, and whether it
    oscillates (the body is not locked to it)."""
    name: str
    equilibrium_range_m: float
    range_m: float
    range_closest_m: float
    range_farthest_m: float
    oscillates: bool


@dataclass(frozen=True)
class OceanTide:
    """The tide of an ocean: each raiser's, the spring and neap ranges (m), the same at the extremes of the orbits (`spring_range_max_m`: every raiser closest; `neap_range_min_m`: the
    largest farthest and the rest closest), the factors used, the liquid, and the reason for a zero (None when the ocean has a free surface)."""
    liquid: Liquid
    love_factor: float
    loading_factor: float
    raisers: tuple
    spring_range_m: float
    neap_range_m: float
    spring_range_max_m: float
    neap_range_min_m: float
    reason: str = None


def equilibrium_range_m(raiser_mass_kg, distance_m, body_mass_kg, body_radius_m):
    """The high-to-low equilibrium tide, m, a raiser of `raiser_mass_kg` at `distance_m` makes on a body: `(3/2) (m / M) (R / d)^3 R`."""
    return 1.5 * raiser_mass_kg / body_mass_kg * (body_radius_m / distance_m) ** 3 * body_radius_m


def love_factor(k2, h2=None):
    """`1 + k2 - h2`: the share of the tide the sea measures against the floor that rises with it. `h2` is `H2_OVER_K2 * k2` when not given."""
    return 1.0 + k2 - (H2_OVER_K2 * k2 if h2 is None else h2)


def loading_factor(liquid_density_kg_m3, mean_density_kg_m3):
    """The thin-layer factor `1 / (1 - 3 rho_liquid / (5 rho_mean))` of the sea's own weight and attraction; `OceanTideError` for a liquid heavy enough (over `LOADING_REFUSE`) that a thin layer is no description."""
    ratio = 3.0 * liquid_density_kg_m3 / (5.0 * mean_density_kg_m3)
    if ratio > LOADING_REFUSE:
        raise OceanTideError(f"a liquid of {liquid_density_kg_m3:g} kg/m3 on a body of mean density {mean_density_kg_m3:g} kg/m3 is not a thin layer "
                             f"(3 rho_liquid / 5 rho_mean = {ratio:.2f} over {LOADING_REFUSE:g}); a global layer that heavy needs the fluid Love numbers")
    return 1.0 / (1.0 - ratio)


def oscillates(spin_rate_rad_s, raiser):
    """Whether a body spinning at `spin_rate_rad_s` (rad/s, negative for retrograde) sees `raiser` move: false when it turns at the raiser's orbital rate, within `LOCK_TOLERANCE`."""
    return abs(spin_rate_rad_s - raiser.orbital_rate_rad_s) > LOCK_TOLERANCE * abs(raiser.orbital_rate_rad_s)


def ranges(oscillating_ranges_m):
    """`(spring, neap)` for the ranges of the raisers that oscillate: the spring range is their sum (all in line), the neap range the largest less the rest (not below zero); one raiser alone has
    the same range at both, and no raiser has none."""
    if not oscillating_ranges_m:
        return 0.0, 0.0
    ordered = sorted(oscillating_ranges_m, reverse=True)
    spring = sum(ordered)
    neap = ordered[0] if len(ordered) == 1 else max(0.0, ordered[0] - sum(ordered[1:]))
    return spring, neap


def ocean_tide(body_mass_kg, body_radius_m, spin_rate_rad_s, raisers, k2, liquid=WATER, h2=None, ice_covered=False):
    """
    The open-ocean equilibrium tide of a world with an ocean of `liquid`, on a body of `body_mass_kg` and `body_radius_m` spinning at `spin_rate_rad_s` (rad/s, signed), raised by `raisers` (the
    star and the moons as `Raiser`s), with the `k2` of its interior (and `h2` if known). An ice-covered ocean has zero tide and says why; a liquid too heavy for the thin layer raises `OceanTideError`.
    """
    mean_density = body_mass_kg / (4.0 / 3.0 * math.pi * body_radius_m ** 3)
    love = love_factor(k2, h2)
    loading = loading_factor(liquid.density_kg_m3, mean_density)
    if ice_covered:
        return OceanTide(liquid, love, loading, (), 0.0, 0.0, 0.0, 0.0, "ice-covered ocean: the lid flexes, there is no free-surface tide")
    tides = []
    for raiser in raisers:
        base = equilibrium_range_m(raiser.mass_kg, raiser.distance_m, body_mass_kg, body_radius_m)
        mean = base * love * loading
        e = raiser.eccentricity or 0.0
        tides.append(RaiserTide(raiser.name, base, mean, mean / (1.0 - e) ** 3, mean / (1.0 + e) ** 3, oscillates(spin_rate_rad_s, raiser)))
    moving = [t for t in tides if t.oscillates]
    spring, neap = ranges([t.range_m for t in moving])
    spring_max = sum(t.range_closest_m for t in moving)
    if len(moving) > 1:
        largest = max(moving, key=lambda t: t.range_m)
        neap_min = max(0.0, largest.range_farthest_m - sum(t.range_closest_m for t in moving if t is not largest))
    else:
        neap_min = moving[0].range_farthest_m if moving else 0.0
    return OceanTide(liquid, love, loading, tuple(tides), spring, neap, spring_max, neap_min)


# --- The stored columns ------------------------------------------------------------------------------------------

def _liquid_of(body):
    """The liquid of the body's sea: water until the surface liquid inventory (GEN.1035) says otherwise."""
    return WATER


def _axis(body):
    axis = tuple(getattr(body, name, None) for name in ("spin_axis_x", "spin_axis_y", "spin_axis_z"))
    return None if any(value is None for value in axis) else axis


def _normal(body):
    return spin.orbit_normal(getattr(body, "orbital_inclination_deg", None) or 0.0, getattr(body, "orbital_ascending_node_deg", None) or 0.0)


def _retrograde(axis, normal):
    """Whether an orbit with this normal goes round the other way from the spin about `axis`."""
    return sum(a * n for a, n in zip(axis, normal)) < 0.0


def _love_k2(body, mass, radius_m, rate):
    """The Love number `k2` of the body whose sea is raised: a gas giant's fixed one, an icy body's from its ice sphere at the body's spin rate, else a rocky sphere's rigid one."""
    if getattr(body, "body_type", None) == "g":
        return tidal_heating.GIANT_K2
    if tidal_heating.is_icy(body):
        return icy_tide.love_k2(mass, radius_m, rate)
    return tidal_heating.rigid_k2(mass, radius_m)


def compute(body, parent=None):
    """
    The `FIELDS` of a planet or moon as a dict (`parent` the planet a moon orbits): a `surface ocean` gets its spring and neap ranges and its liquid, an
    `ice-covered ocean` zero ranges, any other body `None`. Pure: no draws.
    """
    out = dict.fromkeys(FIELDS)
    hydrosphere = getattr(body, "hydrosphere", None)
    star = getattr(body, "star", None)
    rotation = getattr(body, "rotation_period_hours", None)
    axis = _axis(body)
    if hydrosphere not in ("surface ocean", "ice-covered ocean") or star is None or not rotation or rotation <= 0.0 or axis is None:
        return out
    mass, radius_m = body.mass, body.radius * 1000.0
    host = parent if parent is not None else body
    rate = 2.0 * math.pi / (rotation * 3600.0)
    raisers = [Raiser.of(getattr(star, "name", None) or "star", star.mass, host.distance * constants.AU_TO_M, mass, getattr(host, "eccentricity", 0.0) or 0.0,
                         _retrograde(axis, _normal(host)))]
    if parent is not None:
        raisers.append(Raiser.of("planet", parent.mass, body.distance * constants.AU_TO_M, mass, getattr(body, "eccentricity", 0.0) or 0.0,
                                 _retrograde(axis, _normal(body))))
    for moon in getattr(body, "moons", None) or ():
        raisers.append(Raiser.of("moon", moon.mass, moon.distance * constants.AU_TO_M, mass, getattr(moon, "eccentricity", 0.0) or 0.0,
                                 _retrograde(axis, _normal(moon))))
    liquid = _liquid_of(body)
    tide = ocean_tide(mass, radius_m, rate, raisers, _love_k2(body, mass, radius_m, rate), liquid, ice_covered=hydrosphere == "ice-covered ocean")
    out["tide_range_spring_m"], out["tide_range_neap_m"], out["ocean_liquid"] = tide.spring_range_m, tide.neap_range_m, liquid.name
    return out


def generate(body, parent=None):
    """Sets `FIELDS` on a planet or moon (`parent` the planet a moon orbits)."""
    for name, value in compute(body, parent).items():
        setattr(body, name, value)
