# planetgen/physics/moon_composition.py

"""
Where a moon came from and what it is made of (GEN.466, docs/design/small-bodies-composition.md section 5.4 and 5.5)
====================================================================================================================

A moon takes an origin from a gated table (`ORIGIN_TABLE`), and the origin leaves a composition:

* **disk**: a regular moon of a giant planet, formed in its circumplanetary disk. The ice fraction rises with the formation
  distance past the disk's ice line (`ice_line_km`, Sasaki et al. 2010; the ramp of decision S8, `disk_ice_fraction`), and the rock
  turns from anhydrous through partly hydrated to carbonaceous with the same distance. A disk system (`draw_disk_system`) holds the
  disk lifetime, the total moon mass of Canup and Ward 2009 and its partition among the large moons.
* **impact**: a terrestrial planet's moon from a giant impact (Lock et al. 2018): the condensate has little metal, no ice and the
  planet's rock.
* **captured**: a body caught on an irregular orbit, prograde or retrograde (`draw_inclination_deg`): the icy-body rule (`draw_icy_body`) with ice lost by U(0, 0.5).
* **ring**: a moon that formed from a massive ring near the Roche radius: ice-rich rock and ice.
* **binary**: the partner of a dwarf planet: the icy-body rule, with a mass ratio of 5e-4 to 0.11.

Every weight of the origin table is open (no source): the rows a moon supports share the probability equally. The gate is the
parent's kind, whether the parent has a massive ring and the moon's distance in Roche radii. The orbit follows the origin
(GEN.714): a regular moon lies near the parent's equator, a captured one anywhere outside the Kozai-unstable band.

What each draw returns is a `SmallBodyMakeup` (see `small_body_composition`) with `detail["origin"]`. The density of a moon below
`compaction.TAPER_END_KG` is `compacted_density_kg_m3(radius, rock_fraction)`; above it the structure solver of lane I decides.

Use:

    origin = draw_origin("giant")                                    # a Possibility; origin.row.name is "disk" or "captured"
    inclination = draw_inclination_deg(origin.row.name)              # degrees to the parent's equator
    system = draw_disk_system(parent_mass_kg)                        # tau, alpha/f, masses, ice line, f_inf
    moon = draw_moon(origin.row.name, formation_km, system)          # SmallBodyMakeup
"""

import math

from planetgen.generation import possibility_table as pt
from planetgen.physics import compaction, constants, small_body_composition as sb
from planetgen.util import draw

GIANT = "giant"
TERRESTRIAL = "terrestrial"
DWARF_PLANET = "dwarf-planet"
PARENT_KINDS = (GIANT, TERRESTRIAL, DWARF_PLANET)
ORIGINS = ("disk", "impact", "captured", "ring", "binary")

GIANT_ENVELOPE_FRACTION = 0.01
"""float: A parent with at least this envelope mass fraction is a giant (a parametrisation; an envelope of 1 percent is where a planet
has a disk of its own)."""
DWARF_PLANET_MAX_MASS_KG = 5e22
"""float: A parent without an envelope below this mass is a dwarf planet (a stand-in: Pluto and Eris are 1.3e22 and 1.7e22)."""
KOZAI_FORBIDDEN_DEG = (55.0, 130.0)
"""tuple: The inclinations (to the parent's equator, degrees) at which a distant orbit is unstable and no irregular moon is found:
between them the Kozai resonance stretches the eccentricity until the moon escapes or hits something (Nesvorny et al. 2003, Icarus
126, 398; the 55 to 130 degree bounds as used by Jewitt's group for Jupiter, arXiv:0805.3733; Shen and Tremaine 2008 give 60 to
130). Partly verified: the bounds are a recalled range, not a measured edge."""
RING_MAX_ROCHE_RATIO = 3.0
"""float: A moon forms from a ring only within this many Roche radii (a parametrisation)."""
CAPTURED_FULL_WEIGHT_MASS_KG = 1e16
"""float: STAND-IN CM17 (tier T9, L1): a moon of this mass or less keeps the full weight of the `captured` row. The storage floor of the irregular list (CM15), the
size from which the irregular moons of the four giants are all known."""
CAPTURED_SIZE_INDEX = 1.0
"""float: STAND-IN CM17 (tier T1, L0): the cumulative size index of the irregular moons that survive, N(>R) proportional to R^-1 (differential index 1.9 +/- 0.3 at Jupiter and
Saturn above 5 km, CM5; JH07). A mass M goes as R^3, so the weight of a captured moon of mass M falls as (M / floor)^(-index / 3)."""


def parent_kind(envelope_fraction, mass_kg):
    """`giant`, `terrestrial` or `dwarf-planet` for a parent of this envelope mass fraction and mass."""
    if envelope_fraction >= GIANT_ENVELOPE_FRACTION:
        return GIANT
    return DWARF_PLANET if mass_kg < DWARF_PLANET_MAX_MASS_KG else TERRESTRIAL


def _row(name, classes, conditions, rule, source):
    return pt.Row(name=name, classes=classes, weight=pt.OPEN, source=pt.OPEN, state=pt.OPEN_STATE, conditions=conditions,
                  notes=f"{rule}. {source}. The weight is open: the rows a moon supports share the probability equally.")


ORIGIN_TABLE = pt.REGISTRY.register(pt.Table(
    "moon-origin", 2, [
        _row("disk", (GIANT,), (), "regular moon of a circumplanetary disk, composition by the ice-line ramp",
             "Canup and Ward 2009, Sasaki et al. 2010, Cilibrasi et al. 2018"),
        _row("impact", (TERRESTRIAL,), (), "giant-impact moon: little metal, no ice, the parent's rock",
             "Lock et al. 2018 (verified for the rule; the metal range is a parametrisation)"),
        _row("captured", (pt.ANY_CLASS,), (),
             "captured body: the icy-body rule with ice lost by U(0, 0.5)", "Barr and Schwamb 2016 (qualitative); the factor is a parametrisation"),
        _row("ring", (GIANT,), (pt.Condition("has_ring", equals=1), pt.Condition("roche_ratio", maximum=RING_MAX_ROCHE_RATIO)),
             "ring-born moon: ice-rich near the Roche radius", "Crida and Charnoz 2012 (abstract only, recalled for the composition)"),
        _row("binary", (DWARF_PLANET,), (), "partner of a dwarf planet: the icy-body rule, mass ratio 5e-4 to 0.11",
             "Barr and Schwamb 2016 Table 2"),
    ], description="The origin of a moon, gated by the parent's kind, the parent's ring and the distance in Roche radii."))
"""Table: The moon origin rows (every weight open)."""


def draw_origin(kind, has_ring=False, roche_ratio=math.inf):
    """
    The origin of a moon of a parent of `kind`. A `Possibility`: `.row.name` is `disk`, `impact`, `captured`, `ring` or `binary`.
    The origin comes first and the orbit follows from it (`draw_inclination_deg`).

    Raises:
        PossibilityError: No row supports the moon (cannot happen for the three kinds).
    """
    if kind not in PARENT_KINDS:
        raise ValueError(f"the parent kind {kind!r} is none of {PARENT_KINDS}")
    return ORIGIN_TABLE.draw(kind, {"has_ring": 1 if has_ring else 0, "roche_ratio": roche_ratio})


def captured_weight(mass_kg):
    """
    How much the `captured` row weighs against a row of weight 1 for a moon of `mass_kg` (GEN.1081, a pure function; the moon generator weights the row with it).

    Captured bodies are small: the survivors follow a size law N(>R) ~ R^-1 (`CAPTURED_SIZE_INDEX`), so the share of captured bodies above a mass falls as that mass
    to the power of -1/3. A moon of `CAPTURED_FULL_WEIGHT_MASS_KG` or less weighs 1 (no change from the open weights); one of 1e20 kg (R about 200 km) weighs 0.046, one
    of 1e22 kg 0.010 and Triton's 2.1e22 kg 0.0078: rarely captured, with the exchange route of a Triton a separate window (design section 2.2).

    Raises:
        ValueError: The mass is not positive.
    """
    if not mass_kg > 0.0:
        raise ValueError(f"the moon mass {mass_kg!r} is not positive")
    if mass_kg <= CAPTURED_FULL_WEIGHT_MASS_KG:
        return 1.0
    return (mass_kg / CAPTURED_FULL_WEIGHT_MASS_KG) ** (-CAPTURED_SIZE_INDEX / 3.0)


def origin_probabilities(kind, mass_kg, has_ring=False, roche_ratio=math.inf):
    """
    `[(origin name, probability)]` of a moon of `mass_kg` of a parent of `kind`: the origin table's probabilities (`draw_origin`) with the `captured` row multiplied by
    `captured_weight(mass_kg)` and the rest renormalised, in the table's row order. The probabilities add up to 1. Pure: no draw.

    Raises:
        ValueError: The parent kind is not one of `PARENT_KINDS` or the mass is not positive.
    """
    if kind not in PARENT_KINDS:
        raise ValueError(f"the parent kind {kind!r} is none of {PARENT_KINDS}")
    weight = captured_weight(mass_kg)
    table = ORIGIN_TABLE.probabilities(pt.Subject(kind, {"has_ring": 1 if has_ring else 0, "roche_ratio": roche_ratio}))
    scaled = [(row.name, p * weight if row.name == "captured" else p) for row, p in table]
    total = math.fsum(p for _name, p in scaled)
    return [(name, p / total) for name, p in scaled]


def draw_inclination_deg(origin):
    """
    The inclination of a moon's orbit to its parent's equator, degrees, from its origin (GEN.714).

    A regular moon (disk, impact, ring, binary) lies near the equator: uniform on 0 to `constants.MOON_ORBITAL_INCLINATION_MAX_DEG`,
    the generator's long-standing bound (the regular satellites of the giants lie within a degree or so of the equator, the
    Galilean moons within 0.5 degrees, so the bound is wider than what is observed: a stand-in). A captured moon takes an
    isotropic orbit, the cosine of its inclination uniform from -1 to 1, outside the Kozai-unstable band `KOZAI_FORBIDDEN_DEG`,
    so it is prograde below the band or retrograde above it (Nesvorny et al. 2003). The draw is one number either way: the
    uniform number is mapped onto the allowed cosines, so no value is rejected.

    Args:
        origin (str): `disk`, `impact`, `captured`, `ring` or `binary`.
    """
    if origin not in ORIGINS:
        raise ValueError(f"the origin {origin!r} is none of {ORIGINS}")
    if origin != "captured":
        return draw.uniform(0.0, constants.MOON_ORBITAL_INCLINATION_MAX_DEG)
    low, high = (math.radians(angle) for angle in KOZAI_FORBIDDEN_DEG)
    prograde = 1.0 - math.cos(low)          # the cosines from 1 down to cos(low)
    retrograde = 1.0 + math.cos(high)       # the cosines from cos(high) down to -1
    point = draw.random() * (prograde + retrograde)
    cosine = 1.0 - point if point < prograde else math.cos(high) - (point - prograde)
    return math.degrees(math.acos(min(max(cosine, -1.0), 1.0)))


# ---------------------------------------------------------------------------------------------------------------------
# The disk system

JUPITER_MASS_KG = 1.89813e27
JUPITER_RADIUS_KM = 71492.0
DISK_LIFETIME_MYR = (1.0, 10.0)
"""tuple: The disk lifetime, log-uniform (Sasaki used 5 Myr as the reference)."""
ALPHA_OVER_F = (1e-6, 5e-4)
"""tuple: The ratio of the viscosity parameter to the disk's gas-to-solid factor, log-uniform (Canup and Ward 2009)."""
ICE_LINE_VALID_MASS_JUPITER = (0.1, 3.0)
"""tuple: The parent masses, in Jupiter masses, inside which the ice-line formula is not outside what it was fitted for."""
LARGE_MOONS = (1, 5)
"""tuple: The number of large moons of a disk system (a parametrisation)."""


def ice_line_km(parent_mass_kg, lifetime_myr):
    """
    The ice line of a circumplanetary disk, km: `20 R_J (M_p / M_J)^(2/3) (tau / 5 Myr)^(-1/3)` (Sasaki et al. 2010). The formula was
    fitted for a Jupiter-like disk; outside `ICE_LINE_VALID_MASS_JUPITER` it is an extrapolation (partly verified, limit named).
    """
    if not parent_mass_kg > 0.0 or not lifetime_myr > 0.0:
        raise ValueError("the mass and the lifetime must be positive")
    return 20.0 * JUPITER_RADIUS_KM * (parent_mass_kg / JUPITER_MASS_KG) ** (2.0 / 3.0) * (lifetime_myr / 5.0) ** (-1.0 / 3.0)


def total_moon_mass_ratio(alpha_over_f):
    """The mass of all the large moons over the parent's: `2.5e-4 (alpha/f / 3e-5)^(1/3)` (Canup and Ward 2009, Eq. 4.15)."""
    return 2.5e-4 * (alpha_over_f / 3e-5) ** (1.0 / 3.0)


class DiskSystem:
    """
    The circumplanetary disk of a giant planet and the moons it makes.

    Attributes:
        parent_mass_kg, lifetime_myr, alpha_over_f (float): The inputs of the draw.
        ice_line_km (float): The distance of the ice line from the parent.
        extrapolated (bool): True where the parent's mass lies outside `ICE_LINE_VALID_MASS_JUPITER`.
        f_inf (float): The asymptotic ice fraction of the ramp (`small_body_composition.F_INF`).
        moon_masses_kg (list): The masses of the large moons, summing to the total of Canup and Ward.
    """

    def __init__(self, parent_mass_kg, lifetime_myr, alpha_over_f, f_inf, moon_masses_kg):
        self.parent_mass_kg = parent_mass_kg
        self.lifetime_myr = lifetime_myr
        self.alpha_over_f = alpha_over_f
        self.ice_line_km = ice_line_km(parent_mass_kg, lifetime_myr)
        low, high = ICE_LINE_VALID_MASS_JUPITER
        self.extrapolated = not low <= parent_mass_kg / JUPITER_MASS_KG <= high
        self.f_inf = f_inf
        self.moon_masses_kg = list(moon_masses_kg)

    def to_dict(self):
        return {"parent_mass_kg": self.parent_mass_kg, "lifetime_myr": self.lifetime_myr, "alpha_over_f": self.alpha_over_f,
                "f_inf": self.f_inf, "moon_masses_kg": list(self.moon_masses_kg)}

    @classmethod
    def from_dict(cls, data):
        return cls(data["parent_mass_kg"], data["lifetime_myr"], data["alpha_over_f"], data["f_inf"], data["moon_masses_kg"])


DISK_FIELDS = ("moon_disk_lifetime_myr", "moon_disk_alpha_over_f", "moon_disk_f_inf", "moon_disk_ice_line_km")
"""tuple: The columns of a giant planet that keep its circumplanetary disk (DB.70), as stored in `planets` (NULL for a body with no disk)."""
BODY_FIELDS = ("moon_origin", *DISK_FIELDS)
"""tuple: What the moon composition stores in `planets` and `moons` (DB.70): a moon's origin and its giant parent's disk."""


def store_disk(planet, system):
    """Sets `DISK_FIELDS` on the giant `planet` from its `DiskSystem` (the moons' share of the disk's ice line, lifetime and ramp)."""
    planet.moon_disk_lifetime_myr = system.lifetime_myr
    planet.moon_disk_alpha_over_f = system.alpha_over_f
    planet.moon_disk_f_inf = system.f_inf
    planet.moon_disk_ice_line_km = system.ice_line_km


def draw_disk_system(parent_mass_kg):
    """
    The disk system of a giant of this mass. The draws, in order: the lifetime, alpha/f (both log-uniform), the number of large moons,
    the Dirichlet partition of the total moon mass, then the asymptotic ice fraction.
    """
    if not parent_mass_kg > 0.0:
        raise ValueError("the parent mass must be positive")
    lifetime = math.exp(draw.uniform(math.log(DISK_LIFETIME_MYR[0]), math.log(DISK_LIFETIME_MYR[1])))
    alpha = math.exp(draw.uniform(math.log(ALPHA_OVER_F[0]), math.log(ALPHA_OVER_F[1])))
    count = draw.randint(*LARGE_MOONS)
    shares = draw.dirichlet([1.0] * count)
    total = total_moon_mass_ratio(alpha) * parent_mass_kg
    f_inf = sb._truncated_normal(*sb.F_INF)
    return DiskSystem(parent_mass_kg, lifetime, alpha, f_inf, [total * share for share in shares])


# ---------------------------------------------------------------------------------------------------------------------
# The compositions

HYDRATION_START = 0.4
"""float: Rock formed inside this fraction of the ice line is anhydrous."""
HYDRATION_END = 1.0
"""float: Rock formed beyond this fraction of the ice line is carbonaceous (fully hydrated)."""
IMPACT_METAL = (0.01, 0.05)
"""tuple: The metal mass fraction of an impact moon, uniform (a parametrisation: the condensate has no metal phase, iron enters as FeO)."""
IMPACT_WATER_PPM = (1.0, 1000.0)
"""tuple: The water of an impact moon, ppm by mass, log-uniform (Lock et al. 2018)."""
CAPTURE_ICE_LOSS = (0.0, 0.5)
"""tuple: The share of its ice a captured body loses, uniform (a parametrisation)."""
RING_ICE = (0.5, 0.95)
"""tuple: The ice mass fraction of a ring-born moon, uniform (a parametrisation, recalled)."""
BINARY_MASS_RATIO = (5e-4, 0.11)
"""tuple: The mass ratio of a dwarf planet's partner, log-uniform (Barr and Schwamb 2016 Table 2)."""

_ANHYDROUS = sb.TAXON_DATA["S"].fractions
_HYDRATED = {name: value for name, value in sb.TAXON_DATA["C"].fractions.items()}


def hydration(r_over_r_ice):
    """The share of the rock that is carbonaceous (0 anhydrous to 1), rising linearly in the ratio from `HYDRATION_START` to `HYDRATION_END`."""
    return min(max((r_over_r_ice - HYDRATION_START) / (HYDRATION_END - HYDRATION_START), 0.0), 1.0)


def _mix(weight_hydrated):
    mixed = {}
    for name, value in _ANHYDROUS.items():
        mixed[name] = mixed.get(name, 0.0) + value * (1.0 - weight_hydrated)
    for name, value in _HYDRATED.items():
        mixed[name] = mixed.get(name, 0.0) + value * weight_hydrated
    return mixed


def draw_disk_moon(formation_km, system):
    """
    A regular moon that formed `formation_km` from its parent in `system`: ice by the ramp of decision S8 at `formation_km / ice line`,
    the rest rock, anhydrous inside 0.4 of the ice line and carbonaceous beyond it. The formation distance is not the present one
    (moons migrate).

    Returns:
        SmallBodyMakeup: With `detail` holding the origin, the ratio, the hydration and the system's `f_inf`; the density depends
        on the mass (`density_kg_m3`).
    """
    if not formation_km > 0.0:
        raise ValueError("the formation distance must be positive")
    ratio = formation_km / system.ice_line_km
    ice = sb.disk_ice_fraction(ratio, system.f_inf)
    wet = hydration(ratio)
    rock = _mix(wet)
    total = math.fsum(rock.values())
    components = {name: (1.0 - ice) * value / total for name, value in rock.items()}
    if ice > 0.0:
        components["water_ice"] = components.get("water_ice", 0.0) + ice
    low, high = sb.POROSITY_RANGES[sb.ROCKY]
    return sb._makeup("moon", components, draw.uniform(low, high), detail={
        "origin": "disk", "r_over_r_ice": ratio, "hydration": wet, "f_inf": system.f_inf}, size_dependent=True)


def draw_impact_moon():
    """A giant-impact moon: the anhydrous rock of the parent with a small metal share, no ice and a trace of water."""
    metal = draw.uniform(*IMPACT_METAL)
    water = math.exp(draw.uniform(math.log(IMPACT_WATER_PPM[0]), math.log(IMPACT_WATER_PPM[1]))) * 1e-6
    silicate = {name: value for name, value in _ANHYDROUS.items() if sb.group_of(name) == "rock"}
    total = math.fsum(silicate.values())
    components = {name: (1.0 - metal - water) * value / total for name, value in silicate.items()}
    components["kamacite"] = metal
    components["water_ice"] = water
    low, high = sb.POROSITY_RANGES[sb.ROCKY]
    return sb._makeup("moon", components, draw.uniform(low, high), detail={"origin": "impact", "metal": metal, "water": water},
                      size_dependent=True)


def draw_captured_moon():
    """A captured body: the rock and ice of a Kuiper belt object (`draw_icy_body`), with its ice reduced by U(0, 0.5)."""
    body = sb.draw_icy_body(sb.KUIPER_BELT_OBJECT)
    loss = draw.uniform(*CAPTURE_ICE_LOSS)
    components = dict(body.components)
    components["water_ice"] = components.get("water_ice", 0.0) * (1.0 - loss)
    return sb._makeup("moon", components, body.porosity, detail={"origin": "captured", "ice_loss": loss,
                                                                  "rock_fraction": body.detail["rock_fraction"]}, size_dependent=True)


def draw_ring_moon():
    """A moon formed from a ring: ice-rich, U(0.5, 0.95) ice and anhydrous rock."""
    ice = draw.uniform(*RING_ICE)
    rock = _mix(0.0)
    total = math.fsum(rock.values())
    components = {name: (1.0 - ice) * value / total for name, value in rock.items()}
    components["water_ice"] = components.get("water_ice", 0.0) + ice
    return sb._makeup("moon", components, 0.0, detail={"origin": "ring", "ice": ice}, size_dependent=True)


def draw_binary_partner():
    """The partner of a dwarf planet: the icy-body rule, and a mass ratio log-uniform in `BINARY_MASS_RATIO` in `detail["mass_ratio"]`."""
    body = sb.draw_icy_body(sb.DWARF_PLANET)
    ratio = math.exp(draw.uniform(math.log(BINARY_MASS_RATIO[0]), math.log(BINARY_MASS_RATIO[1])))
    return sb._makeup("moon", body.components, body.porosity, detail={"origin": "binary", "mass_ratio": ratio,
                                                                       "rock_fraction": body.detail["rock_fraction"]}, size_dependent=True)


def draw_moon(origin, formation_km=None, system=None):
    """The composition of a moon of this origin (`disk` needs `formation_km` and `system`)."""
    if origin == "disk":
        if formation_km is None or system is None:
            raise ValueError("a disk moon needs its formation distance and its disk system")
        return draw_disk_moon(formation_km, system)
    drawers = {"impact": draw_impact_moon, "captured": draw_captured_moon, "ring": draw_ring_moon, "binary": draw_binary_partner}
    try:
        return drawers[origin]()
    except KeyError:
        raise ValueError(f"the origin {origin!r} is none of disk, {', '.join(drawers)}") from None


def density_kg_m3(makeup, mass_kg):
    """
    The bulk density of a moon of this makeup and mass below `compaction.TAPER_END_KG`: the tapered brittle law on its rock fraction
    (metal counted with the rock). Above that mass the structure solver of GEN.429 decides: `compaction.OutOfRange`.
    """
    low, high = compaction.ROCK_FRACTION_RANGE
    rock = min(max(makeup.rock + makeup.metal, low), high)
    radius_km = compaction.compacted_radius_km(mass_kg, rock)
    return mass_kg / (4.0 / 3.0 * math.pi * (radius_km * 1e3) ** 3)



def mass_from_radius_kg(makeup, radius_km):
    """
    The mass of a moon of this makeup and radius below `compaction.TAPER_END_KG`: the tapered brittle law on its rock fraction (metal
    counted with the rock) gives the density at that radius. Heavier: `compaction.OutOfRange`, and the structure solver of GEN.429
    decides.
    """
    low, high = compaction.ROCK_FRACTION_RANGE
    rock = min(max(makeup.rock + makeup.metal, low), high)
    return compaction.compacted_density_kg_m3(radius_km, rock) * 4.0 / 3.0 * math.pi * (radius_km * 1e3) ** 3
