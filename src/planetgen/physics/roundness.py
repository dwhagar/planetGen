# planetgen/physics/roundness.py

"""
Roundness: does a solid body's own weight beat its strength? (GEN.1158; /mnt/project-files/reports/round-shape-minimum-mass-2026-10-31.md)
==============================================================================================================================================

A solid body is round when the stress its own gravity puts on its inside is bigger than the strength of the material, so the lumps flow until the body settles into a ball.
Lineweaver and Norman 2010 (arXiv:1004.1091, READ) work this out: the pressure half way down a body of radius `R` and density `rho` is `(2 pi / 3) G rho^2 (R^2 - r^2)`, and
a body with a bump on it (the longest axis `f` times the mean radius) is made level when the difference of that pressure under the bump and beside it reaches the yield strength
`sigma`: `(2 pi / 3) G rho^2 R^2 (f^2 - 1) >= sigma` (their Eqs. 7 to 9). Solved for the radius it gives the **round radius**

    R_round = sqrt(sigma / (K G rho^2)),    K = (2 pi / 3) (f^2 - 1)

(with `f = 1.1`, a longest axis at most 10 percent over the mean radius, `K` is 0.44 and a 10 MPa body of 2.5 g/cm3 rounds at 233 km, their Eq. 10's 240 km).

1. **The roundness ratio** (`roundness_ratio`) is the body's mean radius over its round radius. At 1 or more the body is round (`is_round`, `ROUND_RATIO`); below 1 it is not, and the
   ratio says how far from round (0.5 is a body half the round size). The ratio is a smooth number, not a yes or no, and it rises with size and with density.
2. **The yield strength** (`yield_strength_pa`) of a body is not measured for any body but a few: Lineweaver and Norman find that icy moons turn round at a mean radius of about 200 km and
   rocky asteroids at about 300 km (their Fig. 2: Mimas 198 km round, Hyperion 140 km not; Ceres 475, Pallas 275 and Vesta 265 round, Juno 115 not), at bulk densities near 1 and 2.5 g/cm3. The
   strength of ice and of rock is fixed so that the round radius is exactly that at those densities (1.2 MPa and 16.5 MPa). That is a strength during the warm period when the shape was set,
   not today's (they find 14 times more for rock than ice then, against 2 times now).

   **Where real data drops in.** The strengths are not in this module but in `data/roundness_strengths.json` (`TABLE`): one entry per material (ice, rock, metal), either a measured
   `strength_pa` or the round radius and density it was seen at, plus the name of a mixing rule (`MIXING_RULES`: how the ice, rock and metal shares of a body combine) and of a porosity rule
   (`POROSITY_RULES`: what voids do to the strength). The engine's interface is the body's ice fraction, metal fraction and porosity in, and the strength, the round radius and the roundness
   ratio out; the strength research gives new numbers in the file, a new rule in one of the three dicts, and nothing else changes. The numbers in the file now are the strength research's
   (/mnt/project-files/reports/planetary-material-strength-2026-10-31.md, section 7): rock 16 MPa (confirmed three ways, range 4 to 30), warm ice 1.2 MPa and colder ice stronger by a power law
   of its temperature when the shape was set (`temperature_k`, left out for ice that was warm), metal 100 MPa (weak source, range 30 to 300), and a mix that follows the *volume* shares: ice
   strength while rock is under 40 percent of the volume, rock strength over 70 percent, geometric between. Porosity is `none` by default: the research proposes an exponential weakening, but
   with it Proteus (a porous ice-rich moon, 210 km, not round) would come out round, because the weaker strength outweighs the lower density; the rule is in `POROSITY_RULES` as `exponential`
   and is off until the rubble-pile branch (friction plus a little cohesion, no size threshold) is built.
3. **Binding energy** is in `binding_energy.py`.

The ice fraction is the share of the body's mass that is water or ice (`makeup.water + makeup.ice` of a small body). Only solid bodies have a roundness or a strength: a gas giant or a star
is a fluid and always settles into a ball, so nothing here applies to it.

Test (not an input): the ten bodies of the report, Mimas, Miranda, Enceladus, Proteus, Hyperion, Ceres, Vesta, Pallas and Juno fall on the right side; Hygiea (round at about 217 km with
a density of about 1.9 g/cm3, from Vernazza et al. 2019, not opened) falls on the wrong side, because the rock strength is fixed by the 300 km of the larger asteroids.

Use:

    ratio = roundness_ratio(radius_m, density_kg_m3, ice_fraction)           # 1.0 or more is round; metal_fraction, porosity and temperature_k are optional
    sigma = yield_strength_pa(ice_fraction, metal_fraction, porosity, temperature_k)
    ratio = roundness_ratio_from_mass(mass_kg, density_kg_m3, ice_fraction)
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen.physics import constants
from planetgen.util.checks import finite_domain

DATA_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "roundness_strengths.json")
"""str: The strengths of the materials and the rules that mix them."""

ROUND_RATIO = 1.0
"""float: A body is round at a roundness ratio of at least this (the round radius is already the 10 percent sphericity of the table, so one is attainable). A definition, not a stand-in."""

MATERIALS = ("ice", "rock", "metal")

GRAIN_DENSITIES_KG_M3 = {"ice": 920.0, "rock": 3000.0, "metal": 7800.0}
"""dict: The density of each material without voids, kg/m3, for the volume shares of a mix (water ice 917, a typical silicate grain, iron). [Stand-in RN6, tier T2, tweak L0.]"""


class RoundnessError(ValueError):
    """A strength table that is refused."""


@dataclass(frozen=True)
class StrengthTable:
    """The strength of each material in Pa, the sphericity of a round body, and the rules (by name and numbers) that mix materials, warm or cool ice and account for porosity."""

    strengths_pa: dict
    sphericity_factor: float
    mixing: dict
    porosity: dict
    temperature: dict


def stress_coefficient(sphericity_factor):
    """`K = (2 pi / 3)(f^2 - 1)`: the differential stress is `K G rho^2 R^2` for a longest axis `f` times the mean radius (Lineweaver and Norman 2010, Eqs. 7 and 8)."""
    return (2.0 * math.pi / 3.0) * (sphericity_factor ** 2 - 1.0)


def strength_from_round_radius_pa(round_radius_km, density_kg_m3, sphericity_factor):
    """The yield strength at which a body of this density turns round at this mean radius: `K G rho^2 R^2`."""
    return stress_coefficient(sphericity_factor) * constants.G * density_kg_m3 ** 2 * (round_radius_km * 1.0e3) ** 2


def _positive(value, name):
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not value > 0.0 or math.isinf(value):
        raise RoundnessError(f"{name} must be a positive number")
    return float(value)


def load_table(document):
    """The `StrengthTable` of a document in the shape of `data/roundness_strengths.json`, refused with `RoundnessError` where it is wrong."""
    factor = _positive(document.get("sphericity_factor"), "sphericity_factor")
    if factor <= 1.0:
        raise RoundnessError("the sphericity factor is over 1 (the longest axis is longer than the mean radius)")
    materials = document.get("materials")
    if not isinstance(materials, dict) or set(materials) != set(MATERIALS):
        raise RoundnessError(f"'materials' needs exactly {MATERIALS}")
    strengths = {}
    for name in MATERIALS:
        entry = materials[name]
        if "strength_pa" in entry:
            strengths[name] = _positive(entry["strength_pa"], f"{name} strength_pa")
        elif "round_radius_km" in entry:
            strengths[name] = strength_from_round_radius_pa(
                _positive(entry["round_radius_km"], f"{name} round_radius_km"), _positive(entry.get("reference_density_kg_m3"), f"{name} reference_density_kg_m3"), factor)
    for name in MATERIALS:
        if name not in strengths:
            target = materials[name].get("same_as")
            if target not in strengths:
                raise RoundnessError(f"{name} needs a strength or a 'same_as' material that has one")
            strengths[name] = strengths[target]
    mixing, porosity, temperature = document.get("mixing") or {}, document.get("porosity") or {}, document.get("temperature") or {}
    if mixing.get("rule") not in MIXING_RULES:
        raise RoundnessError(f"the mixing rule must be one of {sorted(MIXING_RULES)}")
    if porosity.get("rule") not in POROSITY_RULES:
        raise RoundnessError(f"the porosity rule must be one of {sorted(POROSITY_RULES)}")
    if temperature.get("rule") not in TEMPERATURE_RULES:
        raise RoundnessError(f"the temperature rule must be one of {sorted(TEMPERATURE_RULES)}")
    if mixing["rule"] == "volume-threshold":
        low, high = _positive(mixing.get("skeleton_below"), "skeleton_below"), _positive(mixing.get("skeleton_above"), "skeleton_above")
        if not low < high <= 1.0:
            raise RoundnessError("the volume thresholds need 0 < skeleton_below < skeleton_above <= 1")
    return StrengthTable(strengths, factor, mixing, porosity, temperature)


def _ice_power_law(table, ice_strength_pa, temperature_k):
    """Ice that took its shape cold is stronger: `coefficient (|T in C| / 10)^exponent`, capped at `cap_pa`, and never under the warm value (Schulson and Renshaw 2022, two points of the
    lab strength scaled to a body by 0.1; frozen soil follows the same exponent, Wan et al. 2024)."""
    if temperature_k is None or temperature_k >= constants.CELSIUS_ZERO_K:
        return ice_strength_pa
    rule = table.temperature
    cold = _positive(rule["coefficient_pa"], "coefficient_pa") * ((constants.CELSIUS_ZERO_K - temperature_k) / 10.0) ** _positive(rule["exponent"], "exponent")
    return min(max(ice_strength_pa, cold), _positive(rule["cap_pa"], "cap_pa"))


TEMPERATURE_RULES = {
    "none": lambda table, ice_strength_pa, temperature_k: ice_strength_pa,
    "ice-power-law": _ice_power_law,
}
"""dict: Rules for the strength of ice that took its shape at a temperature, each `(table, warm_ice_strength_pa, temperature_k)` -> Pa. `temperature_k` of None is ice that was warm. [Stand-in RN7, tier T5, tweak L2.]"""


def _volume_threshold(table, ice_fraction, metal_fraction, temperature_k):
    """The strength of a mix by the volume shares of its materials (the strength research, section 7): the solids' own strength (rock and metal, geometric mean weighted by the metal share of
    their volume) once the solids are `skeleton_above` of the volume or more, ice strength while they are `skeleton_below` or less, geometric between. Ice is the load-bearing matrix until the
    solid grains touch (hard grains only strengthen ice at high volume fractions: Durham 1997, abstract)."""
    ice_strength = TEMPERATURE_RULES[table.temperature["rule"]](table, table.strengths_pa["ice"], temperature_k)
    rock_fraction = 1.0 - ice_fraction - metal_fraction
    volumes = {name: share / GRAIN_DENSITIES_KG_M3[name] for name, share in (("ice", ice_fraction), ("rock", max(rock_fraction, 0.0)), ("metal", metal_fraction))}
    solid_volume = volumes["rock"] + volumes["metal"]
    total = solid_volume + volumes["ice"]
    skeleton = solid_volume / total
    metal_share = volumes["metal"] / solid_volume if solid_volume > 0.0 else 0.0
    solid_strength = table.strengths_pa["rock"] ** (1.0 - metal_share) * table.strengths_pa["metal"] ** metal_share
    low, high = table.mixing["skeleton_below"], table.mixing["skeleton_above"]
    weight = min(max((skeleton - low) / (high - low), 0.0), 1.0)
    return ice_strength ** (1.0 - weight) * solid_strength ** weight


MIXING_RULES = {"volume-threshold": _volume_threshold}
"""dict: Rules for the strength of a mix of materials, each `(table, ice_fraction, metal_fraction, temperature_k)` -> Pa. The strength research adds a rule here and names it in the file. [Stand-in RN3, tier T6, tweak L3.]"""

POROSITY_RULES = {
    "none": lambda table, strength, porosity: strength,
    "exponential": lambda table, strength, porosity: strength * math.exp(-_positive(table.porosity["coefficient"], "coefficient") * porosity),
}
"""dict: Rules for what voids do to the strength, each `(table, strength_pa, porosity)` -> Pa. `none` leaves it to the density; `exponential` is `sigma exp(-b phi)` with `b` in the file (4 to 7 in the research). [Stand-in RN5, tier T3, tweak L3.]"""


def _read_table(path=DATA_FILE):
    with open(path, encoding="utf-8") as handle:
        return load_table(json.load(handle))


TABLE = _read_table()
"""StrengthTable: The strengths and rules in use (`data/roundness_strengths.json`). Rows RN1 to RN7 of the stand-in registry."""

SPHERICITY_FACTOR = TABLE.sphericity_factor
"""float: The longest axis of a body counted round is at most this times its mean radius (Lineweaver and Norman 2010, a 10 percent bump; their Eq. 9). [Stand-in RN1, tier T1, tweak L0.]"""

ICE_STRENGTH_PA = TABLE.strengths_pa["ice"]
"""float: The yield strength of ice that was warm when the shape was set, Pa (1.2 MPa), from a round radius of 200 km at 1000 kg/m3. [Stand-in RN2, tier T1, tweak L0.]"""

ROCK_STRENGTH_PA = TABLE.strengths_pa["rock"]
"""float: The yield strength of rock when the shape was set, Pa (16.5 MPa), from a round radius of 300 km at 2500 kg/m3. [Stand-in RN2, tier T1, tweak L0.]"""

METAL_STRENGTH_PA = TABLE.strengths_pa["metal"]
"""float: The yield strength of iron and nickel metal, Pa (100 MPa; the research's range is 30 to 300 MPa). [Stand-in RN4, tier T4, tweak L3.]"""


def _share(value, name):
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"the {name} is a share, 0 to 1")
    return value


@finite_domain()
def yield_strength_pa(ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None):
    """The yield strength of a body, Pa: the mixing rule of the table over its ice and metal mass fractions (ice at the strength of `temperature_k` when its shape was set; None for ice that
    was warm), then the porosity rule over its void share."""
    ice, metal, voids = _share(ice_fraction, "ice fraction"), _share(metal_fraction, "metal fraction"), _share(porosity, "porosity")
    if ice + metal > 1.0 + 1e-9:
        raise ValueError("the ice and metal fractions add up to more than the whole body")
    if temperature_k is not None and not temperature_k > 0.0:
        raise ValueError("the temperature must be above absolute zero")
    strength = MIXING_RULES[TABLE.mixing["rule"]](TABLE, ice, metal, temperature_k)
    return POROSITY_RULES[TABLE.porosity["rule"]](TABLE, strength, voids)


@finite_domain()
def round_radius_m(density_kg_m3, ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None):
    """The mean radius at which a solid body of this bulk density and composition turns round, m: `sqrt(sigma / (K G rho^2))`."""
    if not density_kg_m3 > 0.0:
        raise ValueError("the density must be positive")
    sigma = yield_strength_pa(ice_fraction, metal_fraction, porosity, temperature_k)
    return math.sqrt(sigma / (stress_coefficient(TABLE.sphericity_factor) * constants.G * density_kg_m3 ** 2))


@finite_domain()
def round_mass_kg(density_kg_m3, ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None):
    """The mass at which a solid body of this bulk density and composition turns round, kg: a sphere of `round_radius_m` at that density."""
    return 4.0 / 3.0 * math.pi * round_radius_m(density_kg_m3, ice_fraction, metal_fraction, porosity, temperature_k) ** 3 * density_kg_m3


@finite_domain()
def bulk_density_kg_m3(mass_kg, radius_m):
    """The mean density of a body of this mass and mean radius, kg/m3."""
    if not (mass_kg > 0.0 and radius_m > 0.0):
        raise ValueError("the mass and the radius must be positive")
    return mass_kg / (4.0 / 3.0 * math.pi * radius_m ** 3)


@finite_domain()
def roundness_ratio(radius_m, density_kg_m3, ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None):
    """The body's mean radius over its round radius: 1 or more is round (`is_round`), less is not, and the number says how far."""
    if not radius_m > 0.0:
        raise ValueError("the radius must be positive")
    return radius_m / round_radius_m(density_kg_m3, ice_fraction, metal_fraction, porosity, temperature_k)


@finite_domain()
def roundness_ratio_from_mass(mass_kg, density_kg_m3, ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None):
    """`roundness_ratio` of a body known by its mass and bulk density (the radius of the sphere they make)."""
    if not (mass_kg > 0.0 and density_kg_m3 > 0.0):
        raise ValueError("the mass and the density must be positive")
    radius_m = (3.0 * mass_kg / (4.0 * math.pi * density_kg_m3)) ** (1.0 / 3.0)
    return roundness_ratio(radius_m, density_kg_m3, ice_fraction, metal_fraction, porosity, temperature_k)


def is_round(ratio):
    """Whether a roundness ratio is a round body."""
    return ratio >= ROUND_RATIO
