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
   strength of rock is fixed so that the round radius is exactly 300 km at 2500 kg/m3 (16.5 MPa). That is a strength during the warm period when the shape was set, not today's.

   **Where real data drops in.** The strengths are not in this module but in `data/roundness_strengths.json` (`TABLE`): one entry per material (ice, rock, metal), either a `strength_pa` or the
   round radius and density it was seen at, plus the name of a mixing rule (`MIXING_RULES`: how the ice, rock and metal shares of a body combine), of a porosity rule (`POROSITY_RULES`: what
   voids do to the strength) and of a temperature rule (`TEMPERATURE_RULES`: the strength of ice by the temperature at which the shape was set). The engine's interface is the body's ice fraction,
   metal fraction, porosity and shape-set temperature in, and the strength, the round radius and the roundness ratio out; research gives new numbers in the file, a new rule in one of the three
   dicts, and nothing else changes. The numbers in the file now are the strength research's (/mnt/project-files/reports/planetary-material-strength-2026-10-31.md, sections 7 and 11):
   rock 16.5 MPa (confirmed three ways, range 4 to 30); ice 1 MPa by default (the creep strength of pure ice near 158 K, range 0.4 to 3), and when a shape-set `temperature_k` is given the
   Durham, Kirby and Stern 1997 creep law of ice at 1e-15 per second (0.1 MPa at 200 K, 3.6 at 140 K, 34 at 100 K); metal 100 MPa (range 50 to 250, Pohl and Britt 2020); and a mix that follows
   the *volume* shares, ice strength times `exp(2 f)` for the solids' volume fraction `f` (measured to 0.56), continued to 0.65 and then geometric to the solids' strength at 1.
   Porosity is `none`: meteorites show no clear strength trend with porosity (Pohl and Britt 2020). An `exponential` rule (ceramics-based, a recalled placeholder) is built and off; with it
   Proteus (a porous ice-rich moon, 210 km, not round) would come out round. Not yet in the table: methane clathrate (20 to 30 times the strength of ice), solid CO2, ammonia ices and the
   other ice phases, and the rubble-pile branch (friction plus a little cohesion, a shape-set temperature and no size threshold).
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
    if mixing["rule"] == "exp-volume":
        _positive(mixing.get("exponent"), "exponent")
        tested, packing = _positive(mixing.get("tested_to"), "tested_to"), _positive(mixing.get("close_packing"), "close_packing")
        if not tested <= packing < 1.0:
            raise RoundnessError("the mixing law needs 0 < tested_to <= close_packing < 1")
    if temperature["rule"] == "creep-law":
        _positive(temperature.get("rate_per_s"), "rate_per_s")
        _positive(temperature.get("max_temperature_k"), "max_temperature_k")
        regimes = temperature.get("regimes")
        if not isinstance(regimes, list) or not regimes:
            raise RoundnessError("the creep law needs a list of 'regimes'")
        for regime in regimes:
            if not isinstance(regime.get("log10_a"), (int, float)):
                raise RoundnessError("each creep regime needs 'log10_a'")
            _positive(regime.get("n"), "n")
            _positive(regime.get("e_kj_mol"), "e_kj_mol")
    return StrengthTable(strengths, factor, mixing, porosity, temperature)


def _ice_creep_law(table, ice_strength_pa, temperature_k):
    """The creep strength of pure ice that took its shape at `temperature_k`, Pa: from `rate = A sigma^n exp(-E / R T)` the stress `(rate / A)^(1/n) exp(E / (n R T))` of each regime, the weakest
    governing (Durham, Kirby and Stern 1997 Table 2, the pressure term dropped). None is the default ice, `ice_strength_pa` (1 MPa, the law near 158 K); a temperature above the melting point is
    held at it. Dislocation creep of fine-grained lab ice, so an upper bound at the lowest stresses."""
    if temperature_k is None:
        return ice_strength_pa
    rule = table.temperature
    kelvin = min(temperature_k, _positive(rule["max_temperature_k"], "max_temperature_k"))
    rate = _positive(rule["rate_per_s"], "rate_per_s")
    stresses = [(rate / 10.0 ** regime["log10_a"]) ** (1.0 / regime["n"]) * math.exp(regime["e_kj_mol"] / (regime["n"] * _GAS_KJ_MOL_K * kelvin)) for regime in rule["regimes"]]
    return min(stresses) * 1.0e6


_GAS_KJ_MOL_K = 8.314462618e-3

TEMPERATURE_RULES = {
    "none": lambda table, ice_strength_pa, temperature_k: ice_strength_pa,
    "creep-law": _ice_creep_law,
}
"""dict: Rules for the strength of ice that took its shape at a temperature, each `(table, default_ice_strength_pa, temperature_k)` -> Pa. `temperature_k` of None is the default ice. [Stand-in RN7, tier T5, tweak L1.]"""


def _exp_volume(table, ice_fraction, metal_fraction, temperature_k):
    """The strength of a mix by the volume shares of its materials (strength research section 7): ice strength times `exp(b f)` for the solids' volume fraction `f` up to `tested_to` (0.56, measured:
    Durham and others 1992), the same law continued to `close_packing` (0.65), then geometric to the solids' own strength (rock and metal, geometric mean weighted by the metal share of their
    volume) at `f` of 1. Grains of 1 to 150 microns in ice; the leg beyond 0.56 is an extrapolation."""
    ice_strength = TEMPERATURE_RULES[table.temperature["rule"]](table, table.strengths_pa["ice"], temperature_k)
    rock_fraction = 1.0 - ice_fraction - metal_fraction
    volumes = {name: share / GRAIN_DENSITIES_KG_M3[name] for name, share in (("ice", ice_fraction), ("rock", max(rock_fraction, 0.0)), ("metal", metal_fraction))}
    solid_volume = volumes["rock"] + volumes["metal"]
    skeleton = solid_volume / (solid_volume + volumes["ice"])
    metal_share = volumes["metal"] / solid_volume if solid_volume > 0.0 else 0.0
    solid_strength = table.strengths_pa["rock"] ** (1.0 - metal_share) * table.strengths_pa["metal"] ** metal_share
    rule = table.mixing
    packing = rule["close_packing"]
    if skeleton <= packing:
        return ice_strength * math.exp(rule["exponent"] * skeleton)
    jammed = ice_strength * math.exp(rule["exponent"] * packing)
    weight = (skeleton - packing) / (1.0 - packing)
    return jammed ** (1.0 - weight) * solid_strength ** weight


MIXING_RULES = {"exp-volume": _exp_volume}
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
"""float: The default yield strength of ice at the scale of a body, Pa (1 MPa, the creep strength of pure ice near 158 K; the Lineweaver and Norman fit is 0.4 to 3 MPa). [Stand-in RN2, tier T2, tweak L1.]"""

ROCK_STRENGTH_PA = TABLE.strengths_pa["rock"]
"""float: The yield strength of rock when the shape was set, Pa (16.5 MPa), from a round radius of 300 km at 2500 kg/m3. [Stand-in RN2, tier T1, tweak L0.]"""

METAL_STRENGTH_PA = TABLE.strengths_pa["metal"]
"""float: The yield strength of iron and nickel metal at the scale of a body, Pa (100 MPa; range 50 to 250). [Stand-in RN4, tier T2, tweak L1.]"""


def _share(value, name):
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"the {name} is a share, 0 to 1")
    return value


@finite_domain()
def yield_strength_pa(ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None):
    """The yield strength of a body, Pa: the mixing rule of the table over its ice and metal mass fractions (ice at the creep strength of `temperature_k` when its shape was set; None for the default ice), then the porosity rule over its void share."""
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
