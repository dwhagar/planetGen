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
   **Hot rock.** A rock body that took its shape hot is weaker than cold rock: with a shape-set `temperature_k` the rock strength follows the creep of dry olivine (`TEMPERATURE_RULES`,
   `olivine-creep-law`), between the Hirth and Kohlstedt 2003 and Korenaga and Karato 2008 columns, and never above the cold whole-body strength of 16.5 MPa, which stands for the brittle limit. So
   rock keeps 16.5 MPa up to about 1100 to 1150 K, then falls to 2 to 5 MPa at 1300 K and 1 to 2 MPa at 1400 K. A rubble pile that took its shape hot is held by the same law (its `hot_strength_pa`
   of 2.5 MPa is olivine creep near 1330 K). Intact olivine crystal strength (10 GPa) is never used.

   **Fuzzy boundary.** Real strength is never known to better than a factor of 2 to 3, so every function takes an optional `draw` (0 to 1): the quantile of a body in the range of each strength (`RANGES_PA`, and
   the friction angle, cohesion and hot strength of rubble). The stand-in is the median (0.5; `draw=None` is the same). The caller draws it once per body from `draw.py`; the round radius then moves
   about 1.5 times either way, so the boundary between a round body and a potato is a band.

   **Rubble piles.** A pile of pieces (`structure="rubble"`) holds together by friction and a little cohesion, not by the strength of rock. Friction grows with the pressure, and the pressure at
   the centre grows as the square of the radius, just as the stress that would level a bump does, so with a friction angle over about 5.5 degrees a cold pile never turns round at any size: its
   roundness ratio stays near 0.28. It rounds only if it took its shape weak from heat (`heated`, the table's hot strength of 1 to 4 MPa) or, for an icy pile, at the creep strength of its ice at
   the `temperature_k` its shape was set at. Friction angle 35 degrees (30 to 40) and cohesion 44 to 210 Pa are from the studies of fast-spinning asteroids (1950 DA, Rozitis and others 2014;
   Hirabayashi and Scheeres). The failure stress is a Mohr-Coulomb fit (`RUBBLE_RULES`); a body's intact crystals (olivine, about 10 GPa) are never used for any of this.
3. **Binding energy** is in `binding_energy.py`.

The ice fraction is the share of the body's mass that is water or ice (`makeup.water + makeup.ice` of a small body). Only solid bodies have a roundness or a strength: a gas giant or a star
is a fluid and always settles into a ball, so nothing here applies to it.

Test (not an input): the ten bodies of the report, Mimas, Miranda, Enceladus, Proteus, Hyperion, Ceres, Vesta, Pallas and Juno fall on the right side; Hygiea (round at about 217 km with
a density of about 1.9 g/cm3, from Vernazza et al. 2019, not opened) falls on the wrong side, because the rock strength is fixed by the 300 km of the larger asteroids.

Use:

    ratio = roundness_ratio(radius_m, density_kg_m3, ice_fraction)           # 1.0 or more is round; metal_fraction, porosity and temperature_k are optional
    sigma = yield_strength_pa(ice_fraction, metal_fraction, porosity, temperature_k)
    ratio = roundness_ratio_from_mass(mass_kg, density_kg_m3, ice_fraction)
    ratio = roundness_ratio(radius_m, density_kg_m3, ice_fraction, structure="rubble", heated=False)       # a pile of rubble
    ratio = roundness_ratio(radius_m, density_kg_m3, ice_fraction, draw=u)                                # strength at quantile u of its range (0.5 is the stand-in)
"""

import json
import math
import os
import dataclasses
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
    rubble: dict
    ranges_pa: dict = None
    ice_creep_scale: float = 1.0
    draw: float = 0.5


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
    mixing, porosity, temperature, rubble = document.get("mixing") or {}, document.get("porosity") or {}, document.get("temperature") or {}, document.get("rubble") or {}
    if mixing.get("rule") not in MIXING_RULES:
        raise RoundnessError(f"the mixing rule must be one of {sorted(MIXING_RULES)}")
    if porosity.get("rule") not in POROSITY_RULES:
        raise RoundnessError(f"the porosity rule must be one of {sorted(POROSITY_RULES)}")
    if not temperature or set(temperature) - set(MATERIALS):
        raise RoundnessError(f"'temperature' holds a rule for each of {MATERIALS} it applies to")
    for name, block in temperature.items():
        if block.get("rule") not in TEMPERATURE_RULES:
            raise RoundnessError(f"the temperature rule of {name} must be one of {sorted(TEMPERATURE_RULES)}")
    if rubble.get("rule") not in RUBBLE_RULES:
        raise RoundnessError(f"the rubble rule must be one of {sorted(RUBBLE_RULES)}")
    if rubble["rule"] == "mohr-coulomb":
        if not 0.0 < _positive(rubble.get("friction_angle_deg"), "friction_angle_deg") < 90.0:
            raise RoundnessError("the friction angle lies between 0 and 90 degrees")
        _positive(rubble.get("cohesion_pa"), "cohesion_pa")
        _positive(rubble.get("hot_strength_pa"), "hot_strength_pa")
    if mixing["rule"] == "exp-volume":
        _positive(mixing.get("exponent"), "exponent")
        tested, packing = _positive(mixing.get("tested_to"), "tested_to"), _positive(mixing.get("close_packing"), "close_packing")
        if not tested <= packing < 1.0:
            raise RoundnessError("the mixing law needs 0 < tested_to <= close_packing < 1")
    for name, block in temperature.items():
        if block["rule"] in ("creep-law", "olivine-creep-law"):
            _positive(block.get("rate_per_s"), "rate_per_s")
            regimes = block.get("regimes" if block["rule"] == "creep-law" else "columns")
            if not isinstance(regimes, list) or not regimes or (block["rule"] == "olivine-creep-law" and len(regimes) != 2):
                raise RoundnessError(f"the {block['rule']} of {name} needs its list of regimes (two columns for olivine)")
            if block["rule"] == "creep-law":
                _positive(block.get("max_temperature_k"), "max_temperature_k")
            for regime in regimes:
                if not isinstance(regime.get("log10_a"), (int, float)):
                    raise RoundnessError("each creep regime needs 'log10_a'")
                _positive(regime.get("n"), "n")
                _positive(regime.get("e_kj_mol"), "e_kj_mol")
    ranges = {}
    for name in MATERIALS:
        bounds = materials[name].get("range_pa")
        if bounds is not None:
            if not (isinstance(bounds, list) and len(bounds) == 2 and 0.0 < bounds[0] <= strengths[name] <= bounds[1]):
                raise RoundnessError(f"{name} range_pa needs [low, high] around its strength")
            ranges[name] = (float(bounds[0]), float(bounds[1]))
    for key in ("friction_range_deg", "cohesion_range_pa", "hot_range_pa"):
        bounds = rubble.get(key)
        if bounds is not None and not (isinstance(bounds, list) and len(bounds) == 2 and 0.0 < bounds[0] <= bounds[1]):
            raise RoundnessError(f"rubble {key} needs [low, high]")
    return StrengthTable(strengths, factor, mixing, porosity, temperature, rubble, ranges)


_GAS_KJ_MOL_K = 8.314462618e-3


def _creep_stress_mpa(regime, rate, kelvin):
    """The stress, MPa, at which one flow-law regime `rate = A sigma^n exp(-E / R T)` runs at `rate` per second: `(rate / A)^(1/n) exp(E / (n R T))`. The pressure term is dropped."""
    try:
        return (rate / 10.0 ** regime["log10_a"]) ** (1.0 / regime["n"]) * math.exp(regime["e_kj_mol"] / (regime["n"] * _GAS_KJ_MOL_K * kelvin))
    except OverflowError:       # a few kelvin: the stress that moves it in 1e-15 per second is past any float, and nothing creeps
        return math.inf


def _ice_creep_law(table, rule, ice_strength_pa, temperature_k):
    """The creep strength of pure ice that took its shape at `temperature_k`, Pa: the stress of each regime, the weakest governing (Durham, Kirby and Stern 1997 Table 2, the pressure term
    dropped). None is the default ice, `ice_strength_pa` (1 MPa, the law near 158 K); a temperature above the melting point is held at it. Dislocation creep of fine-grained lab ice, so an upper
    bound at the lowest stresses."""
    if temperature_k is None:
        return ice_strength_pa
    kelvin = min(temperature_k, _positive(rule["max_temperature_k"], "max_temperature_k"))
    rate = _positive(rule["rate_per_s"], "rate_per_s")
    return min(_creep_stress_mpa(regime, rate, kelvin) for regime in rule["regimes"]) * 1.0e6 * table.ice_creep_scale


def _olivine_creep_law(table, rule, rock_strength_pa, temperature_k):
    """The strength of silicate rock that took its shape at `temperature_k`, Pa: dry dislocation creep of olivine, between the two compilations of Korenaga and Karato 2008 Table 3 (Hirth and
    Kohlstedt 2003 as printed there, and their re-fit) by the quantile of the body (`table.draw`; the geometric mean at the median), and never above `rock_strength_pa`, the 16.5 MPa whole-body
    strength of fractured cold rock that stands for the brittle (Byerlee) limit. So rock keeps its cold strength until it is hot enough to creep below it (about 1100 to 1150 K), then weakens:
    2 to 5 MPa at 1300 K. None is cold rock. The pressure term is dropped; intact olivine crystal strength (10 GPa) is never used."""
    if temperature_k is None:
        return rock_strength_pa
    rate = _positive(rule["rate_per_s"], "rate_per_s")
    low, high = (_creep_stress_mpa(column, rate, temperature_k) for column in rule["columns"])
    return min(rock_strength_pa, low ** (1.0 - table.draw) * high ** table.draw * 1.0e6)


TEMPERATURE_RULES = {
    "none": lambda table, rule, strength_pa, temperature_k: strength_pa,
    "creep-law": _ice_creep_law,
    "olivine-creep-law": _olivine_creep_law,
}
"""dict: Rules for the strength of a material that took its shape at a temperature, each `(table, the material's block of the file, default_strength_pa, temperature_k)` -> Pa. `temperature_k` of
None is the default strength. [Stand-in RN7 for ice, tier T5, tweak L1; RN11 for rock, tier T2, tweak L0.]"""


def _material_strength_pa(table, name, temperature_k):
    """The strength of one material at the temperature its shape was set: its default strength put through its temperature rule, if the file gives it one."""
    block = table.temperature.get(name)
    base = table.strengths_pa[name]
    return base if block is None else TEMPERATURE_RULES[block["rule"]](table, block, base, temperature_k)


def _exp_volume(table, ice_fraction, metal_fraction, temperature_k):
    """The strength of a mix by the volume shares of its materials (strength research section 7): ice strength times `exp(b f)` for the solids' volume fraction `f` up to `tested_to` (0.56, measured:
    Durham and others 1992), the same law continued to `close_packing` (0.65), then geometric to the solids' own strength (rock and metal, geometric mean weighted by the metal share of their
    volume) at `f` of 1. Grains of 1 to 150 microns in ice; the leg beyond 0.56 is an extrapolation."""
    ice_strength = _material_strength_pa(table, "ice", temperature_k)
    rock_fraction = 1.0 - ice_fraction - metal_fraction
    volumes = {name: share / GRAIN_DENSITIES_KG_M3[name] for name, share in (("ice", ice_fraction), ("rock", max(rock_fraction, 0.0)), ("metal", metal_fraction))}
    solid_volume = volumes["rock"] + volumes["metal"]
    skeleton = solid_volume / (solid_volume + volumes["ice"])
    metal_share = volumes["metal"] / solid_volume if solid_volume > 0.0 else 0.0
    solid_strength = _material_strength_pa(table, "rock", temperature_k) ** (1.0 - metal_share) * _material_strength_pa(table, "metal", temperature_k) ** metal_share
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


def _mohr_coulomb(table):
    """The strength of a friction-limited pile of rubble as `(b, a)`: the differential stress it fails at, at the centre of a body, is `b + a G rho^2 R^2` Pa. Mohr-Coulomb with a confining
    pressure of the centre's `(2 pi / 3) G rho^2 R^2`: `(2 c cos(phi) + 2 P sin(phi)) / (1 - sin(phi))` for cohesion `c` and friction angle `phi`. The cohesion is `b`, the part that does not
    grow with the body; the friction `a` grows as its square, as the bump's stress does."""
    rule = table.rubble
    phi = math.radians(_positive(rule["friction_angle_deg"], "friction_angle_deg"))
    cohesion = _positive(rule["cohesion_pa"], "cohesion_pa")
    return 2.0 * cohesion * math.cos(phi) / (1.0 - math.sin(phi)), (2.0 * math.pi / 3.0) * 2.0 * math.sin(phi) / (1.0 - math.sin(phi))


RUBBLE_RULES = {"mohr-coulomb": _mohr_coulomb}
"""dict: Rules for the strength of a pile of rubble, each `(table)` -> `(b, a)` for a failure stress of `b + a G rho^2 R^2` Pa. [Stand-in RN8 (friction angle, cohesion) and RN9 (hot strength).]"""


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


def _span(draw, low, mid, high):
    """The value at quantile `draw` (0 to 1) of a range `low` to `high` whose median is `mid`: log-linear from `low` at 0 to `mid` at 0.5 to `high` at 1."""
    if draw <= 0.5:
        return low * (mid / low) ** (draw / 0.5)
    return mid * (high / mid) ** ((draw - 0.5) / 0.5)


def _table_at(draw):
    """The `TABLE`, or a copy of it with every strength moved to one quantile `draw` of its range (None is the stand-in itself, the median). One quantile moves all materials together: a body
    that is strong for its kind is strong in every material, in the friction angle and the cohesion of its rubble, and in the hot strength. A material or number with no range stays as it is."""
    if draw is None:
        return TABLE
    if not isinstance(draw, (int, float)) or isinstance(draw, bool) or not 0.0 <= draw <= 1.0:
        raise ValueError("the strength draw is a number from 0 to 1")
    strengths = {name: _span(draw, RANGES_PA[name][0], value, RANGES_PA[name][1]) if name in RANGES_PA else value for name, value in TABLE.strengths_pa.items()}
    rubble = dict(TABLE.rubble)
    for key, bounds in (("friction_angle_deg", "friction_range_deg"), ("cohesion_pa", "cohesion_range_pa"), ("hot_strength_pa", "hot_range_pa")):
        if bounds in rubble:
            rubble[key] = _span(draw, rubble[bounds][0], rubble[key], rubble[bounds][1])
    return dataclasses.replace(TABLE, strengths_pa=strengths, rubble=rubble, ice_creep_scale=strengths["ice"] / TABLE.strengths_pa["ice"], draw=float(draw))


@finite_domain()
def drawn_strengths_pa(draw):
    """The yield strength of each material, Pa, at quantile `draw` (0 to 1) of its range: `{"ice": ..., "rock": ..., "metal": ...}`. At 0.5 they are the stand-ins."""
    return dict(_table_at(draw).strengths_pa)


ROCK_CREEP_COLUMNS = TABLE.temperature["rock"]["columns"]
"""list: The two flow-law columns for dry olivine dislocation creep (Hirth and Kohlstedt 2003; Korenaga and Karato 2008 re-fit), `log10_a`, `n` and `e_kj_mol` each, between which a body's quantile `draw` sits. [Stand-in RN11, tier T2, tweak L0.]"""

RANGES_PA = TABLE.ranges_pa
"""dict: The range each material's strength is drawn from, `(low, high)` Pa: ice 0.4 to 3 MPa, rock 4 to 30, metal 50 to 250. The quantile `draw` of a body sits log-linearly in it with the stand-in strength at 0.5. [Stand-in RN10, tier T4, tweak L0.]"""


def _share(value, name):
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"the {name} is a share, 0 to 1")
    return value


@finite_domain()
def yield_strength_pa(ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None, draw=None):
    """The yield strength of a body, Pa: the mixing rule of the table over its ice and metal mass fractions (ice at the creep strength of `temperature_k` when its shape was set; None for the default ice), then the porosity rule over its void share."""
    ice, metal, voids = _share(ice_fraction, "ice fraction"), _share(metal_fraction, "metal fraction"), _share(porosity, "porosity")
    if ice + metal > 1.0 + 1e-9:
        raise ValueError("the ice and metal fractions add up to more than the whole body")
    if temperature_k is not None and not temperature_k > 0.0:
        raise ValueError("the temperature must be above absolute zero")
    table = _table_at(draw)
    strength = MIXING_RULES[table.mixing["rule"]](table, ice, metal, temperature_k)
    return POROSITY_RULES[table.porosity["rule"]](table, strength, voids)


@finite_domain()
def round_radius_m(density_kg_m3, ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None, draw=None):
    """The mean radius at which a solid body of this bulk density and composition turns round, m: `sqrt(sigma / (K G rho^2))`."""
    if not density_kg_m3 > 0.0:
        raise ValueError("the density must be positive")
    sigma = yield_strength_pa(ice_fraction, metal_fraction, porosity, temperature_k, draw)
    return math.sqrt(sigma / (stress_coefficient(TABLE.sphericity_factor) * constants.G * density_kg_m3 ** 2))


@finite_domain()
def round_mass_kg(density_kg_m3, ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None, draw=None):
    """The mass at which a solid body of this bulk density and composition turns round, kg: a sphere of `round_radius_m` at that density."""
    return 4.0 / 3.0 * math.pi * round_radius_m(density_kg_m3, ice_fraction, metal_fraction, porosity, temperature_k, draw) ** 3 * density_kg_m3


STRUCTURES = ("monolith", "rubble")
"""tuple: How a solid body holds together: `monolith` is a coherent body with the strength of its material; `rubble` is a pile of pieces held by friction and a little cohesion."""


def _heated_strength_pa(ice_fraction, metal_fraction, porosity, temperature_k, heated, draw=None):
    """The strength a rubble pile has if its shape was set while it was weak from heat, Pa, or None: `heated` is a rubble pile that was hot or partly molten when it took its shape (the table's
    `hot_strength_pa`); an ice-bearing pile with a shape-set `temperature_k` creeps at the strength of its ice and rock mix at that temperature. The weakest governs."""
    found = []
    if heated:
        found.append(_positive(_table_at(draw).rubble["hot_strength_pa"], "hot_strength_pa"))
    if temperature_k is not None:
        at_temperature = yield_strength_pa(ice_fraction, metal_fraction, porosity, temperature_k, draw)
        if ice_fraction > 0.0 or at_temperature < yield_strength_pa(ice_fraction, metal_fraction, porosity, None, draw):      # ice creeps; rock only once hot enough to creep below its cold strength
            found.append(at_temperature)
    return min(found) if found else None


@finite_domain()
def rubble_strength_pa(radius_m, density_kg_m3, ice_fraction=0.0, metal_fraction=0.0, porosity=0.0, temperature_k=None, heated=False, draw=None):
    """The differential stress, Pa, that a pile of rubble of this size and density fails at: the friction-limited stress at its centre (`RUBBLE_RULES`), or the strength it has from the heat its
    shape was set in (`heated`, `temperature_k` of an ice-bearing pile), whichever is smaller. Friction grows with the pressure, so as the body grows its strength grows as fast as the stress that
    would level a bump: a cold pile does not become round by being large."""
    if not (radius_m > 0.0 and density_kg_m3 > 0.0):
        raise ValueError("the radius and the density must be positive")
    ice, metal, voids = _share(ice_fraction, "ice fraction"), _share(metal_fraction, "metal fraction"), _share(porosity, "porosity")
    table = _table_at(draw)
    cohesive, friction = RUBBLE_RULES[table.rubble["rule"]](table)
    failure = cohesive + friction * constants.G * density_kg_m3 ** 2 * radius_m ** 2
    hot = _heated_strength_pa(ice, metal, voids, temperature_k, heated, draw)
    return failure if hot is None else min(failure, hot)


@finite_domain()
def rubble_round_radius_m(density_kg_m3, ice_fraction=0.0, metal_fraction=0.0, porosity=0.0, temperature_k=None, heated=False, draw=None):
    """The mean radius at which a pile of rubble turns round, m: where the stress of a bump, `K G rho^2 R^2`, reaches the smaller of the heat-set strength and the friction-limited failure stress.
    Refused (`ValueError`) for a cold pile, which never rounds with a friction angle above about 5.5 degrees at any size."""
    if not density_kg_m3 > 0.0:
        raise ValueError("the density must be positive")
    ice, metal, voids = _share(ice_fraction, "ice fraction"), _share(metal_fraction, "metal fraction"), _share(porosity, "porosity")
    bump = stress_coefficient(TABLE.sphericity_factor) * constants.G * density_kg_m3 ** 2
    table = _table_at(draw)
    cohesive, friction = RUBBLE_RULES[table.rubble["rule"]](table)
    reaches = []
    hot = _heated_strength_pa(ice, metal, voids, temperature_k, heated, draw)
    if hot is not None:
        reaches.append(math.sqrt(hot / bump))
    if bump > friction * constants.G * density_kg_m3 ** 2:
        reaches.append(math.sqrt(cohesive / (bump - friction * constants.G * density_kg_m3 ** 2)))
    if not reaches:
        raise ValueError("a friction-limited pile of rubble with no heat-set strength never turns round")
    return min(reaches)


@finite_domain()
def bulk_density_kg_m3(mass_kg, radius_m):
    """The mean density of a body of this mass and mean radius, kg/m3."""
    if not (mass_kg > 0.0 and radius_m > 0.0):
        raise ValueError("the mass and the radius must be positive")
    return mass_kg / (4.0 / 3.0 * math.pi * radius_m ** 3)


@finite_domain()
def roundness_ratio(radius_m, density_kg_m3, ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None, structure="monolith", heated=False, draw=None):
    """The body's mean radius over its round radius: 1 or more is round (`is_round`), less is not, and the number says how far. For a `rubble` pile the ratio is the square root of the bump's
    stress over the pile's failure stress (`rubble_strength_pa`), which is the same number in the same sense; it stays under 1 for a cold pile at any size."""
    if not radius_m > 0.0:
        raise ValueError("the radius must be positive")
    if structure not in STRUCTURES:
        raise ValueError(f"unknown structure {structure!r}")
    if structure == "rubble":
        bump = stress_coefficient(TABLE.sphericity_factor) * constants.G * density_kg_m3 ** 2 * radius_m ** 2
        return math.sqrt(bump / rubble_strength_pa(radius_m, density_kg_m3, ice_fraction, metal_fraction, porosity, temperature_k, heated, draw))
    return radius_m / round_radius_m(density_kg_m3, ice_fraction, metal_fraction, porosity, temperature_k, draw)


@finite_domain()
def roundness_ratio_from_mass(mass_kg, density_kg_m3, ice_fraction, metal_fraction=0.0, porosity=0.0, temperature_k=None, structure="monolith", heated=False, draw=None):
    """`roundness_ratio` of a body known by its mass and bulk density (the radius of the sphere they make)."""
    if not (mass_kg > 0.0 and density_kg_m3 > 0.0):
        raise ValueError("the mass and the density must be positive")
    radius_m = (3.0 * mass_kg / (4.0 * math.pi * density_kg_m3)) ** (1.0 / 3.0)
    return roundness_ratio(radius_m, density_kg_m3, ice_fraction, metal_fraction, porosity, temperature_k, structure, heated, draw)


def is_round(ratio):
    """Whether a roundness ratio is a round body."""
    return ratio >= ROUND_RATIO
