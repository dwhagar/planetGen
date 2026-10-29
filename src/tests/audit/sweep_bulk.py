# tests/audit/sweep_bulk.py

"""
The rule set of TEST.210 on the body sweep (`tests.audit.sweep`): the bulk composition of a planet or moon agrees with its mass, radius,
density, gravity and the words that describe it.

Rules (names as `RULES` lists them):

- `fractions-sum`: the five whole-body shares (metal, rock, water, ice, envelope) are not negative and add up to 1.
- `detail-ranges`: the rock's FeO and Al2O3, the envelope's helium and heavy-element shares and the porosity are shares; the envelope
  details exist exactly where there is an envelope and its helium and heavy elements fit in it.
- `porosity-range`: a porosity is a void share of the volume, at least 0 and under 1.
- `density-from-mass-radius`: the stored density is the mass over the volume of the stored radius.
- `gravity-from-mass-radius`: the stored surface gravity is G M / R^2 in Earth gravities.
- `radius-from-structure`: the radius is the one the structure model (`bulk_structure.structure`) gives the stored composition and mass:
  to 1e-6 for a planet, to 4 percent for a moon (a moon's mass follows the tapered brittle law of GEN.442 and its grain density, not the
  structure solver). A body with an envelope is left out: its radius depends on the flux and age it was drawn at, which are not stored.
- `denser-than-grains`: with no envelope the density is at least the zero-pressure mixing density of the stored shares and porosity (the
  pressure of the body's own weight only squeezes), less 2 percent for the two ways a moon's porosity is found.
- `text-materials`: the stored phrase is the writer's phrase for the stored shares, names exactly the materials that hold at least
  `composition_text.MIN_SHARE` of the mass and says "loosely packed" exactly where the porosity reaches `composition_text.POROUS`.
"""

import math

from planetgen.physics import bulk_composition, bulk_structure, composition_text, constants
from tests.audit.sweep import Rule, broken_copy

PLANET_RADIUS_REL = 1e-6
MOON_RADIUS_REL = 0.04
GRAIN_SLACK = 0.02
REL = 1e-9


def _bulk(body):
    """`(BulkFractions, [])` for a body whose stored shares make one, `(None, [reason])` for a body that has none or a bad one."""
    try:
        bulk = bulk_composition.bulk_fractions(body)
    except bulk_composition.CompositionError as error:
        return None, [str(error)]
    return bulk, ([] if bulk is not None else ["no stored composition"])


def _shares(bulk):
    return bulk.shares()


# --- The rules ---------------------------------------------------------------------------------------------------------------


def _check_fractions_sum(entry):
    return _bulk(entry.body)[1]


def _break_fractions_sum(entry):
    return broken_copy(entry.body, metal_mass_fraction=(entry.body.metal_mass_fraction or 0.0) + 0.2)


def _check_detail_ranges(entry):
    bulk, problems = _bulk(entry.body)
    if bulk is None:
        return problems
    for name in ("rock_feo", "rock_al2o3", "envelope_helium", "envelope_z", "porosity"):
        value = getattr(bulk, name)
        if value is not None and not 0.0 <= value <= 1.0:
            problems.append(f"{name} is {value!r}, not a share")
    has_envelope = bulk.envelope > 0.0
    for name in ("envelope_helium", "envelope_z"):
        if (getattr(bulk, name) is not None) != has_envelope:
            problems.append(f"{name} is {getattr(bulk, name)!r} with an envelope share of {bulk.envelope!r}")
    if has_envelope and bulk.envelope_helium is not None and bulk.envelope_z is not None and bulk.envelope_helium + bulk.envelope_z > 1.0:
        problems.append(f"the envelope's helium {bulk.envelope_helium!r} and heavy elements {bulk.envelope_z!r} are more than all of it")
    if bulk.rock > 0.0 and not has_envelope and bulk.rock_feo is None:    # a giant's rock is its core, with no mantle FeO stored
        problems.append(f"a rock share of {bulk.rock!r} with no FeO")
    return problems


def _break_detail_ranges(entry):
    return broken_copy(entry.body, rock_feo_fraction=1.5) if getattr(entry.body, "rock_feo_fraction", None) is not None else None


def _check_porosity_range(entry):
    porosity = getattr(entry.body, "porosity", None)
    if porosity is not None and not 0.0 <= porosity < 1.0:
        return [f"a porosity of {porosity!r}"]
    return []


def _break_porosity_range(entry):
    return broken_copy(entry.body, porosity=1.2)


def _check_density(entry):
    body = entry.body
    volume_m3 = 4.0 / 3.0 * math.pi * (body.radius * 1000.0) ** 3
    expected = body.mass / volume_m3 / 1000.0
    if not math.isclose(body.density, expected, rel_tol=REL):
        return [f"density {body.density!r} g/cm3, mass over volume {expected!r}"]
    return []


def _break_density(entry):
    return broken_copy(entry.body, density=entry.body.density * 1.1)


def _check_gravity(entry):
    body = entry.body
    expected = constants.G * body.mass / (body.radius * 1000.0) ** 2 / constants.EARTH_GRAVITY
    if not math.isclose(body.gravity, expected, rel_tol=1e-6):
        return [f"gravity {body.gravity!r} g, G M / R^2 gives {expected!r}"]
    return []


def _break_gravity(entry):
    return broken_copy(entry.body, gravity=entry.body.gravity * 1.1)


def _check_radius(entry):
    body = entry.body
    bulk, problems = _bulk(body)
    if bulk is None:
        return problems
    if bulk.envelope > 0.0:
        return []
    x_fe = bulk_structure.x_fe_from_feo_mass_fraction(bulk.rock_feo) if bulk.rock_feo else 0.0
    radius_km = bulk_structure.structure(bulk, body.mass, x_fe=x_fe).radius_km
    allowed = MOON_RADIUS_REL if entry.is_moon else PLANET_RADIUS_REL
    if abs(radius_km - body.radius) > allowed * body.radius:
        return [f"radius {body.radius!r} km, the structure of its composition and mass gives {radius_km!r} km"]
    return []


def _break_radius(entry):
    bulk, _ = _bulk(entry.body)
    if bulk is None or bulk.envelope > 0.0:
        return None
    return broken_copy(entry.body, radius=entry.body.radius * 1.5)


def _zero_pressure_density_g_cm3(bulk):
    return bulk_structure.mixing_density_kg_m3(bulk.metal, bulk.rock, bulk.water, bulk.ice, bulk.porosity or 0.0) / 1000.0


def _check_denser_than_grains(entry):
    bulk, problems = _bulk(entry.body)
    if bulk is None:
        return problems
    if bulk.envelope > 0.0:
        return []
    floor = _zero_pressure_density_g_cm3(bulk) * (1.0 - GRAIN_SLACK)
    if entry.body.density < floor:
        return [f"density {entry.body.density!r} g/cm3, under the {floor!r} of its materials at zero pressure"]
    return []


def _break_denser_than_grains(entry):
    bulk, _ = _bulk(entry.body)
    if bulk is None or bulk.envelope > 0.0:
        return None
    return broken_copy(entry.body, density=_zero_pressure_density_g_cm3(bulk) * 0.5)


def _check_text_materials(entry):
    body = entry.body
    bulk, problems = _bulk(body)
    if bulk is None:
        return problems
    phrase = body.composition
    expected = composition_text.text(body)
    if phrase != expected:
        problems.append(f"stored phrase {phrase!r}, the writer gives {expected!r} for the stored shares")
    if not isinstance(phrase, str):
        return problems
    for fraction, share in _shares(bulk).items():
        name = composition_text.rock_name(bulk.carbon_rock) if fraction == "rock" else composition_text.NAMES[fraction]
        if (name in phrase) != (share >= composition_text.MIN_SHARE):
            problems.append(f"{name!r} {'is' if name in phrase else 'is not'} in the phrase with a share of {share!r}")
    loose = bulk.porosity is not None and bulk.porosity >= composition_text.POROUS
    if phrase.endswith("loosely packed") != loose:
        problems.append(f"porosity {bulk.porosity!r} and the phrase {'does' if not loose else 'does not'} end in 'loosely packed'")
    return problems


def _break_text_materials(entry):
    bulk, _ = _bulk(entry.body)
    if bulk is None or not isinstance(getattr(entry.body, "composition", None), str):
        return None
    return broken_copy(entry.body, composition="mostly unobtainium")


RULES = (
    Rule("fractions-sum", _check_fractions_sum, _break_fractions_sum),
    Rule("detail-ranges", _check_detail_ranges, _break_detail_ranges),
    Rule("porosity-range", _check_porosity_range, _break_porosity_range),
    Rule("density-from-mass-radius", _check_density, _break_density),
    Rule("gravity-from-mass-radius", _check_gravity, _break_gravity),
    Rule("radius-from-structure", _check_radius, _break_radius),
    Rule("denser-than-grains", _check_denser_than_grains, _break_denser_than_grains),
    Rule("text-materials", _check_text_materials, _break_text_materials),
)
