# planetgen/generation/field_host.py

"""
The asteroid field a young supernova remnant hosts (GEN.1104, docs/design/systems/interstellar-asteroid-fields.md section 4, source S2)
=========================================================================================================================================

A star of 8 solar masses or more that explodes loses most of its mass at once, and a belt it kept (an A or B star often keeps one) is no longer bound: it leaves the place of the explosion as an
expanding ring of rubble at about 11 km/s. For the years it takes to pass a light year, the remnant of the explosion hosts that field. `hosted_field` draws it for a core-collapse remnant, from a
stream of the remnant's own key, so the same remnant has the same field:

1. The progenitor: its initial mass from the IMF over the mass windows of the kind of core the explosion left (`remnant_progenitors.sample_mass_sol`), the mass of that core (`remnant_mass`, or the
   core the remnant carries), and the age the star had at the explosion (`stellar_evolution.main_sequence_lifetime_gy`).
2. The rule (`field_rates.hosted_belt`): the star kept a belt with the share of its age (`field_rates.belt_share`), the core is under half of its mass (the belt is unbound) and the remnant is younger
   than the belt needs to pass a light year (`field_dispersal.is_retired`).
3. The mass: a belt between `BELT_LIMITS_AU` with the solids of the disk there (`belt_mass.primordial_mass_earth`), ground for the progenitor's life (`field_mass.unbound_belt_field` with the
   rocky size law: the blast strips the ice and keeps the rock) and under `field_mass.MIN_FIELD_MASS_EARTH` no field.
4. The field: an `AsteroidField` whose columns are the rule's (`rms_speed_kms`, `birth_age_yr` the remnant's age, `shape` annulus, `origin` unbound belt, `history` the sentence of where the mass
   came from, `largest_body_kg` the ledger's heaviest body) and whose radius follows from them (`asteroid_field.radius_ly_at`).

A Type Ia remnant leaves no core and no belt: it hosts nothing. Stand-ins are HF1 and HF2 of `docs/design/research/gap-filling-plan.md` section 6.202.

    hosted = hosted_field(remnant, config, "sector-key:3")      # a HostedField, or None
"""

import collections

from planetgen import tuning
from planetgen.generation.phenomena.asteroid_field import RING_RADIUS_AU, AsteroidField, asteroid_field_class, radius_ly_at
from planetgen.physics import belt_ledger, belt_mass, constants, field_mass, field_rates, formation, remnant_mass, remnant_progenitors, stellar_evolution, unit_registry
from planetgen.util import draw

BELT_LIMITS_AU = (0.6 * RING_RADIUS_AU, 1.4 * RING_RADIUS_AU)
"""tuple: The stretch the belt a progenitor kept lies in, AU: 30 to 70 about the ring (`asteroid_field.RING_RADIUS_AU`, the design's 50 AU belt of an A to B star; the width of a Kuiper-type belt). [Stand-in HF1, tier T6, tweak level L1.]"""
SHAPE = "annulus"
"""str: What a hosted field was made as: a ring at the radius of the progenitor's belt (`BELT_LIMITS_AU` middle)."""
ORIGIN = "unbound belt"
"""str: What a hosted field's origin is: the belt a supernova's progenitor kept, unbound by the mass the explosion removed."""
NEUTRON_STAR_SHARE_OF_CORES = 1.0 - tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE
"""float: The share of the cores of core-collapse remnants that are neutron stars, for a remnant that shows no core (the share of the visible ones: `tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE`). [Stand-in HF2, tier T9, tweak level L1.]"""

HostedField = collections.namedtuple("HostedField", ("field", "mass", "speed_ms", "belt_radius_m"))
"""namedtuple: A hosted field: the `AsteroidField`, its `field_mass.FieldMass`, the speed the belt left at and the radius it started from."""

_msun_to_kg = unit_registry.converter("Msun", "kg")


def belt_radius_m():
    """The radius the progenitor's belt orbited at, m: the middle of `BELT_LIMITS_AU`."""
    return 0.5 * sum(BELT_LIMITS_AU) * constants.AU_M


def _core(remnant, stream):
    """The kind of the core a core-collapse remnant has and its mass, kg: its own when it shows one, else a neutron star or a black hole by `NEUTRON_STAR_SHARE_OF_CORES`."""
    core = getattr(remnant, "compact_remnant", None)
    if core is not None:
        kind = "black-hole" if hasattr(core, "event_horizon_radius_km") else "neutron-star"
        return kind, core.mass
    kind = "neutron-star" if stream.random() < NEUTRON_STAR_SHARE_OF_CORES else "black-hole"
    return kind, _msun_to_kg(remnant_mass.draw_mass_sol(kind, None, stream))


def hosted_field(remnant, config, key):
    """
    The `HostedField` of a supernova remnant, or `None`.

    Args:
        remnant (SupernovaRemnant): The remnant (`progenitor_type`, `age_years`, the core it carries).
        config (SystemConfig): The shared settings the field is built with.
        key (str): A key unique to the remnant (its sector and place), which the draws are seeded from.
    """
    if getattr(remnant, "progenitor_type", None) != "core-collapse":
        return None
    stream = draw.Stream(f"hosted-field:{key}")
    kind, core_kg = _core(remnant, stream)
    progenitor_sol = remnant_progenitors.sample_mass_sol(kind, rng=stream)
    if progenitor_sol is None:
        return None
    progenitor_kg = _msun_to_kg(progenitor_sol)
    progenitor_age_yr = stellar_evolution.main_sequence_lifetime_gy(progenitor_sol) * 1.0e9
    radius = belt_radius_m()
    belt = field_rates.hosted_belt(stream.random(), progenitor_age_yr, progenitor_kg, core_kg, radius, remnant.age_years)
    if belt is None:
        return None
    luminosity_w = stellar_evolution.map_luminosity_sol(progenitor_sol) * constants.SOLAR_LUMINOSITY
    primordial = belt_mass.primordial_mass_earth(*BELT_LIMITS_AU, progenitor_kg, formation.snow_line_au(luminosity_w))
    mass = field_mass.unbound_belt_field(primordial, *BELT_LIMITS_AU, progenitor_sol, progenitor_age_yr, core_kg, belt_ledger.ROCKY_LAW, stream)
    if mass is None:
        return None
    with draw.bound(stream):
        field = AsteroidField(config)
    field.mass_kg, field.ledger = mass.mass_kg, mass.ledger
    field.rms_speed_kms, field.birth_age_yr, field.shape = belt.speed_ms / 1000.0, belt.age_yr, SHAPE
    field.origin, field.history, field.largest_body_kg = ORIGIN, mass.sentence, mass.ledger.masses_kg[0] if len(mass.ledger.masses_kg) else None
    field.radius_ly = radius_ly_at(field.birth_age_yr, field.rms_speed_kms, field.shape)
    field.field_class = asteroid_field_class(field.composition_family, field.density, field.radius_ly)
    return HostedField(field, mass, belt.speed_ms, radius)
