# planetgen/generation/collision_field.py

"""
The asteroid field two colliding rogue planets leave (GEN.1103, docs/design/systems/interstellar-asteroid-fields.md section 2, source S1)
==========================================================================================================================================

Two rogue planets that meet in the rubble window of their speeds (`rogue_collision.rubble_speed_window_ms`) leave a largest remnant of a half down to a tenth of the mass and a cloud of rubble that flies
free: an asteroid field. The rate is physical and tiny, about `1e-16` standing fields per sector, so a field of this kind is a field that stands somewhere in the Galaxy once in a long while and in
almost no sector. Boss decided (4:54 PM EDT 2026-10-30) that fields come from physical sources only, so the sector roll's density of the kind (`tuning.PHENOMENON_DENSITY_PC3["asteroid-field"]`) is the
density of these (`standing_density_pc3`), and a field a sector rolls is one of them (`collided_field`); the fields of supernova remnants are hosted (`field_host`), a hand-made field is the admin's.

    standing_density_pc3()                       # collision fields standing per pc3
    field = collided_field(config, "sector-key")  # an AsteroidField, or None where the draw left no field

The draw: a pair of terrestrial rogues (`CHARACTERISTIC_PAIR_EARTH`, with a core and an ice shell), a speed in the rubble window, an impact parameter among the hits that are not grazing (`b_crit`,
the `b = sin(theta)` law `rogue_collision.impact_angle_cdf`), the outcome of `field_composition.collision_field` (the layers of the field, the remnant) and the ledger of its bodies
(`field_mass.collision_field`). The field takes the makeup of its layers (`field_composition.field_makeup`), the speed of its pieces (`rogue_collision.unbound_rms_speed_ms`) and the history sentence of
the pair, the shock of the impact (`peak_shock_pressure_pa`) and the second remnant (`second_remnant_fraction`). Stand-ins CF1 to CF3 of `docs/design/research/gap-filling-plan.md` section 6.211.
"""

import math

from planetgen import tuning
from planetgen.generation.phenomena.asteroid_field import AsteroidField, asteroid_field_class, radius_ly_at
from planetgen.physics import (field_composition, field_dispersal, field_mass, field_rates, rogue_collision, small_body_composition,
                               unit_registry)
from planetgen.util import draw

CHARACTERISTIC_PAIR_EARTH = (1.0, 1.0)
"""tuple: The masses, Earth masses, of the pair whose collisions are counted: two Earths, the design's pair (the rate for a pair of the rogues' mass spectrum is within a factor of a few of it). [Stand-in CF1,
tier T6, tweak L1.]"""
ROGUE_SPEED_MS = 40.0e3
"""float: The characteristic relative speed of two rogues in the thin disk, m/s (30 to 60 km/s). [Stand-in CF2, tier T6, tweak L1.]"""
ROGUE_DISPERSION_MS = 25.0e3
"""float: The one-dimensional speed dispersion of the rogues, m/s (the old disk's stars'). [Stand-in CF2, tier T6, tweak L1.]"""
TERRESTRIAL_SHARE = 0.78 / 0.91
"""float: The share of the rogue planets that are terrestrial (0.78 of 0.91 per pc3, the design's figures, from the microlensing mass function). [Stand-in CF3, tier T3, tweak L1.]"""
ICE_FRACTION = 0.1
"""float: The share of a terrestrial rogue's mass in an ice shell. [Stand-in CF3, tier T9, tweak L0.]"""
LARGEST_FRACTION_TYPICAL = 0.3
"""float: The largest-remnant share of a typical collision of the rubble window (a half down to a tenth), which the lifetime of its field is worked out at. [Stand-in CF3, tier T9, tweak L0.]"""

MAX_ATTEMPTS = 100
"""int: How many draws `build_asteroid_field` makes for a rolled collision field before it gives up (a draw leaves a field about three times in five)."""

_earth_kg = unit_registry.converter("Mearth", "kg")


def _pair():
    """`(m1_kg, m2_kg, r1_m, r2_m)` of the characteristic pair."""
    masses = [_earth_kg(mass) for mass in CHARACTERISTIC_PAIR_EARTH]
    radii = [rogue_collision.combined_radius_m(mass) for mass in masses]
    return masses[0], masses[1], radii[0], radii[1]


def rogue_density_pc3():
    """The number density of terrestrial rogue planets, per pc3: the rogues per star of the `Normal` level (`tuning.ROGUE_LEVEL_CEILINGS`) at the stellar density of the solar neighborhood, times the terrestrial share."""
    return tuning.ROGUE_LEVEL_CEILINGS["Normal"] * tuning.REFERENCE_STELLAR_DENSITY_PC3 * TERRESTRIAL_SHARE


def non_grazing_share():
    """The share of the hits that are not grazing, `b_crit^2` (`rogue_collision.impact_angle_cdf` at the critical impact parameter): the rest are hit-and-runs that leave no field."""
    _m1, _m2, r1, r2 = _pair()
    return rogue_collision.impact_angle_cdf(rogue_collision.critical_impact_parameter(r1, r2))


def lifetime_yr():
    """How long a collision field stands: the time its pieces, at the typical rms speed, take to pass the retirement radius (`field_dispersal.retirement_age_yr`)."""
    m1, m2, _r1, _r2 = _pair()
    return field_dispersal.retirement_age_yr(rogue_collision.unbound_rms_speed_ms(LARGEST_FRACTION_TYPICAL, m1 + m2))


def standing_per_sector(sector_volume_pc3):
    """The collision fields that stand in a sector of `sector_volume_pc3`: the births per year (`field_rates.collision_field_rate_per_sector_per_year`, the hits that are not grazing) times the life of a field."""
    m1, m2, r1, r2 = _pair()
    rate = field_rates.collision_field_rate_per_sector_per_year(rogue_density_pc3(), sector_volume_pc3, m1, m2, r1, r2, ROGUE_SPEED_MS, ROGUE_DISPERSION_MS) * non_grazing_share()
    return field_rates.standing_number(rate, lifetime_yr())


def standing_density_pc3():
    """The collision fields standing per pc3 (`standing_per_sector` of one cubic parsec): the density the sector roll uses for the kind."""
    return standing_per_sector(1.0)


def _impact_parameter(unit, critical):
    """The impact parameter of a hit that is not grazing for `unit` in [0, 1): the law `b^2` conditioned on `b < critical`, so `critical * sqrt(unit)`."""
    return critical * math.sqrt(unit)


def collided_field(config, key):
    """
    The `AsteroidField` a collision of two terrestrial rogues leaves, drawn from the stream of `key`, or `None` where the draw leaves no field (the speed and angle gave a merger, or the unbound mass is
    under `field_mass.MIN_FIELD_MASS_EARTH`). The same key gives the same field.

    Args:
        config (SystemConfig): The shared settings the field is built with.
        key (str): A key unique to the field (its sector and place), which the draws are seeded from.
    """
    stream = draw.Stream(f"collision-field:{key}")
    m1, m2, r1, r2 = _pair()
    low, high = rogue_collision.rubble_speed_window_ms(max(m1, m2), min(m1, m2))
    speed = low + (high - low) * stream.random()
    impact = _impact_parameter(stream.random(), rogue_collision.critical_impact_parameter(r1, r2))
    parents = [field_composition.Parent(mass, field_composition.DEFAULT_CORE_FRACTION, ICE_FRACTION) for mass in (m1, m2)]
    outcome = field_composition.collision_field(parents[0], parents[1], speed, impact)
    if outcome.field is None:
        return None
    mass = field_mass.collision_field(m1, m2, r1, r2, speed, stream, impact)
    if mass is None:
        return None
    porosity = stream.uniform(*small_body_composition.POROSITY_RANGES[small_body_composition.ROCKY])
    with draw.bound(stream):
        field = AsteroidField(config)
    field.makeup = field_composition.field_makeup(outcome.field, porosity)
    field.composition = small_body_composition.describe(field.makeup)
    field.chondrite_group = field.makeup.detail.get("group", "none")
    field.composition_family = "collisional"
    field.mass_kg, field.ledger = mass.mass_kg, mass.ledger
    field.rms_speed_kms, field.shape, field.origin = mass.rms_speed_ms / 1000.0, "ball", "rogue collision"
    field.birth_age_yr = field_dispersal.age_from_unit_yr(stream.random(), mass.rms_speed_ms)
    second = rogue_collision.second_remnant_fraction(outcome.largest_fraction) * (m1 + m2)
    shock_gpa = rogue_collision.peak_shock_pressure_pa(speed) / 1.0e9
    field.history = f"{mass.sentence} The impact's peak shock was about {shock_gpa:.3g} GPa, and a second remnant of {unit_registry.convert(second, 'kg', 'Mearth'):.3g} Earth masses may stand beside the first."
    field.largest_body_kg = float(mass.ledger.masses_kg[0]) if len(mass.ledger.masses_kg) else None
    field.radius_ly = radius_ly_at(field.birth_age_yr, field.rms_speed_kms, field.shape)
    field.field_class = asteroid_field_class(field.composition_family, field.density, field.radius_ly)
    return field
