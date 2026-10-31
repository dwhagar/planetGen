# planetgen/physics/field_mass.py

"""
How much mass an interstellar asteroid field holds, and its stored bodies (GEN.1161; docs/design/systems/interstellar-asteroid-fields.md sections 3, 4 and 6)
=============================================================================================================================================================

A field is a belt with no star (the design's short answer 6). Boss approved the rule at 2026-10-31 3:14 PM EDT: **a field holds what the disc stretch it came from has left, by
the belts' own rule, and its stored bodies add up to that mass**. No field is given a typical mass from a table and none is invented:

1. **A belt unbound by a supernova** (`unbound_belt_field`, source S2, the best physical source). The belt is the leftover of the progenitor's disc stretch after the giants acted on
   it (`belt_origin`, its `mass_earth`), ground by collisions for the progenitor's life (`belt_ledger.solve_belt`: the collisional mass of `belt_mass.belt_mass_earth` and its largest
   body settled together, the belts' own step). The explosion unbinds it (`field_dispersal.unbound_belt_speed_ms`) and from then on nothing grinds it (`field_dispersal`: the
   optical depth falls under 1 within years), so the field's mass is the mass at the explosion and stays that.
2. **Two rogue planets that collided** (`collision_field`, source S1). The disc stretch is the pair themselves: the field holds the part of the combined mass that leaves the largest
   remnant faster than its escape speed (`rogue_collision.unbound_fraction`), so a hit-and-run leaves none, a merger or an erosion a trace, and a super-catastrophic hit under a seventh of the mass. The
   remnant is a rogue planet of its own (`largest_remnant_kg`), never a member of the field.

Either way the stored bodies are the belts' ledger (`belt_ledger.draw_ledger`: every body above `L1_MASS_KG`, more of the next heaviest until `LEDGER_COVERAGE` of the mass is held,
at most `LEDGER_MAX_BODIES`) drawn from the field's own stream. Ranks are stable, so `<field ID>.<rank hex>` (`galaxy.object_uid.designation`) names one body whatever the field's age.

A mass under `MIN_FIELD_MASS_EARTH` is not a field (`is_field`): the line is the belts' (GEN.1156, Boss 2026-10-31) so one rule says what is too thin to count.

The trojan swarms and the captured bodies of a planet are fields of a star's system and take their mass from their own makers (`scatter`, `irregular_population`); they are not made here.

What this does not have: the shock heating of the blast that strips the ice of a belt (`interstellar-asteroid-fields.md` section 4, S2: ice lost, rock kept, [R]) is the host rule's work (GEN.1104),
which passes the belt's icy or rocky size law; the size law of the collision's rubble is the rocky law of a main belt, a stand-in until the largest-fragment reading (LS12) replaces it (FM2).
Nothing calls the module yet: the field columns, the host rule and the member rows (GEN.1103, GEN.1104, GEN.1105) wire it.

Use:

    field = unbound_belt_field(primordial_earth, lower_au, upper_au, progenitor_solar, progenitor_age_yr, remnant_kg, belt_ledger.ROCKY_LAW, stream)
    field = collision_field(mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms, stream)      # None where nothing is unbound
"""

import math
from dataclasses import dataclass

from planetgen.physics import belt_ledger, constants, field_dispersal, rogue_collision

MIN_FIELD_MASS_EARTH = 1.0e-9
"""float: The least mass, Earth masses, a field can hold and be kept: the belts' line (`generation.belt_birth.MIN_PRESENT_MASS_EARTH`, GEN.1156), warm dust a survey could see. [Stand-in FM1, tier T9, tweak L3.]"""

COLLISION_LAW = belt_ledger.ROCKY_LAW
"""SizeLaw: The size law of the rubble of a collision: the rocky main-belt law, until a largest-fragment reading replaces it. [Stand-in FM2, tier T9, tweak L4.]"""

ORIGINS = ("unbound-belt", "collision")
"""tuple: What made a field (`FieldMass.origin`)."""


class FieldMassError(ValueError):
    """A mass, radius, speed or age argument is out of range."""


@dataclass(frozen=True)
class FieldMass:
    """
    What a field holds.

    Attributes:
        origin (str): One of `ORIGINS`.
        mass_kg (float): The field's mass, kg; the stored bodies and `ledger.remainder_kg` add up to it.
        ledger (belt_ledger.Ledger): The stored bodies, heaviest first (rank 0).
        rms_speed_ms (float): The root mean square speed of the bodies at infinity, m/s (the field's dispersion, `field_dispersal`).
        remnant_kg (float): The mass of the planet a collision leaves behind (0 for a belt).
        sentence (str): One plain sentence of where the mass came from.
    """

    origin: str
    mass_kg: float
    ledger: "belt_ledger.Ledger"
    rms_speed_ms: float
    remnant_kg: float
    sentence: str

    @property
    def mass_earth(self):
        """The field's mass, Earth masses."""
        return self.mass_kg / constants.EARTH_MASS_TO_KG


def is_field(mass_kg):
    """Whether this mass is a field: `MIN_FIELD_MASS_EARTH` or more."""
    return mass_kg >= MIN_FIELD_MASS_EARTH * constants.EARTH_MASS_TO_KG


def _positive(**values):
    for name, value in values.items():
        if not value > 0.0:
            raise FieldMassError(f"{name} {value!r} is not positive")


def unbound_belt_field(primordial_earth, lower_au, upper_au, progenitor_solar, progenitor_age_yr, remnant_kg, law, stream):
    """
    The field a supernova leaves of its progenitor's belt, or `None` where the remnant keeps over half the mass (the belt stays bound) or the belt is under `MIN_FIELD_MASS_EARTH`.

    Args:
        primordial_earth (float): The belt's mass after the giants acted on it, Earth masses (`belt_origin.BeltOrigin.mass_earth`).
        lower_au, upper_au (float): The belt's limits, AU.
        progenitor_solar (float): The progenitor's mass, solar masses (what the belt orbited).
        progenitor_age_yr (float): The progenitor's age at the explosion, years: the belt was ground for this long and no longer.
        remnant_kg (float): The remnant's mass, kg.
        law (belt_ledger.SizeLaw): `ICY_LAW` or `ROCKY_LAW` (the host rule passes the one the blast leaves).
        stream: An object with `random()` (a `draw.Stream`) the ledger is drawn from.
    """
    _positive(primordial_earth=primordial_earth, lower_au=lower_au, progenitor_solar=progenitor_solar, remnant_kg=remnant_kg)
    if not upper_au > lower_au or progenitor_age_yr < 0.0:
        raise FieldMassError("the belt needs 0 < lower < upper and the age must not be negative")
    progenitor_kg = progenitor_solar * constants.SOLAR_MASS_TO_KG
    belt_radius_m = 0.5 * (lower_au + upper_au) * constants.AU_M
    speed = field_dispersal.unbound_belt_speed_ms(progenitor_kg, remnant_kg, belt_radius_m)
    if speed <= 0.0:
        return None
    settled = belt_ledger.solve_belt(primordial_earth, lower_au, upper_au, progenitor_solar, progenitor_age_yr / 1.0e6, law, stream)
    mass_kg = settled.ledger.belt_mass_kg
    if not is_field(mass_kg):
        return None
    sentence = (f"A belt of {primordial_earth:.3g} Earth masses between {lower_au:.3g} and {upper_au:.3g} AU, ground by collisions for {progenitor_age_yr / 1.0e6:.3g} million years "
                f"to {settled.mass_earth:.3g} Earth masses, was unbound when its star exploded and leaves at {speed / 1000.0:.3g} km/s.")
    return FieldMass("unbound-belt", mass_kg, settled.ledger, speed, 0.0, sentence)


def collision_field(mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms, stream, impact_parameter=0.0):
    """
    The field two colliding rogue planets leave, or `None` where nothing is unbound (a hit-and-run) or too little to be a field or the unbound mass is under `MIN_FIELD_MASS_EARTH`.

    Args:
        mass1_kg, mass2_kg (float): The two bodies' masses, kg.
        radius1_m, radius2_m (float): Their radii, m.
        speed_ms (float): The impact speed at touching, m/s.
        stream: An object with `random()` the ledger is drawn from.
        impact_parameter (float): `b = sin(theta)` from 0 (head-on) to 1 (tangent); over `rogue_collision.critical_impact_parameter` the hit is a graze.
    """
    _positive(mass1_kg=mass1_kg, mass2_kg=mass2_kg, radius1_m=radius1_m, radius2_m=radius2_m, speed_ms=speed_ms)
    target, projectile = (mass1_kg, mass2_kg) if mass1_kg >= mass2_kg else (mass2_kg, mass1_kg)
    target_r, projectile_r = (radius1_m, radius2_m) if mass1_kg >= mass2_kg else (radius2_m, radius1_m)
    ratio = rogue_collision.specific_impact_energy_j_kg(target, projectile, speed_ms) / rogue_collision.oblique_disruption_energy_j_kg(target, projectile, target_r, projectile_r, impact_parameter)
    largest = rogue_collision.largest_remnant_fraction(ratio)
    if largest >= 1.0:
        return None
    total = target + projectile
    mass_kg = total * rogue_collision.unbound_fraction(largest)
    if not is_field(mass_kg):
        return None
    rms = rogue_collision.unbound_rms_speed_ms(largest, total)
    regime = rogue_collision.regime(ratio, impact_parameter > rogue_collision.critical_impact_parameter(target_r, projectile_r))
    ledger = belt_ledger.draw_ledger(COLLISION_LAW, mass_kg, stream)
    sentence = (f"Two rogue planets of {target / constants.EARTH_MASS_TO_KG:.3g} and {projectile / constants.EARTH_MASS_TO_KG:.3g} Earth masses met at {speed_ms / 1000.0:.3g} km/s ({regime.replace('_', ' ')}); "
                f"{mass_kg / constants.EARTH_MASS_TO_KG:.3g} Earth masses of rubble flew free and a {largest * total / constants.EARTH_MASS_TO_KG:.3g} Earth-mass remnant stayed.")
    return FieldMass("collision", mass_kg, ledger, rms, largest * total, sentence)


def stored_bodies(field):
    """`(rank, mass_kg)` of each stored body of the field, rank from 1 (the heaviest): the rank of `object_uid.designation`."""
    return [(i + 1, mass) for i, mass in enumerate(field.ledger.masses_kg)]


def mass_is_conserved(field, tolerance=1.0e-9):
    """Whether the stored bodies and the remainder add up to the field's mass (the ledger's budget rule), to a relative `tolerance`."""
    return math.isclose(sum(field.ledger.masses_kg) + field.ledger.remainder_kg, field.mass_kg, rel_tol=tolerance)
