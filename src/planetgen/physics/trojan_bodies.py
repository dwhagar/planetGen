# planetgen/physics/trojan_bodies.py

"""
Whole bodies at L4 and L5: the pure tests behind GEN.1151 (docs/design/systems/lagrange-trojans.md, reports/lagrange-trojans-2026-10-30.md; rows 6.125 TB1 to TB4)
=====================================================================================================================================

`physics/trojans.py` gives the mass rule for any trojan. This module asks it for a planet or a moon, and says what a formation draw must supply. It sets no
rate and draws nothing; the generator (GEN.1151, after GEN.1150) calls it.

The cases, by what the host and the secondary are:

* `star_giant`: a star and a planet. The room is `max_trojan_mass(star, giant)` (about 4 percent of the star less the giant).
* `star_substellar`: a star and a brown dwarf or a very light companion under `trojans.MU_CRIT` of the pair; the room is the same sum, small after the companion.
* `star_star`: a companion over `MU_CRIT` of the pair. No room, so neither a Class M nor a Class J world (the stars of a close pair).
* `planet_moon`: a planet and a moon, in the mass of the planet. A Jupiter-mass planet (318 Earth masses) holds about 12 Earth masses of moon plus trojan.
* `moon_moon`: the same rule, a moon as the host of a smaller moon, which is `planet_moon` read with the moon as host.

A world is allowed when the room is at least its mass: `room_earth >= mass_earth`. The classes: a Class J world is 35 Earth masses or more (`physics/data/planet_classes.json`); a
Class M world is terrestrial, so it has the mass of a terrestrial world (`CLASS_M_MASS_EARTH`, a range). `classes_allowed` says which classes the room admits.

The formation gate (inputs named for the planet architecture oligarch chain, `physics/planet_architecture.py`): a trojan grows from the solids in the 1:1 zone of the
secondary, so the draw needs `solids_in_zone_earth` (the solids the oligarch chain leaves in the secondary's co-orbital feeding zone, Earth masses) and the disc lifetime
(`disc_lifetime_myr`; the trojan must grow while gas is present, [R0919]). `can_form` is true when the zone holds `FORMATION_EFFICIENCY` times
the trojan's mass in solids and the disc lasted `MIN_DISC_LIFETIME_MYR`. It is a gate, not a rate: whether a trojan appears is whatever the generated disc allows.

Use:

    room_earth(1.0 * 332946, 318.0)         # a Sun and a Jupiter, in Earth masses
    classes_allowed(332946.0, 318.0)        # {"room_earth": ..., "class_j": True, "class_m": True}
    can_form(1.0, solids_in_zone_earth=12.0, disc_lifetime_myr=3.0)
"""

from planetgen.physics import trojans

CLASS_J_MIN_EARTH = 35.0
"""float: A Class J world is this many Earth masses or more (`physics/data/planet_classes.json`, 'a gas giant: 35 Earth masses or more'). [Stand-in TB1, tier T1, tweak L0.]"""

CLASS_M_MASS_EARTH = (0.1, 10.0)
"""tuple: The masses (Earth masses) a Class M world is drawn from: terrestrial and below the first giant row (`planet_classes.json`, giants from 10 Earth masses). [Stand-in TB2, tier T6, tweak L1.]"""

FORMATION_EFFICIENCY = 0.1
"""float: The share of the solids in the secondary's co-orbital zone that grow into the trojan; the zone must hold `trojan mass / this`. [Stand-in TB3, tier T4, tweak L4.]"""

MIN_DISC_LIFETIME_MYR = 1.0
"""float: The shortest gas-disc life (Myr) in which a trojan of planet mass can grow beside its secondary. [Stand-in TB4, tier T4, tweak L4.]"""

CASES = ("star_giant", "star_substellar", "star_star", "planet_moon", "moon_moon")


def case(host_kind, secondary_kind, host_mass, secondary_mass):
    """The case name for a pair: kinds `"star"`, `"planet"` or `"moon"`; a star pair whose companion is over `trojans.MU_CRIT` of the pair is `star_star`."""
    if host_kind == "star":
        if secondary_kind == "star":
            return "star_substellar" if trojans.mass_parameter(host_mass, secondary_mass) < trojans.MU_CRIT else "star_star"
        return "star_giant"
    if host_kind == "planet" and secondary_kind == "moon":
        return "planet_moon"
    if host_kind == "moon" and secondary_kind == "moon":
        return "moon_moon"
    raise trojans.TrojanError(f"no case for a {secondary_kind} about a {host_kind}")


def room_earth(host_mass_earth, secondary_mass_earth):
    """The heaviest trojan (Earth masses) the pair holds; 0 when the pair is over the limit."""
    return trojans.max_trojan_mass(host_mass_earth, secondary_mass_earth)


def classes_allowed(host_mass_earth, secondary_mass_earth):
    """`{"room_earth", "class_j", "class_m"}`: the room, and whether a Class J (at least `CLASS_J_MIN_EARTH`) and a Class M (at least the smallest `CLASS_M_MASS_EARTH`) world fit."""
    room = room_earth(host_mass_earth, secondary_mass_earth)
    return {"room_earth": room, "class_j": room >= CLASS_J_MIN_EARTH, "class_m": room >= CLASS_M_MASS_EARTH[0]}


def fits(host_mass_earth, secondary_mass_earth, trojan_mass_earth, eccentricity=0.0, inclination_deg=0.0):
    """A whole body of this mass may sit at L4 or L5 of the pair: `trojans.allows` in Earth masses."""
    return trojans.allows(host_mass_earth, secondary_mass_earth, trojan_mass_earth, eccentricity, inclination_deg)


def can_form(trojan_mass_earth, solids_in_zone_earth, disc_lifetime_myr):
    """The formation gate (TB3, TB4): enough solids in the 1:1 zone and a gas disc that lasted long enough."""
    if trojan_mass_earth <= 0.0:
        return False
    return solids_in_zone_earth * FORMATION_EFFICIENCY >= trojan_mass_earth and disc_lifetime_myr >= MIN_DISC_LIFETIME_MYR
