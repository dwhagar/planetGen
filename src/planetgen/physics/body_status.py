# planetgen/physics/body_status.py

"""
The status of a body: planet, dwarf planet, minor body or moon, and the numbers behind the word (GEN.1063; docs/design/systems/belt-evolution.md section 3.2 and decision 2, belt-sources.md)
====================================================================================================================================================================================

Boss's rule is that the status word is a pure function of the body's mass, its composition, its orbit and the neighbours in the ledger. Nothing here reads the Solar System as an input: Earth, Mars, Ceres, Pluto, Eris, Jupiter
and the Trojans are test targets. **Both the IAU word and the numbers are kept** (section 11, decision 2), so either can be used:

1. **Roundness** (`round_mass_kg`, `roundness`). A solid body turns round when its own gravity beats its strength: the potato radius is 200 km for ice and 300 km for rock (3.3e19 and 3.4e20 kg, BR8) and a fully relaxed
   body is 1e20 to 1e21 kg of ice and 1e21 to 1e22 kg of rock (BR9). The shape follows the compaction taper (`compaction.taper_weight`, 1e19 to 1e21 kg for ice, linear in log mass): a roundness of 0 at its
   start and 1 at its end, 0.5 at 1e20 kg, which is the lower edge of the relaxed ice band. Rock turns round one decade later (the ratio of the two BR9 bands), so the window is shifted by a factor of 10. Between ice and rock the
   shift is linear in the ice mass fraction. A body is **round** at a roundness of 0.5 or more.
2. **Soter's mu** (`soter_mu`): the mass of the body over the mass of the other bodies in its zone. A body above `SOTER_MU_PLANET` has cleared its zone (Earth 1.7e6, Ceres 1/3, Pluto 0.08, BR6). The zone is the orbit's
   band of `ZONE_HILL_RADII` Hill radii on each side (Margot's C, BR22), and `zone_mass_kg` adds up what the ledger puts there: stored bodies by their orbits, and the leftover mass of a belt, spread evenly across its width.
3. **Margot's Pi** (`margot_pi`): the body's mass over the mass that clears an orbit in the main-sequence lifetime, `1.9e-4 C^(3/2) (a / 1 AU)^(9/8) (M* / Msun)^(5/2)` Earth masses with `C = 2 sqrt 3` (BR7; Jupiter 4e4, Mars 54).
   Margot's `807` in the design note is `1 / (1.9e-4 C^(3/2))` = 816 with the printed digits. It is the cross-check on Soter's word: it needs no ledger.
4. **Stern and Levison's Lambda** (`stern_levison_lambda`): `k m^2 / P`, with `m` in Earth masses per solar mass of the star, over the orbital period in years. Their `k` was not read (second hand, Soter 2006); the one here is fixed so that Lambda is 1 for a body of exactly Margot's
   clearing mass at 1 AU round the Sun (BR23), which puts Earth at 7e5 and Ceres at 4e-3. Above 1 the body clears its orbit.
5. **The word** (`status`): a body that orbits another planet or dwarf planet is a `moon`. One that orbits the star is a `planet` if it is round and has cleared its zone (mu above `SOTER_MU_PLANET`), a `dwarf planet` if it is round and
   has not, and a `minor body` if it is not round, whatever its zone. Soter's mu decides, Margot's Pi is stored as the second opinion (`BodyStatus.pi_agrees`).

Limits. A zone with no other mass has an infinite mu. A body that orbits a planet has no zone of its own, so a moon carries only its roundness. Margot's Pi and the Lambda take the orbit as circular. Nothing calls this yet: the ledger
and the system wiring (GEN.1069 and later) ride merge trains.

Use:

    others = zone_mass_kg(a_au, mass_kg, star_mass_solar, neighbours, belts)           # neighbours: (a_au, mass_kg); belts: (lower_au, upper_au, mass_kg)
    result = status(mass_kg, ice_fraction, a_au, star_mass_solar, others)               # a BodyStatus, word and numbers
"""

import math
from dataclasses import dataclass

from planetgen.physics import compaction, constants, unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")

PLANET = "planet"
DWARF_PLANET = "dwarf planet"
MINOR_BODY = "minor body"
MOON = "moon"
WORDS = (PLANET, DWARF_PLANET, MINOR_BODY, MOON)

SOTER_MU_PLANET = 100.0
"""float: A body whose mass is more than this many times the mass of the others in its zone is a planet (Soter 2006; the Solar System's gap is mu 0.3 to 1e5). [Stand-in BR6, tier T1, tweak L0.]"""

MARGOT_COEFFICIENT = 1.9e-4
"""float: Margot's clearing mass at 1 AU round the Sun is this times `C^(3/2)` Earth masses (Margot 2015, printed). [Stand-in BR7, tier T1, tweak L0.]"""

MARGOT_C = 2.0 * math.sqrt(3.0)
"""float: Margot's `C`, the width in Hill radii of the zone a planet must clear (2 sqrt 3, printed). [Stand-in BR7, tier T1, tweak L0.]"""

ZONE_HILL_RADII = MARGOT_C
"""float: The half-width of the zone of a body, in its Hill radii: Margot's `C`, so Soter's zone and Margot's agree. [Stand-in BR22, tier T2, tweak L0.]"""

LAMBDA_K = 1.0 / (MARGOT_COEFFICIENT * MARGOT_C ** 1.5) ** 2
"""float: The `k` of Stern and Levison's `Lambda = k m^2 / P` (m in Earth masses per solar mass of the star, P in years), fixed so that Lambda is 1 at Margot's clearing mass at 1 AU round the Sun. [Stand-in BR23, tier T2, tweak L0.]"""

ROUND_ROCK_SHIFT_DECADES = 1.0
"""float: Rock turns round this many decades of mass above ice (the BR9 bands are 1e20 to 1e21 kg for ice and 1e21 to 1e22 kg for rock). [Stand-in BR24, tier T4, tweak L1.]"""

ROUND_THRESHOLD = 0.5
"""float: A body is round at a roundness of at least this (the log midpoint of the taper, 1e20 kg for ice). [Stand-in BR24, tier T4, tweak L1.]"""

_MARGOT_MASS_EARTH_AT_1AU = MARGOT_COEFFICIENT * MARGOT_C ** 1.5


@dataclass(frozen=True)
class BodyStatus:
    """The word and the numbers: `word` is one of `WORDS`; `round_mass_kg` and `roundness` for the body's composition; `mu`, `pi` and `lam` are None for a moon
    (and `lam` None where no period is given); `pi_agrees` says whether Margot's Pi puts the body on the same side as Soter's mu (None for a moon)."""

    word: str
    round_mass_kg: float
    roundness: float
    mu: float | None
    pi: float | None
    lam: float | None
    pi_agrees: bool | None


def round_mass_kg(ice_fraction):
    """The mass at which a body of this ice mass fraction (0 rock, 1 ice) is half way to round: 1e20 kg of ice, 1e21 kg of rock, linear in log mass between."""
    fraction = min(max(ice_fraction, 0.0), 1.0)
    midpoint = 0.5 * (math.log10(compaction.TAPER_START_KG) + math.log10(compaction.TAPER_END_KG))
    return 10.0 ** (midpoint + ROUND_ROCK_SHIFT_DECADES * (1.0 - fraction))


def roundness(mass_kg, ice_fraction):
    """0 for a rubble pile or a potato, 1 for a relaxed body: 1 minus the compaction taper's weight, shifted to the body's composition."""
    if mass_kg <= 0.0:
        raise ValueError("the mass must be positive")
    shift = round_mass_kg(ice_fraction) / 10.0 ** (0.5 * (math.log10(compaction.TAPER_START_KG) + math.log10(compaction.TAPER_END_KG)))
    return 1.0 - compaction.taper_weight(mass_kg / shift)


def hill_radius_au(a_au, mass_kg, star_mass_solar):
    """The Hill radius of a body of this mass at `a_au` round a star of `star_mass_solar`, in AU."""
    return a_au * (mass_kg / (3.0 * star_mass_solar * constants.SOLAR_MASS_TO_KG)) ** (1.0 / 3.0)


def zone_mass_kg(a_au, mass_kg, star_mass_solar, neighbours=(), belts=()):
    """The mass of the other bodies in the zone of a body: `neighbours` are `(a_au, mass_kg)` pairs, each counted whole if its orbit lies in the zone (the body itself, at `a_au` and `mass_kg`, is left out);
    `belts` are `(lower_au, upper_au, mass_kg)`, each counted by the share of its width inside the zone."""
    half = ZONE_HILL_RADII * hill_radius_au(a_au, mass_kg, star_mass_solar)
    low, high = a_au - half, a_au + half
    total = 0.0
    for other_a, other_mass in neighbours:
        if other_a == a_au and other_mass == mass_kg:
            continue
        if low <= other_a <= high:
            total += other_mass
    for lower, upper, belt_mass in belts:
        overlap = min(upper, high) - max(lower, low)
        if overlap > 0.0 and upper > lower:
            total += belt_mass * overlap / (upper - lower)
    return total


def soter_mu(mass_kg, others_kg):
    """Soter's mu: the body's mass over the mass of the others in its zone; infinite when there are none."""
    if mass_kg <= 0.0:
        raise ValueError("the mass must be positive")
    return math.inf if others_kg <= 0.0 else mass_kg / others_kg


def margot_clearing_mass_earth(a_au, star_mass_solar):
    """Margot's mass that clears an orbit at `a_au` round a star of `star_mass_solar` in its main-sequence lifetime, Earth masses."""
    return _MARGOT_MASS_EARTH_AT_1AU * a_au ** (9.0 / 8.0) * star_mass_solar ** 2.5


def margot_pi(mass_kg, a_au, star_mass_solar):
    """Margot's Pi: the body's mass over `margot_clearing_mass_earth`; at least 1 is a planet."""
    return _kg_to_mearth(mass_kg) / margot_clearing_mass_earth(a_au, star_mass_solar)


def stern_levison_lambda(mass_kg, a_au, star_mass_solar):
    """Stern and Levison's Lambda, `k (m / M*)^2 / P` with the mass in Earth masses, the star's in solar masses and the period in years; above 1 the body clears its orbit."""
    period_years = a_au ** 1.5 / math.sqrt(star_mass_solar)
    return LAMBDA_K * (_kg_to_mearth(mass_kg) / star_mass_solar) ** 2 / period_years


def status(mass_kg, ice_fraction, a_au, star_mass_solar, others_kg, parent_is_planet=False):
    """The `BodyStatus` of a body of `mass_kg` and `ice_fraction` at `a_au` round a star of `star_mass_solar`, with `others_kg` the mass of the others in its zone (`zone_mass_kg`).
    A body that orbits a planet or dwarf planet (`parent_is_planet`) is a moon and carries only its roundness."""
    mid = round_mass_kg(ice_fraction)
    shape = roundness(mass_kg, ice_fraction)
    if parent_is_planet:
        return BodyStatus(MOON, mid, shape, None, None, None, None)
    mu = soter_mu(mass_kg, others_kg)
    pi = margot_pi(mass_kg, a_au, star_mass_solar)
    lam = stern_levison_lambda(mass_kg, a_au, star_mass_solar)
    cleared = mu > SOTER_MU_PLANET
    if shape < ROUND_THRESHOLD:
        word = MINOR_BODY
    else:
        word = PLANET if cleared else DWARF_PLANET
    return BodyStatus(word, mid, shape, mu, pi, lam, cleared == (pi >= 1.0))
