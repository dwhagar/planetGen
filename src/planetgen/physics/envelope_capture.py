# planetgen/physics/envelope_capture.py

"""
The hydrogen-helium envelope a rocky planet captured and kept: sub-Neptunes and the radius valley (GEN.1323)
============================================================================================================

A core of a few Earth masses that finishes growing inside the gas disk captures a hydrogen-helium envelope. Close to its star the star's
X-rays and heat strip it again (`envelope_escape`, GEN.822); farther out the planet keeps it, and is a sub-Neptune of 2 to 4 Earth radii.
That is the radius valley ([R1478]: about half the planets of 1 to 4 Earth radii inside 100 days are above 2.0). The envelope
and loss recipes of rows e1 and e2 (GEN.509, GEN.722, GEN.822) were tuned to the valley, but nothing drew an envelope for a rocky planet, so
no generated planet of 1 to 4 Earth radii came out above 2.0 (GEN.1323, found by the validation sweep of GEN.753).

* **The start.** A planet of `MIN_MASS_EARTH` to `MAX_MASS_EARTH` Earth masses draws the envelope it had when the disk went: a log-normal
  fraction of its mass, median `MEDIAN_START` (the middle of the recipe bracket 0.05 to 0.2), width `START_SPREAD` (the survey's), held at
  `MAX_START`. Below `REFERENCE_MASS_EARTH` (the peak of the core masses the valley was tuned on, [R3797]) the median
  falls as the mass over it to the power `MASS_EXPONENT`: a flat start would give every Earth-mass planet an envelope, and the Earth and Venus
  have none. The exponent is a stand-in with no source (registry GEN.725, section 6.206).
* **The loss.** `envelope_escape.envelope_left_earth` over the planet's real mass, light, host and age, with the mean light of its orbit.
  The planet keeps its drawn mass as its total (the masses the survey models are measured ones); the radius is the law of [R2741]
  (`bulk_structure.envelope_radius_earth`, at the 20 mbar level the transit radii and the valley survey use) of the planet's mass and what is left.
* **When.** Twice: a planet is made with the envelope the star's age and its orbit give it then (`planets.apply_envelope`, so a request judged by class
  sees it), and once the orbit and the star's age are final the loss is taken again (`planets.settle_envelope`, in the settle pass of the system, like a
  giant's settled radius, `planets.settle_giant`). A planet with no star keeps the whole envelope it started with (no light, no loss).

Use:

    envelope_capture.draw_start(mass_earth)                    # the starting fraction, from the stream bound now
    envelope_capture.kept(body, start)                         # (envelope fraction, radius_earth) after the loss, or None for a bare core
    envelope_capture.kept(body, start, tides=True)             # the same with the obliquity tide of the planet's own tilt (GEN.1347)
"""

import math

from planetgen.physics import bulk_structure, constants, envelope_escape, helium_fractionation, orbit_elements, super_puff
from planetgen.physics import radiation_field
from planetgen.physics import unit_registry as units
from planetgen.util import draw

HELIUM_LOSS_GATE = 0.5
"""float: The share of its starting envelope a planet must have lost before the helium history is worked out: [R2839] find a planet must lose at least half of its envelope to reach
a helium mass fraction of 0.4 (registry GEN.725 row HF3 of section 6.77)."""

HELIUM_HISTORY_POINTS = 40
"""int: The most steps of the loss history `helium_fractionation` is given (numerical: the loss takes up to thousands of steps)."""

MIN_MASS_EARTH = 1.0
"""float: The lightest planet that captures an envelope (the e1 and e2 rows start at 1 Earth mass; a lighter one cannot hold a nebular gas)."""

MAX_MASS_EARTH = 20.0
"""float: The heaviest rocky planet (`planet_kind` turns a heavier core into a giant, whose envelope its interior holds)."""

MEDIAN_START = 0.1
"""float: STAND-IN EC1 (T9, L3, registry 6.206): the median starting envelope, a fraction of the mass, of a planet of `REFERENCE_MASS_EARTH` or more;
the middle (geometric) of the recipe bracket `f_env0` 0.05 to 0.2 that both valley targets sit inside (GEN.822)."""

MAX_START = bulk_structure.CLOSED_FORM_RANGES["envelope"]["envelope"][1]
"""float: STAND-IN EC5 (T6, L1, registry 6.206): the largest starting envelope a rocky planet draws, 0.2 of its mass: the top of the range the law of Lopez and Fortney was
fitted over [R2741]. The valley survey holds its start at 0.5 (`envelope_escape.MAX_ENVELOPE`); beyond 0.2 a planet of a few Earth masses would be a puffed-up
Saturn of 5 to 7 Earth radii, which no core of that mass captures."""

START_SPREAD = 0.6
"""float: The log-normal width of the starting envelope, in ln (the survey's, `radius_valley.ENVELOPE_SIGMA`; a stand-in)."""

REFERENCE_MASS_EARTH = 3.72
"""float: The core mass the valley survey is centered on ([R3797], 3.72 +0.45 -0.33 Earth masses, `radius_valley.CORE_MASS_MEDIAN_EARTH`)."""

MASS_EXPONENT = 3.0
"""float: STAND-IN EC2 (T7, L3, registry 6.206): below `REFERENCE_MASS_EARTH` the median start goes as (mass / reference) to this power. No source read;
the envelope a core captures grows steeply with its mass (a deeper well holds more of the disk's gas), and the Earth-mass planets of the sample have
none: at 3 the median start of a 1 Earth-mass planet is 0.2 percent of its mass (a bare rock) and of a 2 Earth-mass planet 1.6 percent."""

PUFF_MIN_MASS_EARTH = 2.0
"""float: STAND-IN EC8 (T4, L0, registry 6.206): the lightest planet that can start as a super-puff, the lower edge of Lee and Chiang's class (2 to 6 Earth masses, [R2609]);
the heaviest is `super_puff.START_MAX_MASS_EARTH` (row PF1)."""

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_msun = units.converter("kg", "Msun")
_DAYS_PER_YEAR = 365.25


def in_range(mass_earth):
    """Whether a planet of `mass_earth` Earth masses can capture an envelope."""
    return MIN_MASS_EARTH <= mass_earth <= MAX_MASS_EARTH


def median_start(mass_earth):
    """The median starting envelope, a fraction of the mass, of a planet of `mass_earth` Earth masses."""
    return MEDIAN_START * min(1.0, mass_earth / REFERENCE_MASS_EARTH) ** MASS_EXPONENT


def draw_start(mass_earth):
    """The starting envelope of a planet of `mass_earth` Earth masses: a fraction of its mass, one draw from the stream bound now (0 outside the mass range, no draw)."""
    if not in_range(mass_earth):
        return 0.0
    return min(median_start(mass_earth) * math.exp(START_SPREAD * draw.gauss()), MAX_START)


def period_days(body):
    """The orbital period of a bound planet, days (Kepler's third law from its distance and its star's mass)."""
    mass_solar = _kg_to_msun(body.star.mass)
    return _DAYS_PER_YEAR * math.sqrt(body.distance ** 3 / mass_solar)


def puff_eligible(mass_earth, period_d):
    """Whether a planet can start as a super-puff: a sub-Neptune of `PUFF_MIN_MASS_EARTH` to `super_puff.START_MAX_MASS_EARTH` Earth masses with a period inside
    `super_puff.START_MAX_PERIOD_D` (the Kepler sub-Neptunes the occurrence share counts)."""
    return PUFF_MIN_MASS_EARTH <= mass_earth <= super_puff.START_MAX_MASS_EARTH and period_d <= super_puff.START_MAX_PERIOD_D


def draw_puff_start(mass_earth, period_d):
    """
    The starting envelope of a planet that is a super-puff, or 0.0 when it is not one (GEN.1343, W15; [R2609]: a dust-free atmosphere accreted in the outer disk,
    then migration in). An eligible planet (`puff_eligible`) is one with the probability `super_puff.START_SHARE` (row PF3); it then starts with an
    envelope log-uniform over `super_puff.PUFF_ENVELOPE` (row PF2). Two draws from the stream bound now for an eligible planet, none for any other.
    """
    if not puff_eligible(mass_earth, period_d):
        return 0.0
    is_puff = draw.random() < super_puff.START_SHARE
    low, high = super_puff.PUFF_ENVELOPE
    fraction = math.exp(math.log(low) + draw.random() * (math.log(high) - math.log(low)))
    return fraction if is_puff else 0.0


def tidal_q_prime(body):
    """
    The tidal quality factor `Q'` of `body` (GEN.1347, row PF5): log-uniform over `super_puff.TIDAL_Q_PRIME_RANGE`, from a stream of the planet's own (its mass, formation
    distance and star), so the run's stream does not move and every pass draws the same value.
    """
    from planetgen.physics import planets
    low, high = super_puff.TIDAL_Q_PRIME_RANGE
    stream = draw.Stream(f"tidal-q-prime:{float(body.mass)!r}:{float(planets.formation_distance_au(body))!r}:{float(body.star.mass)!r}")
    with draw.bound(stream):
        return math.exp(math.log(low) + draw.random() * (math.log(high) - math.log(low)))


def _loss(body, start, tides, trace=None):
    """
    The loss of a starting envelope of `start` of the mass of `body` (a rocky planet with its star, mass, orbit and the star's final age):
    `(mass in Earth masses, insolation, flux (the insolation floored at the fit's), envelope left in Earth masses)`. With `tides` the radius the loss sees is the one
    the obliquity tide gave (`kept_with_tide`). `trace`, a list, is filled with the history of the envelope (`envelope_escape.envelope_left_earth`).
    """
    star = body.star
    mass_earth = _kg_to_mearth(body.mass)
    luminosity = star.luminosity / constants.SOLAR_LUMINOSITY
    distance_au = orbit_elements.flux_distance_au(body.distance, getattr(body, "eccentricity", 0.0))
    insolation = luminosity / distance_au ** 2 + radiation_field.companion_insolation_sol(star, distance_au)     # GEN.1441: the other stars' light too
    flux_floor = bulk_structure.CLOSED_FORM_RANGES["envelope"]["flux_earth"][0]       # the law is fitted from 0.1 Earth's light up; a cold or starless planet is read there
    flux = max(insolation, flux_floor)
    step_ratio = None
    if tides and getattr(body, "axial_tilt_deg", None) is not None:
        def step_ratio(mass, fraction, radius):
            return tidal_ratio(body, mass, fraction, flux, radius)
    left = envelope_escape.envelope_left_earth(mass_earth, start, insolation, _kg_to_msun(star.mass), luminosity, star.age, distance_au=distance_au if radiation_field.companions_of(star) else None,
                                               emitters=radiation_field.xuv_emitters(star, distance_au), radius_ratio=step_ratio, trace=trace)      # GEN.1442: the other stars' XUV too
    return mass_earth, insolation, flux, left


def helium_share(body, start, tides=False):
    """
    The helium share of the mass of the hydrogen-helium part of the envelope `body` still has after the loss of a starting envelope of `start` of its mass, when the loss was large enough to
    have changed it (GEN.1527, `helium_fractionation`, W16): hydrogen drags helium out only while the escape flux is large, so a planet that lost
    `HELIUM_LOSS_GATE` of its envelope or more is left richer in helium than the star's gas. `None` for a planet that lost less (its envelope is the star's gas, no history is
    worked out), that kept nothing, or that has no star's gas to start from. The start is the star's own share (`bulk_draw.host_ratios`); the history is the loss of
    `kept_with_tide` run once more with its steps kept.
    """
    if start <= 0.0 or getattr(body, "star", None) is None:
        return None
    from planetgen.physics import bulk_draw      # imported here: bulk_draw is loaded by the planet before this module is used
    mass_earth, insolation, flux, left = _loss(body, start, tides)
    start_earth = min(start, envelope_escape.MAX_ENVELOPE) * mass_earth / (1.0 - min(start, envelope_escape.MAX_ENVELOPE))
    if left <= 0.0 or left > (1.0 - HELIUM_LOSS_GATE) * start_earth:
        return None
    trace = []
    _loss(body, start, tides, trace=trace)
    stride = max(1, len(trace) // HELIUM_HISTORY_POINTS)
    history = trace[::stride]
    if history[-1] != trace[-1]:
        history.append(trace[-1])
    _, y = helium_fractionation.helium_fraction_from_history(mass_earth, history, flux, bulk_draw.host_ratios(body.star).y_envelope)
    return y


def kept(body, start, tides=False):
    """
    The envelope `body` (a rocky planet with its star, mass, orbit and the star's final age) still has after the loss of a starting
    envelope of `start` of its mass, and the radius it gives the planet. With `tides` the radius is multiplied by the obliquity-tide inflation
    (GEN.1347, `super_puff.tidal_inflation_ratio`) of the planet's own tilt (`axial_tilt_deg`) and tidal quality (`tidal_q_prime`), where its tilt is known and the tidal power
    over its starlight reaches `super_puff.TIDAL_MIN_LUMINOSITY_RATIO`; a planet with no tilt yet has none. The inflation is also in the loss: the XUV takes
    the envelope at the radius the planet has (GEN.1480, `envelope_escape.envelope_left_earth` with `radius_ratio`).

    Returns:
        tuple or None: `(fraction of the mass, radius in Earth radii)`, or `None` for a planet that is a bare core (nothing kept).
    """
    result = kept_with_tide(body, start, tides)
    return None if result is None else result[:2]


def kept_with_tide(body, start, tides=False):
    """
    `kept` with the factor the tide put on the radius as a third number: `(fraction of the mass, radius in Earth radii, tide ratio)` (the ratio is 1.0 for a
    planet it did not inflate), or `None` for a bare core. The caller that has to know whether a planet was inflated reads the ratio here and does not work
    out a second, uninflated envelope to compare with (that one would keep more envelope than the tide's loss leaves).
    """
    if start <= 0.0:
        return None
    star = body.star
    mass_earth, insolation, flux, left = _loss(body, start, tides)
    if left <= 0.0:
        return None
    fraction = left / (mass_earth + left)
    radius = super_puff.radius_earth(mass_earth, fraction, flux, star.age)
    ratio = tidal_ratio(body, mass_earth, fraction, flux, radius) if tides else 1.0
    return fraction, radius * ratio, ratio


def tidal_ratio(body, mass_earth, fraction, flux_earth, radius_earth):
    """The inflation of a planet of `mass_earth`, `fraction` envelope, `flux_earth` Earth fluxes and (uninflated) radius `radius_earth` by its obliquity tides: 1 when it has no tilt yet,
    its tilt is under `super_puff.MIN_OBLIQUITY_DEG`, its light is under `super_puff.TIDAL_MIN_FLUX_EARTH` or the tide is too weak against the starlight."""
    tilt = getattr(body, "axial_tilt_deg", None)
    if tilt is None or flux_earth < super_puff.TIDAL_MIN_FLUX_EARTH or not super_puff.covers(mass_earth, fraction):
        return 1.0
    q_term = super_puff.q_term(tidal_q_prime(body), tilt)
    if not math.isfinite(q_term):
        return 1.0
    radius_m = radius_earth * constants.EARTH_RADIUS_KM * 1.0e3
    distance_m = body.distance * constants.AU_TO_KM * 1.0e3        # the tide is worked out for a circular orbit at the semi-major axis ([R3079] eq. 11)
    if super_puff.tidal_luminosity_ratio(body.star.mass, radius_m, distance_m, flux_earth, q_term) < super_puff.TIDAL_MIN_LUMINOSITY_RATIO:
        return 1.0
    return super_puff.tidal_inflation_ratio(mass_earth, fraction, flux_earth, q_term)
