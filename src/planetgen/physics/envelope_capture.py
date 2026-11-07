# planetgen/physics/envelope_capture.py

"""
The hydrogen-helium envelope a rocky planet captured and kept: sub-Neptunes and the radius valley (GEN.1323)
============================================================================================================

A core of a few Earth masses that finishes growing inside the gas disk captures a hydrogen-helium envelope. Close to its star the star's
X-rays and heat strip it again (`envelope_escape`, GEN.822); farther out the planet keeps it, and is a sub-Neptune of 2 to 4 Earth radii.
That is the radius valley (Fulton et al. 2017: about half the planets of 1 to 4 Earth radii inside 100 days are above 2.0). The envelope
and loss recipes of rows e1 and e2 (GEN.509, GEN.722, GEN.822) were tuned to the valley, but nothing drew an envelope for a rocky planet, so
no generated planet of 1 to 4 Earth radii came out above 2.0 (GEN.1323, found by the validation sweep of GEN.753).

* **The start.** A planet of `MIN_MASS_EARTH` to `MAX_MASS_EARTH` Earth masses draws the envelope it had when the disk went: a log-normal
  fraction of its mass, median `MEDIAN_START` (the middle of the recipe bracket 0.05 to 0.2), width `START_SPREAD` (the survey's), held at
  `MAX_START`. Below `REFERENCE_MASS_EARTH` (the peak of the core masses the valley was tuned on, [R3797]) the median
  falls as the mass over it to the power `MASS_EXPONENT`: a flat start would give every Earth-mass planet an envelope, and the Earth and Venus
  have none. The exponent is a stand-in with no source (registry GEN.725, section 6.196).
* **The loss.** `envelope_escape.envelope_left_earth` over the planet's real mass, light, host and age, with the mean light of its orbit.
  The planet keeps its drawn mass as its total (the masses the survey models are measured ones); the radius is the law of [R2741]
  (`bulk_structure.envelope_radius_earth`, at the 20 mbar level the transit radii and the valley survey use) of the planet's mass and what is left.
* **When.** Twice: a planet is made with the envelope the star's age and its orbit give it then (`planets.apply_envelope`, so a request judged by class
  sees it), and once the orbit and the star's age are final the loss is taken again (`planets.settle_envelope`, in the settle pass of the system, like a
  giant's settled radius, `planets.settle_giant`). A planet with no star keeps the whole envelope it started with (no light, no loss).

Use:

    envelope_capture.draw_start(mass_earth)                    # the starting fraction, from the stream bound now
    envelope_capture.kept(body, start)                         # (envelope fraction, radius_earth) after the loss, or None for a bare core
"""

import math

from planetgen.physics import bulk_structure, constants, envelope_escape, orbit_elements
from planetgen.physics import unit_registry as units
from planetgen.util import draw

MIN_MASS_EARTH = 1.0
"""float: The lightest planet that captures an envelope (the e1 and e2 rows start at 1 Earth mass; a lighter one cannot hold a nebular gas)."""

MAX_MASS_EARTH = 20.0
"""float: The heaviest rocky planet (`planet_kind` turns a heavier core into a giant, whose envelope its interior holds)."""

MEDIAN_START = 0.1
"""float: STAND-IN EC1 (T9, L3, registry 6.196): the median starting envelope, a fraction of the mass, of a planet of `REFERENCE_MASS_EARTH` or more;
the middle (geometric) of the recipe bracket `f_env0` 0.05 to 0.2 that both valley targets sit inside (GEN.822)."""

MAX_START = bulk_structure.CLOSED_FORM_RANGES["envelope"]["envelope"][1]
"""float: STAND-IN EC5 (T6, L1, registry 6.196): the largest starting envelope a rocky planet draws, 0.2 of its mass: the top of the range the law of Lopez and Fortney was
fitted over ([R2741]). The valley survey holds its start at 0.5 (`envelope_escape.MAX_ENVELOPE`); beyond 0.2 a planet of a few Earth masses would be a puffed-up
Saturn of 5 to 7 Earth radii, which no core of that mass captures."""

START_SPREAD = 0.6
"""float: The log-normal width of the starting envelope, in ln (the survey's, `radius_valley.ENVELOPE_SIGMA`; a stand-in)."""

REFERENCE_MASS_EARTH = 3.72
"""float: The core mass the valley survey is centred on ([R3797], 3.72 +0.45 -0.33 Earth masses, `radius_valley.CORE_MASS_MEDIAN_EARTH`)."""

MASS_EXPONENT = 3.0
"""float: STAND-IN EC2 (T7, L3, registry 6.196): below `REFERENCE_MASS_EARTH` the median start goes as (mass / reference) to this power. No source read;
the envelope a core captures grows steeply with its mass (a deeper well holds more of the disk's gas), and the Earth-mass planets of the sample have
none: at 3 the median start of a 1 Earth-mass planet is 0.2 percent of its mass (a bare rock) and of a 2 Earth-mass planet 1.6 percent."""

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_msun = units.converter("kg", "Msun")


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


def kept(body, start):
    """
    The envelope `body` (a rocky planet with its star, mass, orbit and the star's final age) still has after the loss of a starting
    envelope of `start` of its mass, and the radius it gives the planet.

    Returns:
        tuple or None: `(fraction of the mass, radius in Earth radii)`, or `None` for a planet that is a bare core (nothing kept).
    """
    if start <= 0.0:
        return None
    star = body.star
    mass_earth = _kg_to_mearth(body.mass)
    luminosity = star.luminosity / constants.SOLAR_LUMINOSITY
    distance_au = orbit_elements.flux_distance_au(body.distance, getattr(body, "eccentricity", 0.0))
    insolation = luminosity / distance_au ** 2
    left = envelope_escape.envelope_left_earth(mass_earth, start, insolation, _kg_to_msun(star.mass), luminosity, star.age)
    if left <= 0.0:
        return None
    fraction = left / (mass_earth + left)
    flux = max(insolation, bulk_structure.CLOSED_FORM_RANGES["envelope"]["flux_earth"][0])     # the law is fitted from 0.1 Earth's light up; a cold or starless planet is read there
    return fraction, bulk_structure.envelope_radius_earth(mass_earth, fraction, flux, star.age)
