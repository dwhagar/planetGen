# planetgen/physics/potassium_depletion.py

"""
Potassium depletion of a rocky body against CI chondrites (GEN.492, row o3 of the atmosphere recipes; docs/design/composition/radiation-model.md R15)
===============================================================================================================================

K/U of a body over K/U of CI chondrites measures how much of the moderately volatile potassium the body lost; uranium is refractory, so the
ratio does not move with magmatic differentiation [R1762]. `radiogenic` drew the depletion log-uniform from 0.03 to 1,
a placeholder. [R1003] (Planetary Science Journal 3:29, Table 1, READ) re-derive K/U for six differentiated bodies, with CI
chondrites from [R0270] (K 550 ppm, U 0.0077 ppm, K/U 71,429):

    bulk silicate Earth 0.136, Mars 0.243, Mercury surface 0.200 (0.055 to 1 Earth masses),
    Moon 0.034, Vesta 0.058, the angrite parent body 0.0018 (0.012 Earth masses down to a planetesimal).

Two classes follow from the gap between Moon and Mercury, and the module draws the log10 depletion from a normal distribution of each class's
mean and predictive spread (`s * sqrt(1 + 1/n)` of the three bodies): planets of 0.03 Earth masses or more, centred on 0.19 with 0.15 dex;
smaller bodies, centred on 0.015 with 0.94 dex. The draw never exceeds 1 (the chondrites of Table 1 run from 0.26 for CV to 1.2 for the
enstatite groups; a body that kept all its potassium has 1).

Limits, kept as registry rows (GEN.725, section 6.61): three bodies per class; the heliocentric trend the earlier recalled rule used
(depletion rising as `a^1.2`) is not in these six bodies (Mercury is more depleted than a trend from Earth to Mars would give), so none is
drawn; polluted white dwarfs say a majority of pollutants (38 of the 50 with sodium measured) are not depleted in sodium, a moderately volatile
element of nearly the same volatility as potassium ([R1832], READ), which the six differentiated bodies do not
show, so the undepleted channel is a falsifier and a later item, not a draw. Nothing in generation calls this yet: the wiring step is
`radiogenic.draw_content` taking the depletion from here with the body's mass, and CI K/U from `ci_k_over_u` (it carried 67,407).
"""

import json
import math
import os
import statistics

from planetgen.physics import elements
from planetgen.physics.elements import DataError
from planetgen.util import draw
from planetgen.util.checks import finite_domain

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "potassium_depletion.json")
"""str: The data file."""

MAX_DEPLETION = 1.0
"""float: A body keeps at most all of the potassium CI chondrites have per uranium."""


def load(path=DATA_FILE):
    """Reads and checks the table at `path`: `(reference K/U, {name: (mass_earth or None, K/U)}, class mass, source)`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the potassium table {path}: {error}") from error

    def positive(value, what):
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise DataError(f"{what} must be a positive number, not {value!r}")
        return float(value)

    reference = positive((document.get("reference") or {}).get("k_over_u"), "the CI K/U")
    bodies = {}
    for raw in document.get("bodies") or []:
        mass = raw.get("mass_earth")
        bodies[raw["name"]] = (None if mass is None else positive(mass, f"the mass of {raw['name']}"), positive(raw.get("k_over_u"), f"K/U of {raw['name']}"))
    if len(bodies) < 4:
        raise DataError("the potassium table has fewer than four bodies")
    return reference, bodies, positive(document.get("class_mass_earth"), "the class mass"), document.get("source", "")


CI_K_OVER_U, BODIES, CLASS_MASS_EARTH, SOURCE = load()
"""CI K/U by weight, the bodies `{name: (mass in Earth masses or None, K/U)}`, the Earth masses that separate the classes, the source."""


def ci_k_over_u():
    """The mass ratio of potassium to uranium in CI chondrites ([R0270]: K 550 ppm, U 0.0077 ppm)."""
    return CI_K_OVER_U


def depletion_of(name):
    """The depletion (K/U over CI K/U) of one of the six bodies in the table."""
    try:
        return BODIES[name][1] / CI_K_OVER_U
    except KeyError:
        raise DataError(f"{name!r} is none of {', '.join(BODIES)}") from None


def is_planet(mass_earth):
    """True for a body of `CLASS_MASS_EARTH` Earth masses or more."""
    return mass_earth >= CLASS_MASS_EARTH


def _members(planet):
    return [math.log10(k_over_u / CI_K_OVER_U) for mass, k_over_u in BODIES.values() if (mass is not None and is_planet(mass)) == planet]


def statistics_of(planet):
    """`(n, mean log10 depletion, predictive spread in dex)` of the planets (`planet` true) or of the smaller bodies."""
    logs = _members(planet)
    spread = statistics.stdev(logs) * math.sqrt(1.0 + 1.0 / len(logs))
    return len(logs), statistics.fmean(logs), spread


@finite_domain()
def median_depletion(mass_earth):
    """The depletion at the centre of the class of a body of `mass_earth` Earth masses."""
    return 10.0 ** statistics_of(is_planet(mass_earth))[1]


@finite_domain()
def draw_depletion(mass_earth, rng=draw):
    """
    The potassium depletion of a body of `mass_earth` Earth masses: log10 normal with its class's mean and predictive spread, never above
    `MAX_DEPLETION`. One `gauss` draw, whatever the class.
    """
    _, mean, spread = statistics_of(is_planet(mass_earth))
    return min(MAX_DEPLETION, 10.0 ** rng.gauss(mean, spread))
