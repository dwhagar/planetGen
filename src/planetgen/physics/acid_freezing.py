# planetgen/physics/acid_freezing.py

"""
The freezing curve of the sulfuric acid hydrates (GEN.1644, registry rows 6.287:AF1 and AF2; replaces the no-freezing rule SV7 of 6.228)
==========================================================================================================================================

An aqueous sulfuric acid droplet does not stay liquid as it cools: at 0.1 MPa the first solid appears at the freezing temperature of
its acid mole fraction, and which solid it is (ice, one of the five hydrates, or the pure acid) follows the same curve. [R4877] fitted the Gibbs energy of the solution and of each solid to the freezing points of Gable, Betz and Maron and of Kunzler and
Giauque, and printed the resulting phase diagram (Figure 1). `physics/data/acid_freezing.json` holds that diagram's stable liquidus, read off
the printed curve to about 1 K and anchored at the exact temperatures of Table 6.

* `freezing_point_k(x)`: the temperature below which a solution of acid mole fraction `x` has a solid in it.
* `solid_phase(x)`: which solid freezes out first there (`"ice"`, `"6.5H"`, `"4H"` (the tetrahydrate), `"3H"`, `"2H"`, `"1H"` (the monohydrate) or `"acid"`).
* `is_frozen(temperature_k, x)`: whether a droplet of that strength at that temperature has begun to freeze.
* `freezing_point_of_weight_percent(w)`: the same from a strength in weight percent, as the cloud step holds it (`sulfuric_vapor.mole_fraction`).

Limits. The curve is the equilibrium liquidus: a droplet can supercool well below it before it freezes (Gable, Betz and Maron measured metastable
eutectics 20 K or more lower), so `is_frozen` says where a droplet *can* be solid, not where it *must* be. It holds at 0.1 MPa; the melting curves
of these solids move by a few hundredths of a kelvin per bar, small against the 1 K of the reading, so the cloud step uses it at any pressure of a cloud.

Pure and unwired: GEN.1644 slice 1 makes the cold acid cloud read it in place of the 185 K floor (row 6.228:SV7).
"""

import json
import os

import numpy as np

from planetgen.physics import elements, sulfuric_vapor

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "acid_freezing.json")


class AcidFreezingError(ValueError):
    """A composition outside 0 to 1, or a rule file that is refused."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    xs = [row[0] for row in document["points"]]
    if xs != sorted(xs) or xs[0] != 0.0 or xs[-1] != 1.0 or len(set(xs)) != len(xs):
        raise AcidFreezingError("the points must run from mole fraction 0 to 1 in strictly rising order")
    return document


DOCUMENT = _load()
VERSION = DOCUMENT["version"]
ANCHORS = DOCUMENT["anchors"]
"""dict: The solid -> its Table 6 reference temperature, K (the melting point of ice, the pure acid, the monohydrate and the tetrahydrate; the
dihydrate, trihydrate and 6.5 hydrate are incongruent, so those are reference temperatures of the fit, not points of the liquidus)."""
EUTECTICS = {name: tuple(point) for name, point in DOCUMENT["eutectics"].items()}
"""dict: Name -> (acid mole fraction, temperature K) of each stable eutectic."""
_X = np.array([row[0] for row in DOCUMENT["points"]])
_T = np.array([row[1] for row in DOCUMENT["points"]])
_PHASE = [row[2] for row in DOCUMENT["points"]]
LOWEST_K = float(_T.min())
"""float: The lowest freezing temperature of the diagram, the ice and 6.5 hydrate eutectic, K. A solution below it is solid at every strength."""
HIGHEST_K = float(_T.max())
"""float: The highest freezing temperature, the melting point of the pure acid, K."""


def _check(x):
    if not 0.0 <= x <= 1.0:
        raise AcidFreezingError(f"an acid mole fraction is 0 to 1, got {x!r}")


def freezing_point_k(acid_mole_fraction):
    """
    The temperature, K, below which a solution of `acid_mole_fraction` (acid moles over acid plus water moles) has a solid in it at equilibrium.

    Raises:
        AcidFreezingError: The fraction is outside 0 to 1.
    """
    _check(acid_mole_fraction)
    return float(np.interp(acid_mole_fraction, _X, _T))


def solid_phase(acid_mole_fraction):
    """The solid that freezes out first from a solution of `acid_mole_fraction`: `"ice"`, `"6.5H"`, `"4H"`, `"3H"`, `"2H"`, `"1H"` or `"acid"`."""
    _check(acid_mole_fraction)
    index = int(np.searchsorted(_X, acid_mole_fraction, side="left"))
    return _PHASE[min(index, len(_PHASE) - 1)]


def is_frozen(temperature_k, acid_mole_fraction):
    """Whether a droplet of that strength at `temperature_k` is at or below its freezing temperature (a solid can be in it)."""
    return temperature_k <= freezing_point_k(acid_mole_fraction)


def freezing_point_of_weight_percent(weight_percent_acid):
    """The freezing temperature, K, of a solution of `weight_percent_acid` percent sulfuric acid (0 to 100)."""
    if not 0.0 <= weight_percent_acid <= 100.0:
        raise AcidFreezingError(f"a weight percent is 0 to 100, got {weight_percent_acid!r}")
    return freezing_point_k(sulfuric_vapor.mole_fraction(weight_percent_acid))
