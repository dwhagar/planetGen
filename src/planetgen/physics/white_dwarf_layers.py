# planetgen/physics/white_dwarf_layers.py

"""
The Layers of a Helium-Core or Oxygen-Neon White Dwarf (GEN.743)
================================================================

The downloaded cooling sequences (`datasets.white_dwarf_models`, [R0308]) are carbon-oxygen from 0.20 to 1.30 solar masses, so a white
dwarf the core-type rule (`stellar_abundances.white_dwarf_core`) calls helium-core (under 0.45) or oxygen-neon (over 1.05) has no layers
in them. This module gives those two their layers as the total mass fractions of the star, in the keys the carbon-oxygen layers
(`white_dwarf_models.layer_mass_fractions`) use, plus `neon` and `trace`:

* An **oxygen-neon** core is about 55 percent 16O and 30 percent 20Ne by mass, with 22Ne, 23Na and 24Mg as the rest [R0627: Sect. 2.2];
  the helium layer is the 3.24e-4, 1.82e-4, 0.78e-4 and 0.21e-4 solar masses of the 1.10, 1.16, 1.23 and 1.29 solar mass models, and the
  hydrogen of an atmosphere of hydrogen is the 1e-6 solar masses those models adopt as the most they expect. The log of the helium
  layer runs straight between the models and is held outside them; the layers are then fractions of the star's own mass.
* A **helium-core** star has a core of helium alone, whatever the layer on it: a hydrogen envelope that is a large share of a very
  light star's mass (the stand-in WL5) and nothing else, so the helium fraction is what the hydrogen leaves.

A helium-atmosphere star of either kind carries the thin hydrogen layer of the sequences' own thin case.
The radius of both is still the nearest carbon-oxygen sequence scaled by the zero-temperature relation (`stellar_composition`).
"""

import math

from planetgen import tuning


def _log_interpolated(points, x):
    """`y` at `x` along the line through `(x, log y)` of `points` (ascending in x), held at the first and last."""
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x <= x1:
            weight = (x - x0) / (x1 - x0)
            return math.exp(math.log(y0) * (1.0 - weight) + math.log(y1) * weight)
    return points[-1][1]    # pragma: no cover (the loop returns first)


def helium_core(mass_sol, envelope="thick"):
    """The layers of a helium-core white dwarf of `mass_sol`: hydrogen (thick: the envelope fraction of `WD_HELIUM_CORE_HYDROGEN_FRACTIONS`
    at the mass; thin: `WD_THIN_HYDROGEN_FRACTION`) and helium for the rest."""
    hydrogen = (_log_interpolated(tuning.WD_HELIUM_CORE_HYDROGEN_FRACTIONS, mass_sol)
                if envelope == "thick" else tuning.WD_THIN_HYDROGEN_FRACTION)
    return {"hydrogen": hydrogen, "helium": 1.0 - hydrogen, "carbon": 0.0, "oxygen": 0.0, "neon": 0.0, "trace": 0.0}


def oxygen_neon(mass_sol, envelope="thick"):
    """The layers of an oxygen-neon white dwarf of `mass_sol`: the hydrogen and helium layer masses over the star's mass and the core split
    of `WD_ONE_CORE_MASS_FRACTIONS` for the rest."""
    hydrogen = (tuning.WD_ONE_HYDROGEN_MASS_SOL / mass_sol if envelope == "thick" else tuning.WD_THIN_HYDROGEN_FRACTION)
    helium = _log_interpolated(tuning.WD_ONE_HELIUM_MASS_BY_WD_MASS_SOL, mass_sol) / mass_sol
    core = 1.0 - hydrogen - helium
    split = tuning.WD_ONE_CORE_MASS_FRACTIONS
    return {"hydrogen": hydrogen, "helium": helium, "carbon": 0.0, "oxygen": core * split["oxygen"], "neon": core * split["neon"],
            "trace": core * split["trace"]}
