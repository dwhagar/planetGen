# planetgen/physics/comet_activity.py

"""
How near the light a comet must come to turn active (GEN.1491, docs/design/systems/comet-orbital-realism.md "Activity")
========================================================================================================================

A comet shows a coma when the light it gets sublimates its water ice fast enough. The old rule fixed that at a perihelion of 3 AU (`tuning.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU`): the
distance at which the Sun's light does it. Round another star the same light is at another distance: an airless body in the starlight of a star of luminosity `L` takes the same flux
(and so the same temperature, and the same sublimation) at `sqrt(L / Lsun)` times the Sun's distance. So a comet's perihelion is read as the distance from a Sun at which it would be
lit just as much (`solar_equivalent_au`), and the threshold, chances and floor stay the Sun's. The companions' light is already in the distance it is given (`radiation_field.lit_distance_au`,
GEN.1444). A host of the Sun's luminosity, alone, gives the numbers of before.

    solar_equivalent_au(lit_distance_au, luminosity_w)

The scaling is the flux of an airless body, not the disk's: the disk's temperature law (`comet_ice_budget.disk_temperature_k`, `L^0.25 r^-0.62`) is for the dust of a passively heated disk
and moves the ice line of water as `L^0.4`, which is not what a bare nucleus in sunlight feels. The Sun's 3 AU is the old stand-in, unchanged (row CA1); the equal-flux scaling is row CA2.
"""

import math

from planetgen.physics import constants


def _scale(luminosity_w):
    """The distance round a star of this luminosity that is lit like 1 AU round the Sun: `sqrt(L / Lsun)`; 1.0 for a star with no luminosity given."""
    if not luminosity_w or luminosity_w <= 0.0:
        return 1.0
    return math.sqrt(luminosity_w / constants.SOLAR_LUMINOSITY)


def solar_equivalent_au(lit_distance_au, luminosity_w):
    """The distance round the Sun that is lit like `lit_distance_au` round a star of `luminosity_w`: the distance over `sqrt(L / Lsun)`. The distance itself for the Sun's luminosity or none given."""
    return lit_distance_au / _scale(luminosity_w)

