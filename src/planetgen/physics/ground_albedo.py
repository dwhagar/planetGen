# planetgen/physics/ground_albedo.py

"""
The Bond albedo of a dry body's ground (GEN.684 part A; GEN.462)
================================================================

A body with no surface water takes its Bond albedo from its ground, not from a range drawn for every dry rock: the reflectance of the dominant dry ground in the three visible
bands of the materials table (`surface_composition.band_albedo`, the table's 26 rows from [R3425] and the USGS spectral library), times `VISIBLE_TO_BOND`. The ground type is the
body's own as far as it is known before its air is solved (`dry_ground_type`: carbon rock, a stripped metal body, else the basalt the melt of a mantle makes), and the grain
is a powder on a surface older than `surface_composition.RANGES["young_surface_gy"]`, else a slab.

Space weathering darkens a surface that has no air over it: `A(t) = A_fresh [1 - d (1 - exp(-t / tau_w))]` with `t` the age of the surface (the body's age), `tau_w` log-uniform
in 1 to 2 Gyr and `d` uniform in 0.3 to 0.7 (a parametrisation of the design note, placeholders to fit against lunar and asteroid reflectance; the exponential decay of the
unweathered fraction is [R4713]). A body whose air is solved keeps the fresh albedo, because a surface under an atmosphere is not darkened this way; a body that ends with none
takes `weathered_albedo` (`planets.calculate_atmospheric_conditions`). Basaltic powder gives a fresh 0.219, and 0.125 after 3 Gyr at the middle of the ranges: the Moon is 0.11 and
Mercury 0.088. A body with surface water keeps its own draw (`climate_draws.albedo_range`) and the ice law.

Use:

    albedo = fresh_bond_albedo(body)                     # before the air is solved
    albedo = weathered_albedo(albedo, surface_age_gy=4.5)     # when the body ends with no air
"""

import math

from planetgen.physics import atmosphere_colour, surface_composition
from planetgen.util import draw

TAU_W_GY = (1.0, 2.0)
"""tuple: The space-weathering time, log-uniform (the published range is 1 Myr to 2 Gyr; this is the design note's placeholder)."""
DEPTH = (0.3, 0.7)
"""tuple: The darkening amplitude `d`, uniform (a placeholder of the design note)."""

VISIBLE_TO_BOND = 1.0
"""float: The factor from the mean laboratory reflectance of the three visible bands to a Bond albedo (stand-in GA1 of the registry, tweak L1). A reflectance factor of a powder is neither a geometric
nor a Bond albedo and the conversion over the whole spectrum is an open gap; one makes the fresh basalt of the table (0.219) a Moon of 0.117 after 3.5 Gyr at the middle of the weathering ranges,
against the measured 0.11."""

CARBON = "carbon"
METAL = "metal"
BASALTIC = "basaltic"


def dry_ground_type(body):
    """The dominant dry ground of a body before its air is solved: `carbon` for a carbon rock, `metal` for a body whose mantle is stripped (`surface_composition.RANGES["metal_stripped"]`),
    else `basaltic`, which the table of dominant grounds and the melt of a mantle give for an ordinary rock (stand-in GA2: the melt's own type needs the heat flow, which comes after)."""
    if getattr(body, "carbon_rock", False):
        return CARBON
    if (getattr(body, "metal_mass_fraction", None) or 0.0) >= surface_composition.RANGES["metal_stripped"]:
        return METAL
    return BASALTIC


def surface_age_gy(body):
    """The age of the body's surface, Gyr: the body's age (the star's; a rogue body takes the reference age)."""
    from planetgen import tuning
    return getattr(getattr(body, "star", None), "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY


def fresh_bond_albedo(body):
    """
    The Bond albedo of the body's dry ground before any weathering: `VISIBLE_TO_BOND` times the mean of the table's three visible bands for the ground type `dry_ground_type` names (powder, or slab on
    a surface younger than `surface_composition.RANGES["young_surface_gy"]`); a ground the table has no spectrum for (the metal) takes `atmosphere_colour.STAND_IN_GROUND_ALBEDO`.
    """
    ground = dry_ground_type(body)
    grain = surface_composition.SLAB if surface_age_gy(body) < surface_composition.RANGES["young_surface_gy"] else surface_composition.POWDER
    composition = surface_composition.SurfaceComposition(ground_type=ground, mix={ground: 1.0}, components={}, grain=grain)
    bands, covered = surface_composition.band_albedo(composition)
    values = [bands[band] for band in surface_composition.BANDS]
    if not covered > 0.0 or any(value is None for value in values):
        return atmosphere_colour.STAND_IN_GROUND_ALBEDO
    return VISIBLE_TO_BOND * sum(values) / len(values)


def weathering_factor(surface_age_gy, tau_w_gy, depth):
    """`1 - d (1 - exp(-t / tau_w))`: the share of the fresh albedo a surface of this age keeps (1 for a fresh surface)."""
    if surface_age_gy < 0.0 or not tau_w_gy > 0.0 or not 0.0 <= depth < 1.0:
        raise ValueError("the age must not be negative, the time positive and the amplitude in [0, 1)")
    return 1.0 - depth * (1.0 - math.exp(-surface_age_gy / tau_w_gy))


def draw_weathering():
    """`(tau_w_gy, depth)` for one body: the time log-uniform in `TAU_W_GY`, then the amplitude uniform in `DEPTH`."""
    tau = math.exp(draw.uniform(math.log(TAU_W_GY[0]), math.log(TAU_W_GY[1])))
    return tau, draw.uniform(*DEPTH)


def weathered_albedo(fresh, age_gy):
    """`fresh` after `age_gy` of space weathering, with the time and the amplitude drawn for the body (`draw_weathering`)."""
    tau, depth = draw_weathering()
    return fresh * weathering_factor(age_gy, tau, depth)
