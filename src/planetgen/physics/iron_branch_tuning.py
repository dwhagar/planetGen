# planetgen/physics/iron_branch_tuning.py

"""
The simulated sample the iron branch is tuned to (GEN.512, docs/design/composition/iron-branch-tuning.md)
==========================================================================================================

The published statistic is how many of a homogeneous sample of small rocky planets are more iron-rich than their star, not what the planets
are made of: each planet's core mass fraction (CMF) is inferred from a measured mass and radius, and called iron-rich when it exceeds the
star-matched value by 2 sigma of that inference ([R3546]: 1 of 21; [R0541]: the share within 1 sigma of the star, 15 of 20).
The sampler here draws what the harness needs on the same footing:

1. a host, with its abundances from the generator's own star model (`bulk_draw`'s ratios at a drawn [Fe/H]), and the dry partition
   that star-matched iron gives (the star's CMF, `bulk_draw.dry_partition`);
2. the planet's branch (`bulk_draw.branch_ratios`: scatter about the star, or enriched) and its true CMF;
3. a mass log-uniform over `MASS_RANGE_EARTH`, a period log-uniform over `PERIOD_RANGE_D`, the radius the true CMF gives
   (`bulk_structure.radius_fit_earth`); a planet above `MAX_RADIUS_EARTH` is not a rocky sample member and is dropped;
4. the measurement the sample would make (`mass_radius_sample.observe`: its mass and radius errors) and the CMF inferred from it;
5. `sigma`: the CMF inference's 1-sigma error from the two errors (the CMF at the measured radius moved by its error, and the same for
   the mass, added in quadrature); `outlier` (the measured CMF is more than 2 sigma above the star-matched value) and `within`
   (it is within 1 sigma of it).

The precisely-measured cut is the `mass-radius-sample` selection model, so the sampler returns the planet with its measurement and the harness
keeps the planets that pass. The planets are not any the generator keeps.
"""

import math

from planetgen.physics import bulk_draw, bulk_structure, mass_radius_sample, stellar_composition
from planetgen.util import detmath, draw

FE_H_SIGMA = 0.2
"""float: The spread of the hosts' [Fe/H] (normal about 0); a host sample of FGK stars (the homogeneous samples are 0.1 to 0.3 dex wide); the statistic
is relative to each host, so the width matters little."""
MASS_RANGE_EARTH = (1.0, 10.0)
PERIOD_RANGE_D = (0.5, 30.0)
MAX_RADIUS_EARTH = 1.8
"""float: The samples' rocky cut ([R0541]: 1 to 1.8 Earth radii)."""
HOST_MASS_RANGE_SOLAR = (0.8, 1.2)
OUTLIER_SIGMAS = 2.0
FEO = 0.08
LIGHT = 0.07
"""float: The mantle FeO and core light-element shares the inference assumes (mid-range of `bulk_draw.RANGES`)."""


def _inferred_cmf(mass_earth, radius_earth, x_fe):
    cmf = bulk_structure.cmf_from_mass_radius(mass_earth, radius_earth, x_fe, 1.0, high=bulk_draw.MAX_CMF)
    if cmf is None:
        cmf = bulk_draw.MAX_CMF if radius_earth < bulk_structure.radius_fit_earth(mass_earth, bulk_draw.MAX_CMF, x_fe, 1.0) else 0.0
    return cmf


def planet(enriched_weight, scatter_dex):
    """One planet and the host it orbits: the quantities `mass_radius_sample.observe` needs, plus `true_cmf` and `star_cmf`; None above the rocky cut."""
    fe_h = draw.gauss(0.0, FE_H_SIGMA)
    host_mass = math.exp(draw.uniform(math.log(HOST_MASS_RANGE_SOLAR[0]), math.log(HOST_MASS_RANGE_SOLAR[1])))
    key = stellar_composition.abundance_key(mass=host_mass, distance_ly=draw.uniform(20000.0, 30000.0))
    ratios = bulk_draw._ratios(fe_h, 4.0, key)
    star_part = bulk_draw.dry_partition(ratios, FEO, LIGHT, False)
    branch = "iron-enriched" if draw.uniform(0.0, 1.0) < enriched_weight else "star-matched"
    mine = bulk_draw.branch_ratios(ratios, branch, scatter_dex, FEO, LIGHT, False, 0.0)
    true_cmf = bulk_draw.dry_partition(mine, FEO, LIGHT, False).metal
    mass = math.exp(draw.uniform(math.log(MASS_RANGE_EARTH[0]), math.log(MASS_RANGE_EARTH[1])))
    period = math.exp(draw.uniform(math.log(PERIOD_RANGE_D[0]), math.log(PERIOD_RANGE_D[1])))
    x_fe = bulk_structure.x_fe_from_feo_mass_fraction(star_part.rock_feo)
    radius = bulk_structure.radius_fit_earth(mass, min(true_cmf, 0.8), x_fe, 1.0)
    if radius > MAX_RADIUS_EARTH:
        return None
    return {"mass_earth": mass, "radius_earth": radius, "period_d": period, "host_mass_solar": host_mass,
            "host_radius_solar": host_mass ** 0.9, "true_cmf": true_cmf, "star_cmf": star_part.metal, "x_fe": x_fe,
            "branch": branch}


def classify(item):
    """`item` (a measured planet from `mass_radius_sample.observed_sampler`) with `sigma`, `outlier` and `within`."""
    mass, radius, x_fe = item["measured_mass_earth"], item["measured_radius_earth"], item["x_fe"]
    cmf = _inferred_cmf(mass, radius, x_fe)
    by_radius = abs(_inferred_cmf(mass, radius * (1.0 + item["radius_error_share"]), x_fe) - cmf)
    by_mass = abs(_inferred_cmf(mass * (1.0 + item["mass_error_share"]), radius, x_fe) - cmf)
    sigma = max(detmath.hypot(by_radius, by_mass), 1e-6)
    excess = cmf - item["star_cmf"]
    return dict(item, measured_cmf=cmf, sigma=sigma, outlier=excess > OUTLIER_SIGMAS * sigma, within=abs(excess) < sigma)


def sampler(enriched_weight, scatter_dex):
    """A sampler for `tuning_harness.simulate`: a measured planet with its classification, or None above the rocky cut."""
    observed = mass_radius_sample.observed_sampler(lambda _i: planet(enriched_weight, scatter_dex))

    def sample(i):
        item = observed(i)
        return None if item is None else classify(item)
    return sample
