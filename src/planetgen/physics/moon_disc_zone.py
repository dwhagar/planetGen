# planetgen/physics/moon_disc_zone.py

"""
Where the regular moons of a giant form (GEN.1244)
==================================================

The regular moons of a giant planet are made in its circumplanetary disc, from the gas and the solids that fall into it from the planet's Hill sphere. Two edges bound the zone they form in:

* **The outer edge, the centrifugal radius.** Infalling gas lands where a circular orbit has the specific angular momentum of the inflow's average, the centrifugal radius `r_c`, and the infall of gas and of solids is confined inside it (the flux is constant at `r < r_c` and zero beyond;
  solids "remain zero for r > r_c", [R3922] Appendix A). The paper adopts `r_c = 30` planet radii from three-dimensional hydrodynamic runs of the flow within the Hill sphere, the same for Jupiter and Saturn, and the four large moons it
  forms sit at 15 to 26 planet radii (Table 1). The moons are made inside `r_c`; they are not made at tens of the planet's radii past it, where the old draw put them (46 to 130 radii, outside the radiation belts that reach about 20).
* **The inner edge, the corotation radius.** Magnetic coupling truncates the disc where its gas corotates with the planet's spin ([R3922] section 2.4, the Takata and Stevenson picture): `r_co = (G M / omega^2)^(1/3)`, 2.25 planet radii for the Jupiter
  of that paper. A planet that spins so slowly that this lies outside `r_c` has no cavity (the disc reaches in to the planet), and one that spins so fast that it lies inside the body's own reach has the usual clearance (`planets.moon_orbit_bounds_km`).

`zone_km` gives the two edges of a giant, bounded by the orbits a moon can keep (`low_km`, `high_km`); an empty zone has `inner >= outer`. The disc's own diffusion reaches about 150 planet radii but carries no solids there, so it is not an edge for moons.
Stand-ins are MD1 and MD2 of `docs/design/research/gap-filling-plan.md` section 6.167.
"""

import math

from planetgen.physics import constants

CENTRIFUGAL_RADII = 30.0
"""float: The centrifugal radius of the disc, in planet radii: where the infalling gas and solids land, and so the outer edge of the zone the regular moons form in ([R3922] Appendix A, adopted from three-dimensional hydrodynamic runs;
the same for every giant). Source tier T1 (the paper's own value, not re-derived), tweak level L0."""


def corotation_km(mass_kg, spin_rate_rad_s):
    """The corotation radius, km from the planet's centre: `(G M / omega^2)^(1/3)`; `math.inf` for a planet that does not spin."""
    omega = abs(spin_rate_rad_s or 0.0)
    if omega <= 0.0:
        return math.inf
    return (constants.G * mass_kg / omega ** 2) ** (1.0 / 3.0) / 1000.0


def zone_km(mass_kg, radius_km, spin_rate_rad_s, low_km, high_km):
    """
    The zone in which a giant's disc moons form: `(inner_km, outer_km)` from the planet's centre. The outer edge is the centrifugal radius, `CENTRIFUGAL_RADII` planet radii; the inner edge is the corotation radius of the
    birth spin when that lies inside the outer edge, else the disc reaches in to the clearance `low_km`. Both are kept inside the orbits a moon can hold (`low_km` to `high_km`, from `planets.moon_orbit_bounds_km`).
    """
    outer = min(high_km, CENTRIFUGAL_RADII * radius_km)
    corotation = corotation_km(mass_kg, spin_rate_rad_s)
    inner = low_km if corotation >= outer else max(low_km, corotation)
    return inner, outer
