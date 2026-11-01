# planetgen/physics/belt_mass.py

"""
How much is left in a planetesimal belt: the disc's solids ground down by collisions (GEN.477 E2, part 1; docs/design/composition/small-bodies-composition.md section 14).

The belt density of `AsteroidBelt` was a uniform draw of "dense", "sparse" or "typical" and the belt's mass was stated as zero (NAV.62). Here the mass is derived.

1. **What formed.** The solids of the disc between the belt's limits, `planet_mass.solids_surface_density_gcm2` integrated over the annulus: the primordial mass,
   a few Earth masses at 2 to 3 AU for a Sun-like star.
2. **What collisions leave.** [R4785] (MNRAS 380, 1737, [R4785], section 2.2, READ) work out a planetesimal belt in collisional
   equilibrium: the mass falls as `M(t) = M0 / (1 + t / tc(0))` (their Eq. 14) and `tc` depends on the mass, so the late-time mass is independent of the initial mass
   (their Eq. 13 and 15): `M_max = 3.8e-6 rho r^3.5 (dr/r) Dc / (M*^0.5 t) * (12q - 20) [1 + 1.25 (e/I)^2]^-0.5 / ((18 - 9q) G(q, Xc))`, with
   `G(11/6, Xc) = Xc^-0.5 + 0.67 Xc^-1.5 + 0.2 Xc^-2.5 - 1.87`, `Xc = 1e-3 (r Q*D / e^2)^(1/3)`, `q = 11/6`, in Earth masses for `rho` in kg/m^3, `r` and `dr` in AU,
   `Dc` (the largest planetesimal, in equilibrium) in km, `M*` in solar masses and `t` in Myr. `belt_mass_earth` is `1 / (1 / M0 + t / A)` with `A = M_max t`, so a young or
   light belt keeps what formed and an old belt keeps `M_max` whatever formed: in a belt a few AU wide, the disc's thousands of times the mass.
3. **The numbers it needs from the belt.** `Dc` and the disruption threshold `Q*D` of the largest bodies: for a gravity-bound body the generator's own
   [R2633] threshold of `collision_outcome.catastrophic_energy`, with the strength of 200 J/kg of Wyatt et al.'s nominal 2 km planetesimals as a
   floor; the particles' eccentricity from the belt's width, `e = (upper - lower) / (upper + lower)` (the extent of the orbits of a belt that stirs itself), and `e/I = 1`.

The Solar System's belt is a test and never an input: at 4.5 Gyr, 2.2 to 3.4 AU, with Ceres-size bodies, the model leaves 0.3 to 3 thousandths of an Earth mass; the
inferred mass is 0.6 thousandths ([R2471], as Wyatt et al. quote it, whose own estimate for their belt is 0.4 thousandths).

Limits. One size, `Dc`, stands for the whole top of the cascade; `M_max` goes as `Dc` to a power near 2.3 here (the strength grows with size), so the mass is only as good as the largest
body. `Dc` is each belt's own: `belt_ledger.belt_state` draws the largest body of the belt from the size law and settles it with this mass (`belt_ledger.solve_belt`, GEN.1061 and GEN.1000;
no stand-in largest body is left). A belt stirred by a giant has a larger `e` than its width gives. The grains lost to radiation forces and the belt's early mass loss by sweeping are not
in the law.

Wired (GEN.477, GEN.1000): `AsteroidBelt` takes its mass (`asteroid_belts.mass_kg`, a stated zero before) and its density word from `belt_ledger.belt_state`, which replaced the uniform draw of the
word. The belt's composition thresholds still take the Dirichlet draw's shares; the mass is not yet an input there.
"""

import math

from planetgen.physics import collision_outcome, constants, planet_mass
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")
_mearth_to_kg = units.converter("Mearth", "kg")
from planetgen.util.checks import finite_domain

GRAIN_DENSITY_KG_M3 = 2700.0
"""float: The planetesimals' bulk density in Wyatt et al.'s model (their section 2.3, the nominal 2700 kg/m^3). [Stand-in BM1, tier T1, tweak L0.]"""

STRENGTH_FLOOR_J_PER_KG = 200.0
"""float: The disruption threshold of the km-size planetesimals Wyatt et al. take as nominal (their Fig. 2), a floor under the gravity-bound threshold. [Stand-in BM1, tier T1, tweak L0.]"""

DENSE_FACTOR = 3.0
"""float: A belt with this many times the reference mass per area is called dense and one with its inverse sparse; between them it is typical. [Stand-in BM2, tier T9, tweak L4.]"""

_Q = 11.0 / 6.0
_SQRT_COEFFICIENT = (12.0 * _Q - 20.0) / (18.0 - 9.0 * _Q)     # 2 / 1.5 for the size law q = 11/6


def _g(xc):
    return xc ** -0.5 + 0.67 * xc ** -1.5 + 0.2 * xc ** -2.5 - 1.87


@finite_domain()
def disruption_threshold_j_per_kg(largest_km, density_kg_m3=GRAIN_DENSITY_KG_M3):
    """`Q*D` of the largest planetesimals: the Leinhardt and Stewart threshold of two bodies of that diameter, or `STRENGTH_FLOOR_J_PER_KG` if that is larger."""
    if not largest_km > 0.0:
        raise ValueError("the largest planetesimal must have a positive size")
    mass = 4.0 / 3.0 * math.pi * (largest_km * 500.0) ** 3 * density_kg_m3
    return max(collision_outcome.catastrophic_energy(mass, mass), STRENGTH_FLOOR_J_PER_KG)


@finite_domain()
def primordial_mass_earth(lower_au, upper_au, star_mass_kg, snow_line_au, metal_ratio=1.0, scatter=1.0, steps=64):
    """The solids of the disc between two distances, Earth masses: `2 pi a Sigma(a)` integrated over the annulus (log-spaced steps)."""
    if not 0.0 < lower_au < upper_au:
        raise ValueError("the belt needs 0 < lower < upper")
    au_cm = constants.AU_M * 100.0
    total_g, log_low, log_step = 0.0, math.log(lower_au), (math.log(upper_au) - math.log(lower_au)) / steps
    for i in range(steps):
        a_low, a_high = math.exp(log_low + i * log_step), math.exp(log_low + (i + 1) * log_step)
        a_mid = math.sqrt(a_low * a_high)
        sigma = planet_mass.solids_surface_density_gcm2(a_mid, star_mass_kg, snow_line_au, metal_ratio, scatter)
        total_g += math.pi * (a_high ** 2 - a_low ** 2) * au_cm ** 2 * sigma
    return total_g / _mearth_to_kg(1000.0)


@finite_domain()
def collisional_constant(lower_au, upper_au, star_mass_solar, largest_km, density_kg_m3=GRAIN_DENSITY_KG_M3):
    """`A = M_max t`, Earth masses times Myr, of a belt between two distances whose largest planetesimal is `largest_km` across ([R4785: Eq. 15])."""
    r = 0.5 * (lower_au + upper_au)
    dr = upper_au - lower_au
    e = (upper_au - lower_au) / (upper_au + lower_au)
    strength = disruption_threshold_j_per_kg(largest_km, density_kg_m3)
    xc = 1.0e-3 * (r * strength / e ** 2) ** (1.0 / 3.0)
    return (3.8e-6 * density_kg_m3 * r ** 3.5 * (dr / r) * largest_km / star_mass_solar ** 0.5
            * _SQRT_COEFFICIENT * (1.0 + 1.25) ** -0.5 / _g(xc))


@finite_domain()
def belt_mass_earth(primordial_earth, lower_au, upper_au, star_mass_solar, age_myr, largest_km):
    """What is left of a belt of `primordial_earth` Earth masses after `age_myr`: `1 / (1 / M0 + t / A)` ([R4785: Eq. 14] with the mass-dependent lifetime); all of it if `Dc` is so large that `A` has no meaning (`Xc` over 1)."""
    if not (primordial_earth > 0.0 and age_myr >= 0.0):
        raise ValueError("the primordial mass must be positive and the age zero or more")
    a = collisional_constant(lower_au, upper_au, star_mass_solar, largest_km)
    if a <= 0.0:       # G(Xc) <= 0 past Xc = 1: the largest bodies are too strong for the belt's collisions to disrupt, the cascade never starts and the belt keeps what formed (GEN.1061)
        return primordial_earth
    return 1.0 / (1.0 / primordial_earth + age_myr / a)


@finite_domain()
def surface_density_earth_per_au2(mass_earth, lower_au, upper_au):
    """The belt's mass over its area, Earth masses per square AU."""
    return mass_earth / (math.pi * (upper_au ** 2 - lower_au ** 2))


def density_word(mass_per_area, reference_mass_per_area):
    """`dense`, `sparse` or `typical`: the belt's mass per area against a reference (`belt_ledger.belt_state` gives the median belt of the same age), `DENSE_FACTOR` each way."""
    if not (mass_per_area > 0.0 and reference_mass_per_area > 0.0):
        raise ValueError("both surface densities must be positive")
    ratio = mass_per_area / reference_mass_per_area
    if ratio >= DENSE_FACTOR:
        return "dense"
    return "sparse" if ratio <= 1.0 / DENSE_FACTOR else "typical"
