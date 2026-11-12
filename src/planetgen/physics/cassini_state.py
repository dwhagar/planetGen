# planetgen/physics/cassini_state.py

"""
The tilt of a locked or resonant planet from its real precession constant (GEN.1536)
===================================================================================

A planet that turns once an orbit (locked) or in a resonance such as Mercury's 3:2 has its spin axis in the first **Cassini state**: the axis precesses about the orbit's normal at the
precession constant `V`, the normal precesses about the Laplace plane at `s` (negative: nodal regression), and the axis settles in the plane of the two with an obliquity

    eps = i |s| / (V + |s|)          (small `i`, the inclination of the orbit to the Laplace plane)

which is 0 where the planet's shape holds the axis to the orbit (`V >> |s|`) and `i` where it does not. The planet's orbit carries several nodal modes (`obliquity_variation`), each of
them drives the axis on its own in this linear theory, and the tilt is the quadrature sum over the modes `k`: `eps^2 = sum_k (i_k |s_k| / (V + |s_k|))^2`.

**The spin side, `V`.** Peale's form for a body in the p:2 resonance with its orbit (as [R5229] write it: `eps = -c Omega' sin(i) / (c Omega' cos(i) + 2 n C22 (7/2 e - 123/16 e^3) / ... )`), divided
by the body's moment `C`, gives

    V = (n / c) [ J2 (1 - e^2)^(-3/2) + 2 C22 (7/2 e - 123/16 e^3) ]

for the orbital mean motion `n`, `c = C / M R^2`, and the body's **stored figure** `J2` and `C22`. Mercury's MESSENGER figure (`J2` 5.03e-5, `C22` 8.09e-6, `c` 0.346, [R5229]) gives
`V` of 1000 arcsec per year against its nodal regression of 4.7, and with the inclination of 8.6 degrees a tilt of 2.4 arcmin against the 2.04 measured. (For a body locked 1:1 the `C22` term
is taken as 0: its coefficient is of the order of the eccentricity, which a locked body has tidally worn down.)

**The figure.** The first version of GEN.1536's idea took the fluid Darwin-Radau figure `J2 = k_f q / 3`, which gives a Mercury-like body a `V` of 5 arcsec per year and a tilt of 4 to 7
degrees against the real 0.03: the real figure of a body that has despun is not the fluid figure of its present spin. Table 2 of [R5228] sets each body's observed `J2` against the
hydrostatic one it would have (the rotational part `k_f q_omega / 3` and the tidal part `(k_f / 2) q_n (1 - e^2)^(-3/2)`, `q_x = x^2 R^3 / G M`):

| body | observed `J2` / hydrostatic | state |
|---|---|---|
| Io, Europa, Ganymede, Titan, Rhea, Callisto | 0.95 to 1.02 | locked, icy or tidally hot: the figure follows the spin |
| Earth, Mars | 1.0 to 1.07 | free, not despun |
| Moon | 22.7 | locked: a fossil bulge from a much faster spin |
| Venus | 86 | slow, despun: fossil |
| Mercury | 94 | 3:2: fossil |

so a body whose interior can relax (an icy body, or a rocky one heated above `RELAXED_HEAT_FLUX_W_M2`) has the hydrostatic figure, and a cold rocky body that has despun has a figure frozen
in at an earlier, faster spin: `FOSSIL_FACTOR_RANGE` times the hydrostatic `J2` (drawn log-uniform, a stand-in between the Moon and Mercury), with the observed share `FOSSIL_C22_OVER_J2` of it in
`C22` (Mercury 0.161, the Moon 0.110) in place of the hydrostatic 0.3. The figure is derived, not stored in a column: it is a function of the body and one draw from a stream of its own.

**What it does not do.** Moons (`tidal_spin.cassini_tilt_deg` still stands for them: their nodal regression comes from the planet's bulge and the star, not from planet-planet modes, GEN.1536
part 2); a Cassini state 2 (the Moon's own, reached through a different tidal path); the exact, not small-angle, state; a body in a resonance above 3:2 uses the 3:2 coefficient.
"""

import math

from planetgen.physics import bulk_draw, constants, obliquity_variation, tidal_spin

FOSSIL_FACTOR_RANGE = (22.7, 94.2)
"""tuple: The observed `J2` over the hydrostatic one of a cold rocky body that has despun, drawn log-uniform: the Moon's 203.8e-6 over 8.968e-6 and Mercury's 50.3e-6 over 0.534e-6 (Table 2 of [R5228]; Venus
between, 86). Source tier T3 (three bodies), tweak level L1 (the bracket is the data; the log-uniform spread inside it is the stand-in)."""
FOSSIL_C22_OVER_J2 = 0.133
"""float: `C22 / J2` of a fossil figure: the geometric mean of Mercury's 0.161 and the Moon's 0.110 [R5228: Table 2], against the hydrostatic 3/10. Source tier T2, tweak level L0."""
HYDROSTATIC_C22_OVER_J2 = 0.3
"""float: `C22 / J2` of a synchronously turning hydrostatic body, 3/10 [R5228]. Source tier T1, tweak level L0."""
RELAXED_HEAT_FLUX_W_M2 = 0.5
"""float: The internal heat flux above which a rocky body keeps its interior relaxed, so that its figure follows its spin (Io, 2.5; Venus and Mercury, under 0.1, are fossil). Source tier T3, tweak level L1."""


def hydrostatic_figure(mass_kg, radius_m, moment_factor, spin_rate_rad_s, mean_motion_rad_s, eccentricity):
    """`(J2, C22)` of the fluid figure of a body that turns at `spin_rate_rad_s` about a primary that pulls it with `mean_motion_rad_s`: `k_f q_omega / 3 + (k_f / 2) q_n (1 - e^2)^(-3/2)` for `J2`, and 0.3 of it for `C22`."""
    if not (mass_kg > 0.0 and radius_m > 0.0 and 0.0 <= eccentricity < 1.0):
        raise ValueError("mass and radius must be positive and the eccentricity between 0 and 1")
    k_f = obliquity_variation.fluid_love_number(moment_factor)
    scale = radius_m ** 3 / (constants.G * mass_kg)
    j2 = k_f * (spin_rate_rad_s ** 2 * scale / 3.0 + 0.5 * mean_motion_rad_s ** 2 * scale * (1.0 - eccentricity ** 2) ** -1.5)
    return j2, HYDROSTATIC_C22_OVER_J2 * j2


def frozen_figure(hydrostatic_j2, relaxed, u):
    """
    `(J2, C22)` the body keeps: the hydrostatic one where its interior relaxes (`relaxed`), else the fossil one, `FOSSIL_FACTOR_RANGE` times the hydrostatic `J2` log-uniform for the unit draw `u`,
    with `FOSSIL_C22_OVER_J2` of it in `C22`.
    """
    if relaxed:
        return hydrostatic_j2, HYDROSTATIC_C22_OVER_J2 * hydrostatic_j2
    low, high = FOSSIL_FACTOR_RANGE
    j2 = hydrostatic_j2 * low * (high / low) ** u
    return j2, FOSSIL_C22_OVER_J2 * j2


def precession_rad_yr(j2, c22, moment_factor, mean_motion_rad_s, eccentricity, resonant):
    """`V = (n / c) [J2 (1 - e^2)^(-3/2) + 2 C22 (7/2 e - 123/16 e^3)]`, radians per year; the `C22` term only for a resonant (not a 1:1) body."""
    figure = j2 * (1.0 - eccentricity ** 2) ** -1.5
    if resonant:
        figure += 2.0 * c22 * (3.5 * eccentricity - 123.0 / 16.0 * eccentricity ** 3)
    return max(0.0, mean_motion_rad_s / moment_factor * figure * constants.SECONDS_PER_YEAR)


def tilt_deg(modes, precession):
    """
    The obliquity of a body in the first Cassini state, degrees: the quadrature sum over the nodal modes of `i_k |s_k| / (V + |s_k|)`, each mode with its frequency `s_k` (radians per year) and amplitude
    `i_k` (radians) as `obliquity_variation.Mode` has them. 0 for a body with no mode (an orbit in the plane).
    """
    total = 0.0
    for mode in modes:
        speed = abs(mode.frequency_rad_yr)
        if speed <= obliquity_variation.MODE_FLOOR_RAD_YR:
            continue
        total += (mode.amplitude_rad * speed / (precession + speed)) ** 2
    return math.degrees(math.sqrt(total))


def is_relaxed(body, heat_flux_w_m2):
    """Whether the interior of the body keeps the figure of its spin: an icy body, or a rocky one heated above `RELAXED_HEAT_FLUX_W_M2`."""
    return bulk_draw.origin_of(body) == bulk_draw.ICY or heat_flux_w_m2 >= RELAXED_HEAT_FLUX_W_M2


def locked_tilt_deg(body, heat_flux_w_m2, modes, mean_motion_rad_s, turns_per_orbit, u):
    """
    The Cassini-state tilt of a locked or resonant `body` (mass, radius in kilometres, eccentricity), degrees, for the nodal `modes` of its orbit, its mean motion, the spin it turns at
    (`turns_per_orbit` of them each orbit) and the unit draw `u` of the fossil figure.
    """
    moment_factor = tidal_spin.moment_factor(body)
    eccentricity = getattr(body, "eccentricity", 0.0) or 0.0
    spin_rate = turns_per_orbit * mean_motion_rad_s
    j2_hydro, _ = hydrostatic_figure(body.mass, body.radius * 1000.0, moment_factor, spin_rate, mean_motion_rad_s, eccentricity)
    j2, c22 = frozen_figure(j2_hydro, is_relaxed(body, heat_flux_w_m2), u)
    precession = precession_rad_yr(j2, c22, moment_factor, mean_motion_rad_s, eccentricity, turns_per_orbit > 1.0)
    inclination = abs(getattr(body, "orbital_inclination_deg", 0.0) or 0.0)
    ceiling = min(inclination, 180.0 - inclination)      # `x` is at most 1 (GEN.1165): the quadrature sum of the modes' amplitudes can pass the orbit's inclination of this moment, the tilt does not
    return min(tilt_deg(modes, precession), ceiling)
