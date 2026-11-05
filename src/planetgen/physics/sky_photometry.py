# planetgen/physics/sky_photometry.py

"""
The photometry of the system sky, on arrays (VIEW.20): how bright a star looks from another system, and how much
of its light the dust of the galaxy takes on the way.

This is the lighter 8.2 version of the photometry library (VIEW.11) and the dust model (VIEW.8): the sky catalogue
(`galaxy/sky.py`, VIEW.23) and the glow band (VIEW.27) read it. Reddening, the desaturation of faint stars and the
dust of nebulae stay with VIEW.11 and VIEW.8 in 9.1, which can reuse these functions.

* **Apparent magnitude** `m = M_bol,sun - 2.5 log10(L / L_sun) - BC_V(T) + 5 log10(d / 10 pc) + A_V`
  (sky-view.md 2.4): `absolute_visual_magnitude` and `apparent_magnitude`, on numpy arrays. L is the bolometric
  luminosity the generator stores (`luminosity_w` over the Sun's), not a V-band one.
* **The bolometric correction for any temperature** (`bolometric_correction`):
  - 4,100 K to 50,000 K: Flower's (1996) fit with Torres's (2010) coefficients, the same polynomials as
    `photometry.bolometric_correction_v` ([R4407: Table 1], VIEW.9), on arrays.
  - **Under 4,100 K (the cool-star rule VIEW.9 left open):** the mean main-sequence table of Pecaut and Mamajek
    ([R3469: Table 5]; E. Mamajek's table version 2022.04.16), from K7V (4,100 K, -0.93) to L2V (2,060 K,
    -6.62), linear in log T between rows and held at the last row below 2,060 K. Flower's fit "breaks down
    completely for the M dwarfs" [R4407] and is 0.6 mag too faint at 3,500 K against the table (-2.31
    against -1.7); the table is the better source from where it starts, and it meets the fit within 0.1 mag at
    4,100 K. The table is for dwarfs: a cool giant (a larger, brighter star of the same temperature) has a smaller
    correction in the red, so its V magnitude here is somewhat too faint; a stand-in, registered as SKP2.
  - **Over 50,000 K:** the blackbody tail. A star's bolometric flux goes as T^4 and its flux at the V wavelength
    (550 nm) as the Planck function there, which for a hot star is the Rayleigh-Jeans side, proportional to T, so the
    correction falls by 7.5 magnitudes a decade of temperature in the limit. The step from the fit's last value is
    `-2.5 log10[(T^4 expm1(x)) / (T0^4 expm1(x0))]` with `x = c2 / (550 nm T)` (7.2 magnitudes a decade from 50,000 K, 7.5
    at a million K). A hot white dwarf, a Wolf-Rayet star or a stripped helium star is faint in V for its
    luminosity. Stand-in SKP3.
* **Galactic dust** (`galactic_extinction_av`): the smooth dust disc of [R1158],
  `rho(R, z) = rho_sun exp(-(R - R_sun) / h_R) sech^2(z / h_z(R))`, an exponential disc of scale length
  `h_R` = 2.26 kpc, a sech^2 layer whose scale height is 134 pc inside 4.4 kpc and rises 14.8 pc per kpc beyond (188 pc at
  the solar circle), from their Eq. 13 and 14 and Table 1 (the paper's own numbers; the farthest-out values read from its
  summary of Table 1, to confirm), `R_sun` = 8 kpc. Their normalisation is far-infrared, not A_V, so `rho_sun` is set
  so that the plane at the solar circle dims by 1.0 mag per kpc (the plan's stand-in, SKP1; the literature has 0.7 to
  1.0 near the Sun and about 1.8 in the plane over a few kpc), and the line of sight is integrated numerically. The
  paper's warp, spiral arms and central hole are left out (stars through the inner 4 kpc are dimmed too much); the
  caller scales the lengths for a galaxy of another size.

Nothing here is stored: a magnitude is worked out from stored quantities when a sky is built.
"""

import numpy as np

from planetgen.physics import photometry

COOL_BC_TABLE_K_MAG = (
    (2060.0, -6.62), (2160.0, -6.48), (2270.0, -6.25), (2350.0, -6.13), (2380.0, -5.86), (2420.0, -5.78),
    (2570.0, -5.65), (2630.0, -5.32), (2680.0, -4.99), (2740.0, -4.62), (2810.0, -4.13), (2930.0, -3.58),
    (3060.0, -3.11), (3110.0, -2.84), (3210.0, -2.51), (3270.0, -2.28), (3430.0, -1.93), (3470.0, -1.78),
    (3560.0, -1.62), (3620.0, -1.50), (3660.0, -1.42), (3770.0, -1.29), (3850.0, -1.15), (3930.0, -1.07),
    (3990.0, -1.03), (4100.0, -0.93),
)
"""tuple: `(effective temperature K, BC_V mag)` of the main-sequence rows K7V to L2V of Pecaut and Mamajek ([R3469: Table 5]), as
E. Mamajek's table of 2022.04.16 gives them. Same zero point as Flower's fit (the G2V row,
5,770 K, is -0.085 against the fit's -0.08 at 5,772 K)."""

COOL_BELOW_K = COOL_BC_TABLE_K_MAG[-1][0]
"""float: Under this temperature the table answers, from it up the fit."""

FIT_TOP_K = photometry.BOLOMETRIC_CORRECTION_RANGE_K[1]
"""float: Over this temperature the Rayleigh-Jeans tail answers."""

V_WAVELENGTH_M = 550e-9
"""float: The wavelength, m, the hot tail takes the V band at."""

DUST_SOLAR_RADIUS_PC = 8000.0
DUST_SCALE_LENGTH_PC = 2260.0
DUST_HEIGHT_PC = 134.4
DUST_HEIGHT_FLARE_PC_PER_KPC = 14.8
DUST_FLARE_START_PC = 4400.0
"""float: The solar radius, the dust disc's scale length, its scale height at the centre of the flare and the flare, from
[R1158]: `h_z(R) = 134.4 pc + 14.8 pc/kpc (R - 4.4 kpc)` for R over 4.4 kpc."""

DUST_AV_PER_KPC_AT_SUN = 1.0
"""float: Magnitudes of V extinction per kpc in the plane at the solar circle: the stand-in of the sky plan (SKP1)."""

EXTINCTION_STEPS = 512
"""int: Samples along a line of sight. The sech^2 layer is 134 pc thick, so 60 pc steps at 30 kpc still sample it."""


def _log10(value):
    return np.log10(np.asarray(value, dtype=float))


def _flower(log_t):
    """Flower's fit on an array of log10 T already inside its range."""
    out = np.empty_like(log_t)
    start = -np.inf
    for upper, coefficients in photometry.FLOWER_COEFFICIENTS:
        mask = (log_t >= start) & (log_t < upper)
        if mask.any():
            x = log_t[mask]
            out[mask] = sum(coefficient * x ** power for power, coefficient in enumerate(coefficients))
        start = upper
    return out


def _planck_term(temperature_k):
    """`log10(T^4 expm1(c2 / (lambda_V T)))`: how much bolometric flux there is per unit of flux at the V wavelength, up to a constant."""
    t = np.asarray(temperature_k, dtype=float)
    x = photometry.SECOND_RADIATION_CONSTANT_M_K / (V_WAVELENGTH_M * t)
    return 4.0 * np.log10(t) + np.log10(np.expm1(x))


def bolometric_correction(temperature_k):
    """
    `BC_V` (magnitudes, `M_V = M_bol - BC_V`) of stars of effective temperature `temperature_k` (an array or a number), for
    every temperature: Flower's fit from `COOL_BELOW_K` to `FIT_TOP_K`, the main-sequence table under it, the
    blackbody tail over it (the module text).

    Raises:
        ValueError: A temperature that is not positive and finite.
    """
    t = np.asarray(temperature_k, dtype=float)
    if not np.all(np.isfinite(t)) or np.any(t <= 0):
        raise ValueError("a temperature must be positive and finite")
    out = np.empty_like(t)
    cool, hot = t < COOL_BELOW_K, t > FIT_TOP_K
    middle = ~(cool | hot)
    if cool.any():
        rows = np.array(COOL_BC_TABLE_K_MAG)
        out[cool] = np.interp(np.log10(t[cool]), np.log10(rows[:, 0]), rows[:, 1])
    if middle.any():
        out[middle] = _flower(np.log10(t[middle]))
    if hot.any():
        top = float(_flower(np.array([np.log10(FIT_TOP_K)]))[0])
        out[hot] = top - 2.5 * (_planck_term(t[hot]) - _planck_term(FIT_TOP_K))
    return out if out.ndim else float(out)


def absolute_visual_magnitude(luminosity_sol, temperature_k):
    """`M_V = M_bol,sun - 2.5 log10(L / L_sun) - BC_V(T)` for bolometric luminosities in solar units."""
    lum = np.asarray(luminosity_sol, dtype=float)
    if np.any(lum <= 0) or not np.all(np.isfinite(lum)):
        raise ValueError("a luminosity must be positive and finite")
    return photometry.SUN_ABSOLUTE_BOLOMETRIC_MAG - 2.5 * np.log10(lum) - bolometric_correction(temperature_k)


def apparent_magnitude(luminosity_sol, temperature_k, distance_pc, a_v=0.0):
    """
    The V magnitude of stars seen from `distance_pc` (parsecs) through `a_v` magnitudes of dust:
    `M_V + 5 log10(d / 10 pc) + A_V`. Arrays or numbers, broadcast together.

    Raises:
        ValueError: A distance, luminosity or temperature that is not positive and finite.
    """
    d = np.asarray(distance_pc, dtype=float)
    if np.any(d <= 0) or not np.all(np.isfinite(d)):
        raise ValueError("a distance must be positive and finite")
    return absolute_visual_magnitude(luminosity_sol, temperature_k) + 5.0 * np.log10(d / 10.0) + a_v


def dust_av_per_pc(x_pc, y_pc, z_pc, *, av_per_kpc_at_sun=DUST_AV_PER_KPC_AT_SUN, scale_length_pc=DUST_SCALE_LENGTH_PC,
                   height_pc=DUST_HEIGHT_PC, flare_pc_per_kpc=DUST_HEIGHT_FLARE_PC_PER_KPC,
                   flare_start_pc=DUST_FLARE_START_PC, solar_radius_pc=DUST_SOLAR_RADIUS_PC):
    """The dust disc's extinction in magnitudes per parsec at galactocentric `(x, y, z)` in parsecs (arrays or numbers)."""
    r = np.hypot(np.asarray(x_pc, dtype=float), np.asarray(y_pc, dtype=float))
    h_z = height_pc + flare_pc_per_kpc * np.maximum(r - flare_start_pc, 0.0) / 1000.0
    layer = 1.0 / np.cosh(np.asarray(z_pc, dtype=float) / h_z) ** 2
    return av_per_kpc_at_sun / 1000.0 * np.exp(-(r - solar_radius_pc) / scale_length_pc) * layer


def galactic_extinction_av(observer_pc, direction, distance_pc, *, steps=EXTINCTION_STEPS, **dust):
    """
    The V extinction `A_V` (magnitudes) of the dust disc along lines of sight from `observer_pc` (galactocentric `(x, y,
    z)` in parsecs, z up) towards stars at `distance_pc` in the directions `direction`: an array of shape `(n, 3)`
    (they are made unit vectors) with `n` distances, or one direction and one distance. The integral of
    `dust_av_per_pc` by the midpoint rule over `steps` steps of each line. Keyword arguments (`av_per_kpc_at_sun`,
    `scale_length_pc`, `height_pc`, ...) go to `dust_av_per_pc`, to scale the law for another galaxy.

    Returns:
        numpy.ndarray: `A_V` for each line (a float for one).
    """
    observer = np.asarray(observer_pc, dtype=float).reshape(3)
    u = np.atleast_2d(np.asarray(direction, dtype=float))
    d = np.atleast_1d(np.asarray(distance_pc, dtype=float))
    if np.any(d < 0) or not np.all(np.isfinite(d)):
        raise ValueError("a distance must be zero or more and finite")
    norm = np.linalg.norm(u, axis=1)
    if np.any(norm == 0):
        raise ValueError("a direction must not be the zero vector")
    u = u / norm[:, None]
    u = np.broadcast_to(u, (max(len(u), len(d)), 3))
    d = np.broadcast_to(d, (len(u),))
    fractions = (np.arange(steps) + 0.5) / steps
    total = np.empty(len(u))
    chunk = max(1, 200_000 // steps)
    for start in range(0, len(u), chunk):
        stop = start + chunk
        s = d[start:stop, None] * fractions[None, :]
        points = observer[None, None, :] + s[:, :, None] * u[start:stop, None, :]
        rate = dust_av_per_pc(points[..., 0], points[..., 1], points[..., 2], **dust)
        total[start:stop] = rate.sum(axis=1) * d[start:stop] / steps
    return float(total[0]) if np.ndim(distance_pc) == 0 and np.ndim(direction) == 1 else total
