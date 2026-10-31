# planetgen/physics/radiation_belts.py

"""
Trapped radiation belts of a magnetised body and the dose at an orbit inside them (GEN.749)
=============================================================================================

A body with a dipole field holds the charged particles of its surroundings in belts around it. This module gives the belts of a body from its
dipole moment, its radius, the standoff of its magnetopause (which the star's wind sets, `magnetism.wind_pressure_pa`) and its spin, and the dose
rate a spacecraft or a moon takes at a distance `L` (equatorial distance in body radii) behind a thin shield.

**Geometry.** In a dipole the field at `L` on the equator is `B_eq / L^3`, `B_eq = mu0 M / (4 pi R^3)`; a particle bounces along the field line of
shell `L` and drifts around the body, and it stays while the field is steady. The belts have

* an **outer edge** `L_max`, the smaller of the wind's limit and the spin's. The wind's is `WIND_EDGE_FRACTION` (0.75) of the magnetopause standoff
  `R_mp = 9.75 R_earth (M / M_earth)^(1/3) (p / p_earth)^(-1/6)` (Chapman and Ferraro: the field's pressure against the wind's), which puts Earth's at
  7.3 (the outer belt ends at about 40 000 km). The spin's is `ROTATION_EDGE_FACTOR` (9) times the synchronous distance `L_sync = (G m P^2 /
  (4 pi^2 R^3))^(1/3)`, beyond which the plasma no longer turns with the body and its field is stretched into a disc (`L` of 20 for Jupiter, about
  16 for Saturn). A slow rotator is bounded by the wind; a fast one by its spin.
* an **inner edge** `L_in`: the surface (`L` of 1) or, for a body with an atmosphere, the height at which the air stops a bouncing particle,
  `ATMOSPHERE_EDGE_SCALE_HEIGHTS` (45) scale heights up (Earth's belts begin at about 400 km, `L` of 1.06).
* an **outer belt** of electrons of MeV energies that peaks at `L_pk`, `WIND_PEAK_FRACTION` (0.62) of `L_max` for a wind-limited body (Earth's peak
  is at 4 to 5) and `ROTATION_PEAK_FRACTION` (0.25) of it for a spin-limited one (Jupiter's is inside Io's orbit, 5.9, of 20), and an **inner belt**
  of protons of tens of MeV, which cosmic-ray albedo neutrons make where they decay in the field, at `INNER_PEAK_FRACTION` (0.35) of that (Earth's
  at 1.5 to 2), between the two a slot where little is trapped for long.

**Intensity.** The dose rate behind `SHIELD_G_CM2` of aluminium (100 mil, 0.69 g/cm2, the shield a spacecraft part is rated against) at the peak of each belt is
Earth's (`EARTH_OUTER_PEAK_GY_YR`, `EARTH_INNER_PEAK_GY_YR`) times `(M_dip / M_earth)^INTENSITY_EXPONENT`: a stronger field holds particles of higher energy
and more of them, because a particle diffusing inward gains energy as the field grows (the first adiabatic invariant), and the belts of Jupiter are some thousands of
times Earth's at the peak. The profile of each belt is a rise `((L - L_in) / (L_pk - L_in))^p` to the peak (`INNER_RISE_EXPONENT` 2 for the protons, which reach down to the air, `OUTER_RISE_EXPONENT` 4 for the electrons, which leave a slot between the belts), a fall `(L_pk / L)^FALL_EXPONENT` beyond it and
a cut `1 - (L / L_max)^EDGE_EXPONENT` to nothing at the edge (`profile`).

**What a moon in them takes (GEN.983).** A planet's own belts begin above its air, so its ground takes nothing from them. A moon that orbits inside the belts of its
parent (`moon_dose_msv_yr`) takes the dose at its orbit's `L`: the sum of the belts' rates, each in Sv by the weight of its particles (`BELT_WEIGHT`: 1 for electrons, 2
for protons, ICRP 103's radiation weights), seen from the side of the sky the moon does not block (`SURFACE_VIEW_FRACTION`, 0.5), and cut by its own air by `exp(-column /
STOP_COLUMN_G_CM2)`: the particles are electrons of MeV and protons of tens of MeV, which a few g/cm2 stop. The moon's own field and its minimagnetosphere are not
taken (the ones a moon draws are weak beside the parent's belts). A moon deep in the belts of a Jupiter-strength field takes thousands of Sv/yr, a moon outside the belts nothing.

**What it does not do.** It is a model of the main belts of a dipole: no tilt or offset of the field (the South Atlantic Anomaly), no multipole body (the dipole
share only), no moons or rings that sweep the belts clean (Saturn's), no plasma source of a moon (Io's torus), no solar-cycle or storm variation, no heavy ions,
and the shield is aluminium of one thickness. The intensities are Earth's scaled, not derived; the stand-ins and their tiers are in the gap-filling plan (section 6.15).
"""

import math
from dataclasses import dataclass

from planetgen.physics import air_column, constants, magnetism

VACUUM_PERMEABILITY = magnetism.VACUUM_PERMEABILITY
EARTH_MOMENT_A_M2 = magnetism.EARTH_MOMENT_A_M2

WIND_EDGE_FRACTION = 0.75
ROTATION_EDGE_FACTOR = 9.0
WIND_PEAK_FRACTION = 0.62
ROTATION_PEAK_FRACTION = 0.25
INNER_PEAK_FRACTION = 0.35
ATMOSPHERE_EDGE_SCALE_HEIGHTS = 45.0
INNER_RISE_EXPONENT = 2.0
OUTER_RISE_EXPONENT = 4.0
FALL_EXPONENT = 6.0
EDGE_EXPONENT = 6.0
INTENSITY_EXPONENT = 0.8
EARTH_OUTER_PEAK_GY_YR = 3000.0
EARTH_INNER_PEAK_GY_YR = 3000.0
SHIELD_G_CM2 = 0.69
"""float: 100 mil of aluminium, g/cm2 (2.54 mm at 2.7 g/cm3)."""
BELT_WEIGHT = {"electrons": 1.0, "protons": 2.0}
"""dict: Sv per Gy of each belt's particles: the radiation weighting factors of ICRP 103 (photons and electrons 1, protons 2)."""
SURFACE_VIEW_FRACTION = 0.5
"""float: The share of the belt's dose a ground takes: the body blocks half the sky."""
STOP_COLUMN_G_CM2 = 2.0
"""float: The air column, g/cm2, that cuts the belts' dose by a factor of e (MeV electrons go about 1 g/cm2 and tens of MeV protons 1 to 3)."""
EARTH_STANDOFF_RP = magnetism.EARTH_STANDOFF_RP
EARTH_WIND_PA = magnetism.EARTH_WIND_PA
EARTH_RADIUS_M = constants.EARTH_RADIUS_KM * 1000.0


@dataclass(frozen=True)
class Belt:
    """One belt: its species, where it starts, peaks and ends (`L`, body radii on the equator) and its peak dose rate, Gy/yr behind `SHIELD_G_CM2` of aluminium."""

    name: str
    species: str
    inner_l: float
    peak_l: float
    outer_l: float
    peak_gy_yr: float
    rise_exponent: float

    def dose_gy_yr(self, l_shell):
        """The belt's dose rate at `l_shell`, Gy/yr behind `SHIELD_G_CM2` of aluminium (0 outside it)."""
        return self.peak_gy_yr * profile(l_shell, self.inner_l, self.peak_l, self.outer_l, self.rise_exponent)


def profile(l_shell, inner_l, peak_l, outer_l, rise_exponent=OUTER_RISE_EXPONENT):
    """The belt's radial shape, 1 at its peak: a rise from the inner edge, a fall beyond the peak and a cut to nothing at the outer edge."""
    if l_shell <= inner_l or l_shell >= outer_l or peak_l <= inner_l:
        return 0.0
    if l_shell <= peak_l:
        return ((l_shell - inner_l) / (peak_l - inner_l)) ** rise_exponent
    return (peak_l / l_shell) ** FALL_EXPONENT * (1.0 - (l_shell / outer_l) ** EDGE_EXPONENT) / (1.0 - (peak_l / outer_l) ** EDGE_EXPONENT)


def equatorial_field_t(moment_a_m2, radius_m):
    """The dipole's field at the equator of the surface, tesla."""
    return VACUUM_PERMEABILITY * moment_a_m2 / (4.0 * math.pi * radius_m ** 3)


def standoff_m(moment_a_m2, wind_pressure_pa):
    """The magnetopause standoff, m: Earth's 9.75 radii scaled by the cube root of the moment and the sixth root of the wind's pressure."""
    return (EARTH_STANDOFF_RP * EARTH_RADIUS_M * (moment_a_m2 / EARTH_MOMENT_A_M2) ** (1.0 / 3.0)
            * (wind_pressure_pa / EARTH_WIND_PA) ** (-1.0 / 6.0))


def synchronous_l(mass_kg, radius_m, rotation_s):
    """The distance, body radii, at which a circular orbit turns once with the body: `(G m P^2 / (4 pi^2 R^3))^(1/3)`."""
    return (constants.G * mass_kg * rotation_s ** 2 / (4.0 * math.pi ** 2 * radius_m ** 3)) ** (1.0 / 3.0)


def outer_edge(standoff_radii, sync_radii):
    """`(L_max, spin_limited)`: the wind's limit or the spin's, whichever is nearer."""
    wind = WIND_EDGE_FRACTION * standoff_radii
    spin = ROTATION_EDGE_FACTOR * sync_radii
    return (spin, True) if spin < wind else (wind, False)


def belts_for(moment_a_m2, radius_m, mass_kg, rotation_s, wind_pressure_pa, scale_height_m=None):
    """
    The belts of a dipole of `moment_a_m2` on a body of `radius_m` and `mass_kg` turning in `rotation_s`, in a wind of `wind_pressure_pa` (all SI), as a
    tuple of `Belt`s (the inner proton belt, the outer electron belt), or `()` with no field or a field the wind squeezes onto the body. `scale_height_m`
    is the air's, for the inner edge; `None` for a body without air.
    """
    if not moment_a_m2 or moment_a_m2 <= 0 or radius_m <= 0 or wind_pressure_pa <= 0:
        return ()
    standoff_radii = standoff_m(moment_a_m2, wind_pressure_pa) / radius_m
    sync = synchronous_l(mass_kg, radius_m, rotation_s) if rotation_s and rotation_s > 0 else float("inf")
    outer, spin_limited = outer_edge(standoff_radii, sync)
    inner = 1.0 + (ATMOSPHERE_EDGE_SCALE_HEIGHTS * scale_height_m / radius_m if scale_height_m else 0.0)
    if outer <= inner * 1.2:
        return ()
    peak = (ROTATION_PEAK_FRACTION if spin_limited else WIND_PEAK_FRACTION) * outer
    peak = max(peak, inner * 1.5)
    scale = (moment_a_m2 / EARTH_MOMENT_A_M2) ** INTENSITY_EXPONENT
    proton_peak = max(INNER_PEAK_FRACTION * peak, inner * 1.2)
    return (
        Belt("inner", "protons", inner, proton_peak, peak, EARTH_INNER_PEAK_GY_YR * scale, INNER_RISE_EXPONENT),
        Belt("outer", "electrons", inner, peak, outer, EARTH_OUTER_PEAK_GY_YR * scale, OUTER_RISE_EXPONENT),
    )


def air_scale_height_m(body):
    """The air's scale height, `R T / (M g)` at the body's surface temperature, m, or `None` for a body with no air (or no temperature)."""
    column = air_column.column_of(body)
    temperature = getattr(body, "surface_temperature", None)
    if column is None or not temperature:
        return None
    return constants.R * temperature / (column.mean_molar_mass_kg_mol * column.gravity_m_s2)


def belts_of_body(body, star_distance_au=None):
    """
    The belts of a planet or moon from what it stores: its dipole moment, radius (km) and mass (kg),
    rotation period, air scale height, and the star's wind at its distance (`magnetism.wind_pressure_pa`). `()` for a body with no field, no star or
    no distance.
    """
    moment = getattr(body, "magnetic_moment_a_m2", None)
    distance = star_distance_au or getattr(body, "_star_distance_au", None)
    if not moment or not distance or getattr(body, "star", None) is None:
        return ()
    rotation_h = getattr(body, "rotation_period_hours", None)
    return belts_for(moment, body.radius * 1000.0, body.mass,
                     rotation_h * 3600.0 if rotation_h else None, magnetism.wind_pressure_pa(body.star, distance), air_scale_height_m(body))


def dose_gy_yr(belts, l_shell):
    """The dose rate at `l_shell`, Gy/yr behind `SHIELD_G_CM2` of aluminium: the sum of the belts' (0 for none)."""
    return sum(belt.dose_gy_yr(l_shell) for belt in belts)


def dose_at_orbit_gy_yr(belts, orbit_radius_m, body_radius_m):
    """The dose rate of a circular equatorial orbit of `orbit_radius_m` around a body of `body_radius_m`, Gy/yr behind the standard shield."""
    return dose_gy_yr(belts, orbit_radius_m / body_radius_m)


def air_cut(column_g_cm2):
    """The share of the belts' dose that gets through an air of `column_g_cm2` g/cm2 (1 for an airless body)."""
    return math.exp(-max(column_g_cm2, 0.0) / STOP_COLUMN_G_CM2)


def moon_exposure_msv_yr(parent, orbit_au):
    """
    The dose, mSv/yr, an airless moon at `orbit_au` from `parent` takes from the parent's belts: each belt's dose behind the standard shield at the moon's `L` in
    Sv by its weight, seen from half the sky. `0.0` when the parent has no belts or the moon is outside them.
    """
    belts = belts_of_body(parent)
    if not belts or not orbit_au or orbit_au <= 0 or not parent.radius:
        return 0.0
    l_shell = orbit_au * constants.AU_M / (parent.radius * 1000.0)
    gray_per_year = sum(BELT_WEIGHT[belt.species] * belt.dose_gy_yr(l_shell) for belt in belts)
    return gray_per_year * 1000.0 * SURFACE_VIEW_FRACTION
