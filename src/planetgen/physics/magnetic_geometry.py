"""
Where a magnetic dipole points, and the surface field of a star (GEN.1350)
==========================================================================

The generator stored the size of a body's magnetic field (`magnetic_moment_a_m2`, `magnetic_field_gauss`) but not its direction. This module draws the direction
from the physics that decides it, and a surface dipole field for the stars the sources let it. Nothing here moves the draws of the run: every number comes from a
stream keyed on the body alone (`draw.Stream`, the way `oddities` and `chemically_peculiar` do it), so a body made again is the same and no other draw shifts.

**Planets, moons and rogue planets (`of_body`).** `dipole_tilt_deg` is the angle between the dipole's axis and the spin axis, folded into 0 to 90 degrees (the axis, not
the direction of the field), and `dipole_offset_class` is `centred` or `offset` (the dipole's centre more than `OFFSET_CLASS_RP` of a radius from the body's centre).
Dynamo simulations of a thick, fast-rotating shell give a dipole close to the spin axis ([R0784]: Earth 11 degrees, Jupiter 9.6, Ganymede 4, Saturn under 0.1, Mercury under 3
degrees); a thin shell round a stable interior, the geometry of the ice giants, gives a dipole at a large angle that is off centre ([R4196]: Uranus 59 degrees and 0.3 radii,
Neptune 47 and 0.55). The model:

* a **dipole-dominated dynamo** (the dipole holds at least half the field, `magnetism.generate_field`'s share): the tilt is Rayleigh with scale `AXIAL_TILT_SIGMA_DEG`;
  a centred dipole, except the small iron-rich body of the Mercury stand-in, which is `offset`;
* a **multipolar dynamo** (the dipole under half): the axis is isotropic, `cos(tilt)` uniform, and the dipole is `offset`;
* an **ice giant** (the conducting layer is the water of the mantle, a thin shell over a stable core): isotropic and `offset`;
* a **gas giant or brown dwarf** (the metallic hydrogen layer): Rayleigh, centred.

A body with no field has neither. A rogue planet and a moon take the same rules as a planet.

**Neutron stars (`neutron_star_inclination_deg`).** The angle between the magnetic and the spin axis. A neutron star is born with an isotropic axis; the dipole's torque
brings the two axes together (`tan(alpha)` falls as `exp(-t / tau)`, `ALIGNMENT_TIMESCALE_YR` after [R5143]); a recycled millisecond pulsar has been spun up by accretion and is
isotropic again. Folded into 0 to 90 degrees.

**Stars (`star_field_gauss`).** The polar field of the large-scale dipole, gauss. A cool dwarf or subdwarf with a convective envelope (below the Kraft break) follows the
empirical trend of the large-scale field with the Rossby number, `B ~ Ro^-1.40` ([R4544, R5144] for the fully convective M dwarfs), with `Ro = P_rot / tau_conv` and the
turnover time of [R4780]; an Ap or Bp dwarf has the polar field `chemically_peculiar` draws for its kind. Any other star has none stored.

Stand-ins are DP1 to DP7 of `docs/design/research/gap-filling-plan.md`.
"""

import math

from planetgen.physics import chemically_peculiar, spin
from planetgen.util import draw

DIPOLE_DOMINATED_SHARE = 0.5
"""float: A dynamo whose dipole holds at least this share of its field is dipole-dominated (`magnetism.generate_field` calls a lower one multipolar). [Stand-in DP1, tier T5, tweak level L0.]"""

AXIAL_TILT_SIGMA_DEG = 6.0
"""float: The Rayleigh scale, degrees, of the tilt of a dipole-dominated dynamo's axis from the spin axis: a median of 7 degrees, 18 percent of bodies over 11 (Earth), Jupiter's
9.6, Ganymede's 4, Saturn's 0 and Mercury's under 3 inside the draw. [Stand-in DP2, tier T4, tweak level L2.]"""

OFFSET_CLASS_RP = 0.15
"""float: The dipole's offset, in radii, over which the class is `offset`: Earth 0.07 and Jupiter 0.1 lie under it, Mercury 0.19, Uranus 0.3 and Neptune 0.55 over it. [Stand-in DP3, tier T4, tweak level L0.]"""

ALIGNMENT_TIMESCALE_YR = 7.0e6
"""float: The e-folding time, years, in which the dipole of a neutron star aligns with its spin (`tan(alpha)` falls as `exp(-t / tau)`). [Stand-in DP4, tier T3, tweak level L1.]"""

B_SUN_POLAR_GAUSS = 1.5
"""float: The polar dipole field, gauss, of the Sun at its activity-cycle mean (1 to 3 gauss over the cycle). [Stand-in DP5, tier T2, tweak level L1.]"""

ROSSBY_SLOPE = -1.40
"""float: The slope of the large-scale field against the Rossby number ([R4544]: -1.40 +/- 0.14). [Stand-in DP6, tier T2, tweak level L0.]"""

FIELD_SCATTER_DEX = 0.35
"""float: The scatter, dex, of a cool dwarf's dipole field about the trend. [Stand-in DP7, tier T4, tweak level L1.]"""

MAX_FIELD_GAUSS = 5000.0
"""float: The cap, gauss, of a cool dwarf's dipole (the strongest large-scale fields of the M dwarfs are about 1 to 2 kG; [R5144]). [Stand-in DP7, tier T4, tweak level L1.]"""

SUN_MASS_SOL = 1.0
SUN_PERIOD_DAYS = 25.4
"""float: The Sun's rotation period, days, which gives the Rossby number the trend is anchored at."""

MIN_STAR_MASS_SOL = 0.08
MAX_STAR_MASS_SOL = 1.3
"""float: The masses, solar, over which the turnover time of [R4780] is fitted."""

DWARF_CLASSES = ("V", "VI")
"""tuple: The luminosity classes the trend covers: a main-sequence dwarf or a subdwarf (an old, metal-poor dwarf is as convective as a young one)."""

OFFSET_KINDS = ("ice-giant", "multipolar", "mercury")
"""tuple: The kinds of dynamo whose dipole is `offset`."""


def _fold(angle_deg):
    """An axis angle: 0 to 180 degrees folded into 0 to 90."""
    return 180.0 - angle_deg if angle_deg > 90.0 else angle_deg


def _isotropic_deg(stream):
    """The angle, degrees, of an axis drawn isotropically from a fixed one, folded into 0 to 90."""
    return _fold(math.degrees(math.acos(1.0 - 2.0 * stream.random())))


def _rayleigh_deg(stream, sigma_deg):
    return _fold(min(sigma_deg * math.sqrt(-2.0 * math.log(1.0 - stream.random())), 180.0))


def of_body(kind, key):
    """
    `(dipole_tilt_deg, dipole_offset_class)` of a body whose dynamo is of `kind`: `axial` (a dipole-dominated rocky core or a gas giant), `multipolar`, `ice-giant` or
    `mercury` (the iron-rich small body whose field `magnetism` keeps at one value). `key` names the body's stream.
    """
    stream = draw.Stream(f"dipole-tilt:{key}")
    tilt = _isotropic_deg(stream) if kind in ("multipolar", "ice-giant") else _rayleigh_deg(stream, AXIAL_TILT_SIGMA_DEG)
    return tilt, ("offset" if kind in OFFSET_KINDS else "centred")


def neutron_star_inclination_deg(pulsar_type, age_yr, key):
    """The angle, degrees, between a neutron star's magnetic and spin axes: isotropic at birth, aligned by the dipole's torque over `age_yr` (not for a recycled
    millisecond pulsar), folded into 0 to 90."""
    alpha = _isotropic_deg(draw.Stream(f"dipole-neutron-star:{key}"))
    if pulsar_type == "millisecond":
        return alpha
    return math.degrees(math.atan(math.tan(math.radians(min(alpha, 89.999))) * math.exp(-max(age_yr, 0.0) / ALIGNMENT_TIMESCALE_YR)))


def convective_turnover_days(mass_sol):
    """The convective turnover time, days, of a dwarf of this mass: `log tau = 1.16 - 1.49 log M - 0.54 (log M)^2` [R4780]."""
    log_mass = math.log10(mass_sol)
    return 10.0 ** (1.16 - 1.49 * log_mass - 0.54 * log_mass * log_mass)


def rossby_number(mass_sol, period_hours):
    """`Ro = P_rot / tau_conv` of a dwarf."""
    return period_hours / 24.0 / convective_turnover_days(mass_sol)


def covered_cool_dwarf(yerkes_class, mass_sol, temperature_k, period_hours):
    """Whether the Rossby trend covers a star: a main-sequence dwarf or subdwarf below the Kraft break, of a mass the turnover time is fitted over, with a rotation period."""
    return bool(yerkes_class in DWARF_CLASSES and period_hours and temperature_k is not None and MIN_STAR_MASS_SOL <= mass_sol <= MAX_STAR_MASS_SOL
                and spin.kraft_weight(temperature_k) == 0.0)


def star_field_gauss(yerkes_class, mass_sol, mass_kg, radius_km, temperature_k, age_gy, period_hours):
    """
    The polar field of a star's dipole, gauss, or `None` for a star the sources do not cover: a cool main-sequence dwarf with a convective envelope that has a rotation
    period (the trend with the Rossby number), an Ap or Bp dwarf (the polar field of its kind), no other.
    """
    if yerkes_class not in DWARF_CLASSES or None in (mass_kg, radius_km, temperature_k, age_gy):
        return None
    peculiar = chemically_peculiar.of_star(yerkes_class, mass_kg, radius_km, temperature_k, age_gy)
    if peculiar is not None:
        return peculiar.field_gauss          # None for an Am star: no detected field
    if not covered_cool_dwarf(yerkes_class, mass_sol, temperature_k, period_hours):
        return None
    anchor = rossby_number(SUN_MASS_SOL, SUN_PERIOD_DAYS * 24.0)
    trend = B_SUN_POLAR_GAUSS * (rossby_number(mass_sol, period_hours) / anchor) ** ROSSBY_SLOPE
    scatter = 10.0 ** draw.Stream(f"dipole-star:{mass_kg!r}:{temperature_k!r}:{age_gy!r}").gauss(0.0, FIELD_SCATTER_DEX)
    return min(trend * scatter, MAX_FIELD_GAUSS)
