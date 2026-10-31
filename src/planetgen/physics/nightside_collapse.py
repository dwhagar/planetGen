# planetgen/physics/nightside_collapse.py

"""
Whether the air of a tidally locked rocky planet freezes onto its dark side, and how much is left (GEN.1044, oddity W19 of
docs/design/systems/weird-worlds.md; docs/design/composition/nightside-collapse-model.md).

On a locked planet the day side is lit for ever and the night side never. If the air is thin it carries too little heat round,
the night side drops below the condensation temperature of the main gas, and the gas freezes onto the dark side until what is
left is the vapour pressure at the night side's temperature: the air is lost to a nightside cap and the day side is bare.
This module answers that, as a pure function of the pressure, the starlight, the gravity and the gas. It draws nothing and
changes nothing a world stores; the wiring is GEN.1051.

1. **Critical pressure of CO2** (`fit_critical_pressure_pa`). The empirical fit of Wordsworth 2015 (arXiv 1412.5575, Eq. 43 and
   Table 2): `log10(p_crit / 1 Pa)` is a cubic in the stellar flux, from a grid of 3D climate runs of a pure CO2 air over
   0.2 to 3.0 times Earth's flux and 0.01 to 10 bar, for one and ten Earth masses, surface albedo 0.2. About 0.14 bar at Earth's
   flux for an Earth-mass planet, falling to about 0.03 bar at three times it. The paper calls it a conservative upper limit.
   Joshi, Haberle and Reynolds 1997 found 0.03 bar with a grey gas (the 0.1 bar figure is Wordsworth's own).
2. **Night-side temperature** (`night_side_temperature_k`). The thin-air radiative result of the same paper (Eq. 10, 25): the
   night side receives the share `tau / (2 - tau)` of a quarter of the absorbed flux, `tau = kappa p / (g cos(alpha))`,
   kappa 5e-5 m2/kg, cos(alpha) 0.5, no warmer than the planet's own equilibrium temperature from tau of 1 up [stand-in NC3].
3. **What is left** (`residual_pressure_pa`). Under the critical pressure the loss feeds on itself (less gas, a colder night side,
   a lower vapour pressure), so the air does not settle part way; what stays is at most the vapour pressure at the night side's
   temperature (`physics/saturation.py`), a tiny share. Because the grey estimate of 2 puts the night side warmer than the 3D runs do,
   the CO2 case scales the night side's temperature by the one factor that makes the fit's critical pressure the point where
   the balance first loses gas, so the two agree at the threshold [stand-in NC4].
4. **Other gases** (`collapse` with another formula). The paper leaves N2 and CO to later work. For those the critical
   pressure is where the same balance first loses gas (no scaling): N2 freezes at 126 K and under, so it holds at any pressure on a
   world at Earth's flux and collapses on a much colder one [stand-in NC5].

Tests against measurement (never inputs, `tests/test_nightside_collapse.py`): the fit at Earth's flux is 0.14 bar, 0.1 bar of CO2
loses almost all of it, 1 bar of CO2 or air keeps it, a Mars-like 600 Pa of CO2 on a locked world loses almost all of it.

Limits. The fit is for a pure CO2 air over a flat rocky world with no ocean and no clouds, one star (AD Leo's spectrum), and
no obliquity or eccentricity; a non-condensing background gas such as N2 makes the real collapse pressure lower (the paper
calls the fit conservative, and its footnote 4 says neglecting such gases is conservative). Below the lowest flux the
paper studied (273.2 W/m2) it says the air collapses at every pressure it ran, but its own fit gives about 2.9 bar there; the fit
is held at its edges and the result marked extrapolated. The albedo enters only by the absorbed flux, scaling the flux the fit
takes by `(1 - A) / 0.8`. Heat of condensation, a night-side ocean and obliquity are not modelled.
"""

import dataclasses
import math

from planetgen.physics import constants
from planetgen.physics import saturation
from planetgen.util.checks import finite_domain

SOLAR_CONSTANT_W_M2 = 1366.0
"""float: Earth's stellar flux, W/m2, the unit of the paper's flux range (Wordsworth 2015, Table 1)."""

FIT_FLUX_RANGE_W_M2 = (0.2 * SOLAR_CONSTANT_W_M2, 3.0 * SOLAR_CONSTANT_W_M2)
"""tuple: The stellar flux range of the paper's runs, 273.2 to 4098 W/m2 (Table 1: 0.2 to 3.0 times 1366). The fit is held at either
end outside it and the result is marked extrapolated. [Stand-in NC1, tier T1, tweak L0.]"""

FIT_ALBEDO = 0.2
"""float: The surface albedo of every run behind the fit (Table 1). [Stand-in NC1.]"""

FIT_COEFFICIENTS = {
    1.0: (-0.84e-10, 0.73e-6, -0.0022, 6.01),
    10.0: (-0.81e-10, 0.63e-6, -0.0017, 6.34),
}
"""dict: Earth masses -> `(c1, c2, c3, c4)` of `log10(p_crit / 1 Pa) = c1 F^3 + c2 F^2 + c3 F + c4` with F in W/m2 (Wordsworth
2015, Table 2, the units of F being the figure's axis; the Earth-flux, one-Earth-mass value comes out at 0.14 bar as the abstract says).
Between one and ten Earth masses the coefficients are linear in the log of the mass and held outside. [Stand-in NC1, tier T1,
tweak L0; the interpolation is NC2, tier T7, tweak L1.]"""

KAPPA_M2_KG = 5.0e-5
"""float: The grey absorption coefficient of the paper's thin-air estimate, m2/kg (Eq. 10). [Stand-in NC3, tier T1, tweak L0.]"""

COS_ALPHA = 0.5
"""float: The mean slant of the paper's two-stream estimate (Eq. 10). [Stand-in NC3.]"""

PRESSURE_FLOOR_PA = 1.0e-12
"""float: The residual below which the gas is counted as all frozen (a vapour pressure of a ten-trillionth of a bar is none)."""

COLLAPSED_BELOW_RETAINED = 0.5
"""float: A world that keeps under this share of its gas is `collapsed`; the share itself is the number a caller should use.
[Stand-in NC6, tier T9, tweak L0.]"""


@dataclasses.dataclass(frozen=True)
class Collapse:
    """What a locked planet's air does at its night side."""
    retained_fraction: float
    """float: The share of the gas that stays in the air, 1 when the night side holds it, near 0 when it all froze."""
    residual_pressure_pa: float
    """float: The pressure of the gas that stays, Pa."""
    critical_pressure_pa: float
    """float: The pressure under which the gas starts to freeze, Pa."""
    collapsed: bool
    """bool: Whether `retained_fraction` is under `COLLAPSED_BELOW_RETAINED`."""
    source: str
    """str: `wordsworth-2015-fit` (CO2) or `radiative-balance` (any other gas)."""
    extrapolated: bool
    """bool: Whether the flux or the mass is outside what the fit was made over."""


def _fit_coefficients(mass_earth):
    """`(coefficients, extrapolated)` for a planet of `mass_earth`: linear in the log of the mass between the two tabulated, held outside."""
    low, high = sorted(FIT_COEFFICIENTS)
    if mass_earth <= low:
        return FIT_COEFFICIENTS[low], mass_earth < low
    if mass_earth >= high:
        return FIT_COEFFICIENTS[high], mass_earth > high
    t = math.log(mass_earth / low) / math.log(high / low)
    return tuple(a + (b - a) * t for a, b in zip(FIT_COEFFICIENTS[low], FIT_COEFFICIENTS[high])), False


@finite_domain()
def fit_critical_pressure_pa(flux_w_m2, mass_earth=1.0, albedo=FIT_ALBEDO):
    """
    The CO2 collapse pressure of Wordsworth 2015's fit, Pa, and whether it had to be held at the edge of what the paper ran.

    Args:
        flux_w_m2 (float): The stellar flux at the planet, W/m2, above zero.
        mass_earth (float): The planet's mass in Earth masses, above zero (the fit has one and ten).
        albedo (float): The surface albedo, from 0 to under 1; the flux the fit takes is `flux (1 - albedo) / (1 - 0.2)`.

    Returns:
        tuple: `(pressure_pa, extrapolated)`.
    """
    if flux_w_m2 <= 0.0 or mass_earth <= 0.0 or not 0.0 <= albedo < 1.0:
        raise ValueError(f"flux {flux_w_m2}, mass {mass_earth} and albedo {albedo} must be above zero, above zero and in [0, 1)")
    effective = flux_w_m2 * (1.0 - albedo) / (1.0 - FIT_ALBEDO)
    low, high = FIT_FLUX_RANGE_W_M2
    held = min(max(effective, low), high)
    coefficients, mass_held = _fit_coefficients(mass_earth)
    c1, c2, c3, c4 = coefficients
    return 10.0 ** (((c1 * held + c2) * held + c3) * held + c4), (held != effective) or mass_held


@finite_domain()
def night_side_temperature_k(pressure_pa, gravity_ms2, flux_w_m2, albedo=FIT_ALBEDO):
    """
    The night side's temperature under a thin air, kelvin (Wordsworth 2015, Eq. 10 and 25): `tau = kappa p / (g cos(alpha))`,
    the night side radiates the share `tau / (2 - tau)` of `(1 - A) F / 4`, held at the planet's equilibrium temperature from
    `tau` of 1 up.
    """
    if gravity_ms2 <= 0.0 or flux_w_m2 <= 0.0:
        raise ValueError("gravity and flux must be above zero")
    absorbed = (1.0 - albedo) * flux_w_m2 / 4.0
    if pressure_pa <= 0.0:
        return 0.0
    tau = KAPPA_M2_KG * pressure_pa / (gravity_ms2 * COS_ALPHA)
    share = 1.0 if tau >= 1.0 else tau / (2.0 - tau)
    return (share * absorbed / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25


def _vapour_pa(formula, temperature_k):
    """The vapour pressure of `formula` at `temperature_k`, Pa (infinite above its critical temperature)."""
    return saturation.pressure_pa(formula, max(temperature_k, 1.0e-3))


def _loses_gas(formula, pressure_pa, gravity_ms2, flux_w_m2, albedo, scale):
    """Whether a pure `formula` air at `pressure_pa` is over the vapour pressure at its (scaled) night side."""
    night = scale * night_side_temperature_k(pressure_pa, gravity_ms2, flux_w_m2, albedo)
    return pressure_pa > _vapour_pa(formula, night)


def _balance_critical_pa(formula, gravity_ms2, flux_w_m2, albedo, scale, top_pa=1.0e8):
    """
    The pressure under which the balance loses gas, Pa: the lower edge of the window where the night side holds the gas
    (vapour pressure over the pressure), found by stepping up a decade at a time from the floor and bisecting; `top_pa` where no
    pressure under it holds, 0 where the floor itself does.
    """
    p = PRESSURE_FLOOR_PA
    if not _loses_gas(formula, p, gravity_ms2, flux_w_m2, albedo, scale):
        return 0.0
    while p < top_pa:
        if not _loses_gas(formula, p * 10.0, gravity_ms2, flux_w_m2, albedo, scale):
            low, high = p, p * 10.0
            for _ in range(80):
                mid = math.sqrt(low * high)
                low, high = (mid, high) if _loses_gas(formula, mid, gravity_ms2, flux_w_m2, albedo, scale) else (low, mid)
            return high
        p *= 10.0
    return top_pa


def _calibration(gravity_ms2, flux_w_m2, albedo, critical_pa):
    """The factor on the grey night-side temperature that makes the balance first lose CO2 at the fit's `critical_pa`. [Stand-in NC4, tier T9, tweak L2.]"""
    cold = saturation.condensation_temperature_k("CO2", critical_pa)
    night = night_side_temperature_k(critical_pa, gravity_ms2, flux_w_m2, albedo)
    if cold is None or night <= 0.0:
        return 1.0
    return cold / night


@finite_domain()
def residual_pressure_pa(formula, pressure_pa, gravity_ms2, flux_w_m2, albedo=FIT_ALBEDO, scale=1.0):
    """
    The pressure the gas keeps once it has frozen onto the night side, Pa: its vapour pressure at the night side's temperature
    under the air it started with, at most that air's pressure. The loss cools the night side further, so this is an upper bound
    and a collapsed world's real residue is smaller still; it is a tiny share of the air for any gas well under its critical pressure.
    """
    if pressure_pa <= 0.0:
        return 0.0
    night = scale * night_side_temperature_k(pressure_pa, gravity_ms2, flux_w_m2, albedo)
    return min(pressure_pa, _vapour_pa(formula, night))


@finite_domain()
def collapse(formula, pressure_pa, gravity_ms2, flux_w_m2, mass_earth=1.0, albedo=FIT_ALBEDO):
    """
    What a pure `formula` air of `pressure_pa` does at the night side of a tidally locked planet.

    Args:
        formula (str): The main gas, a species of `physics/saturation.py` (`CO2`, `N2`, ...); the air is taken to be that gas alone.
        pressure_pa (float): The air's pressure before any of it froze, Pa.
        gravity_ms2 (float): The surface gravity, m/s2.
        flux_w_m2 (float): The stellar flux at the planet, W/m2.
        mass_earth (float): The planet's mass in Earth masses (selects the CO2 fit's coefficients).
        albedo (float): The surface albedo.

    Returns:
        Collapse: Of a planet that does not lock, the caller should not ask.
    """
    if pressure_pa < 0.0:
        raise ValueError("a pressure cannot be negative")
    if pressure_pa == 0.0:
        return Collapse(1.0, 0.0, 0.0, False, "radiative-balance", False)
    if formula == "CO2":
        critical, extrapolated = fit_critical_pressure_pa(flux_w_m2, mass_earth, albedo)
        scale, source = _calibration(gravity_ms2, flux_w_m2, albedo, critical), "wordsworth-2015-fit"
    else:
        scale, source, extrapolated = 1.0, "radiative-balance", False
        critical = _balance_critical_pa(formula, gravity_ms2, flux_w_m2, albedo, scale)
    if pressure_pa >= critical:
        return Collapse(1.0, pressure_pa, critical, False, source, extrapolated)
    left = min(pressure_pa, residual_pressure_pa(formula, pressure_pa, gravity_ms2, flux_w_m2, albedo, scale))
    retained = left / pressure_pa
    return Collapse(retained, left, critical, retained < COLLAPSED_BELOW_RETAINED, source, extrapolated)
