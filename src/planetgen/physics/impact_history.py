# planetgen/physics/impact_history.py

"""
Late accretion, impact delivery and impact erosion of an air (GEN.739; docs/design/composition/impact-erosion-and-late-delivery.md)
===================================================================================================================================

After the giant impacts a rocky planet goes on accreting leftover planetesimals (the late veneer and the late accretion). They bring
volatiles and they erode the air that is already there. The balance for each volatile species `i` is the one of Sakuraba, Kurokawa and
Genda 2019 (arXiv:1805.07094, Eq. 1 and 2), in the cumulative impactor mass `S`:

    d m_i / dS = (1 - zeta) x_i - eta m_i / m_air

where `x_i` is the mass fraction of the species in the impactor, `zeta` the share of the impactor's vapour that escapes (Shuvalov 2009),
and `eta` the air mass an impactor removes per unit of its own mass.

* *How much accretes.* `late_accreted_fraction`: 0.7 weight percent of an Earth mass times `(M / M_earth)^(-1/3)` (the geometric cross
  section over the mass), log-normal with the width of the N-body scatter. Walker 2009 (0.5 to 0.8 weight percent for the Earth from the
  highly siderophile elements) and Table 1 of Brasser et al. 2016 (arXiv:1609.01785: Mercury 2.7, Venus 0.8, Earth 0.7, Mars 1.4) give the
  scale and the law.
* *What it brings.* `VENEER_VOLATILES`, the carbonaceous reference of Sakuraba et al. times a drawn scale (0.01 to 1; 0.1 is the
  enstatite-chondrite-like veneer that studies find for the Earth).
* *How hard it hits.* `impact_speed_ms`: `sqrt(v_esc^2 + v_inf^2)` with `v_inf = 0.6 v_orbit` (Table 1 of Brasser et al. gives 0.51 to 0.71).
* *What it removes.* `erosion_per_impactor_mass` (Svetsov's Eq. 5 and 6 as written by Sakuraba et al., averaged over a size distribution
  whose mass is flat in the diameter between 10^-1.5 and 10^3 km, Bottke et al. 2005, as they use it).

The gas that a species can hold is capped by the surface (`caps_pa`): an ocean holds the water above its vapour pressure, the
carbonate cycle the CO2 above its balance, an ice cap the rest; the excess is kept, not eroded (Sakuraba et al. section 2.3). The caller
gives the caps; the module knows no climate.

Nothing here knows the generator: it takes numbers and returns masses. The supply law (`surface_pressure.impact_inventory`) turns them into the
element budget the mixture draw reads.
"""

import math

from planetgen.physics import constants as physical_constants
from planetgen.physics.mixture_base import bracket
from planetgen.util import detmath, draw

G = physical_constants.G
"""float: Gravitational constant, SI."""

K_BOLTZMANN = 1.380649e-23
"""float: Boltzmann constant, J/K."""

AMU_KG = 1.66053907e-27
"""float: Atomic mass unit, kg."""

MOLAR_MASS = {"H2O": 18.015, "CO2": 44.010, "N2": 28.013}
"""dict: Species -> molar mass, g/mol."""

SPECIES = tuple(MOLAR_MASS)

LATE_ACCRETED_PERCENT_EARTH = 0.7
"""float: Mass accreted after the last giant impact, weight percent of an Earth mass planet at 1 AU of a solar mass star (Brasser et al.
2016 Table 1: 0.7 plus or minus 0.5 from the N-body runs; Walker 2009: 0.5 to 0.8 from the highly siderophile elements). Tier 1 for the
Solar System, a stand-in for any other system (the leftover planetesimal mass of 1e-3 Earth masses is the Solar System's)."""

LATE_ACCRETED_MASS_EXPONENT = -1.0 / 3.0
"""float: The accreted fraction goes as the cross section over the mass, `R^2 / M` with `R` going as `M^(1/3)`. Checked against the four means of
Brasser et al. 2016 Table 1 (Mercury, Venus, Earth, Mars: 1.9, 0.75, 0.7 and 1.5 predicted against 2.7, 0.8, 0.7 and 1.4). Tier 5 (derived)."""

LATE_ACCRETED_LN_SIGMA = 0.65
"""float: The scatter of the accreted fraction, the natural-log standard deviation (the Earth's 0.7 plus or minus 0.5 of Brasser et al. is
0.65 in the log). Tier 2, tweak level 2 (tuned to a published scatter)."""

LATE_ACCRETED_PERCENT_MAX = 3.0
"""float: The cap of a draw, weight percent (about the largest total of Brasser et al.'s runs for a terrestrial planet)."""

VENEER_VOLATILES = {"H2O": 0.10, "CO2": 0.02, "N2": 0.002}
"""dict: Mass fraction of a carbonaceous-chondrite-like impactor: 10 percent water, 2 percent CO2 and
0.2 percent N2 (Grady and Wright 2003), as used by Sakuraba et al. 2019 section 2.4 (read second-hand: tier 2, Confident Unverified)."""

VENEER_VOLATILE_SCALE = (0.01, 1.0, "log")
"""tuple: The scale of `VENEER_VOLATILES`, drawn: 1 is carbonaceous, 0.1 enstatite-chondrite-like (the nominal value of Sakuraba et al.), 0.01 their lower value. Stand-in, tier 6, tweak level 1."""

V_INF_OVER_V_ORBIT = 0.6
"""float: The approach speed of the leftover planetesimals over the orbital speed. Brasser et al. 2016 Table 1 give impact speeds of 34.2, 24.3,
21.1 and 13.4 km/s for Mercury, Venus, Earth and Mars, which with the escape speeds are 0.71, 0.63, 0.60 and 0.51 of the orbital speed. Tier 5."""

DIAMETER_RANGE_KM = (10.0 ** -1.5, 1000.0)
"""tuple: The impactor diameters, km: the main-belt size distribution of Bottke et al. 2005 (cumulative number as D^-3, so the mass per unit
diameter is flat) between the limits Sakuraba et al. 2019 use. Second-hand, tier 2."""

IMPACTOR_DENSITY = 3320.0
"""float: kg/m^3, Shuvalov 2009 (as used by Sakuraba et al.)."""

TARGET_DENSITY = 2630.0
"""float: kg/m^3, the planet's crust material in Shuvalov 2009."""

SVETSOV_K = 5
"""int: The expansion exponent of the vapour plume for an adiabatic index of 13/11 (Svetsov 2007)."""

SVETSOV_GAMMA = 13.0 / 11.0
"""float: The adiabatic index of the vapour plume (Svetsov 2007)."""

ANGLE_FACTOR = 5.0 + math.pi - 4.0 / 3.0
"""float: The correction for the distribution of the impact angle (Svetsov 2007, Eq. 7 of Sakuraba et al.), 6.81."""

IMPACTOR_ESCAPE_COEFFICIENT = 0.07
"""float: The upper branch of Shuvalov's escaping share of the impactor, `0.07 (rho_t / rho_imp) (V / v_esc)` (Eq. 8 of Sakuraba et al.); the
other branch of the minimum is smaller only at low erosion power, so this is the upper bound (tier 2, second-hand)."""

STEPS = 400
"""int: Steps of the integration over the accreted mass."""


# --- how much accretes, how fast, what it brings --------------------------------------------------------------------------------------

def late_accreted_percent_mean(mass_earth):
    """The mean mass accreted after the last giant impact, weight percent: `0.7 (M / M_earth)^(-1/3)`."""
    return LATE_ACCRETED_PERCENT_EARTH * mass_earth ** LATE_ACCRETED_MASS_EXPONENT


def late_accreted_fraction(mass_earth, constants):
    """
    The mass a rocky planet accretes after its last giant impact, as a fraction of its own mass.

    Args:
        mass_earth (float): The planet's mass, Earth masses.
        constants (dict): The design constants drawn so far; this adds `late_accreted_percent`.

    Returns:
        float: The fraction (not percent).
    """
    percent = min(late_accreted_percent_mean(mass_earth) * math.exp(draw.gauss(0.0, LATE_ACCRETED_LN_SIGMA)), LATE_ACCRETED_PERCENT_MAX)
    constants["late_accreted_percent"] = percent
    return percent / 100.0


def veneer_volatiles(constants):
    """`{species: mass fraction}` of the impactors, the reference times a drawn scale (recorded under `veneer_volatile_scale`)."""
    scale = bracket(constants, "veneer_volatile_scale", VENEER_VOLATILE_SCALE)
    return {species: share * scale for species, share in VENEER_VOLATILES.items()}


def impact_speed_ms(mass_kg, radius_m, orbital_speed_ms):
    """The mean speed at which the leftover planetesimals hit: `sqrt(v_esc^2 + (0.6 v_orbit)^2)`."""
    escape = math.sqrt(2.0 * G * mass_kg / radius_m)
    return detmath.hypot(escape, V_INF_OVER_V_ORBIT * orbital_speed_ms)


def impactor_escape_share(speed_ms, escape_ms):
    """Shuvalov's share of the impactor's vapour that leaves, upper branch (`IMPACTOR_ESCAPE_COEFFICIENT`)."""
    return min(IMPACTOR_ESCAPE_COEFFICIENT * TARGET_DENSITY / IMPACTOR_DENSITY * speed_ms / escape_ms, 1.0)


# --- Svetsov's erosion efficiency -----------------------------------------------------------------------------------------------------

def _plume_fraction(x_lower):
    """`int_x^1 (1 - x^2)^k dx / int_0^1 (1 - x^2)^k dx` for `k = 5`, in closed form (the integrand is a polynomial)."""
    def antiderivative(x):
        # (1 - x^2)^5 = 1 - 5x^2 + 10x^4 - 10x^6 + 5x^8 - x^10
        return x - 5.0 * x ** 3 / 3.0 + 2.0 * x ** 5 - 10.0 * x ** 7 / 7.0 + 5.0 * x ** 9 / 9.0 - x ** 11 / 11.0
    total = antiderivative(1.0)
    return (total - antiderivative(min(max(x_lower, 0.0), 1.0))) / total


def erosion_efficiency(density_kg_m3, scale_height_m, diameter_m, speed_ms, escape_ms):
    """
    The air mass one impactor removes over its own mass (Eq. 5 and 6 of Sakuraba et al. 2019).

    Args:
        density_kg_m3 (float): The air's density at the ground.
        scale_height_m (float): Its scale height.
        diameter_m (float): The impactor's diameter.
        speed_ms (float): The entry speed.
        escape_ms (float): The planet's escape speed.

    Returns:
        float: `eta`, the ejected air mass per impactor mass.
    """
    rho, h, d = density_kg_m3, scale_height_m, diameter_m
    if rho <= 0.0:
        return 0.0
    rho_i = IMPACTOR_DENSITY
    slowing = rho / rho_i * (h / d + 2.0 * h ** 2 * math.sqrt(rho) / (0.75 * d ** 2 * math.sqrt(rho_i)))
    ground = speed_ms * math.exp(-slowing)
    plume = ground * math.sqrt(4.0 * SVETSOV_GAMMA / (SVETSOV_GAMMA - 1.0))
    x_lower = min(1.0, escape_ms / plume) if plume > 0.0 else 1.0
    column = 2.0 * h / d + 16.0 * h ** 2 * math.sqrt(rho) / (3.0 * d ** 2 * math.sqrt(rho_i)) + 16.0 * h ** 3 * rho / (d ** 3 * rho_i)
    return 3.0 * rho / (4.0 * rho_i) * column * _plume_fraction(x_lower) * ANGLE_FACTOR


def _diameter_grid_m(bins=48):
    """Impactor diameters (m) spaced evenly in the logarithm, and their mass weights (mass per unit diameter is flat, so the weight goes as D)."""
    lo, hi = (math.log(d_km * 1000.0) for d_km in DIAMETER_RANGE_KM)
    grid = [math.exp(lo + (hi - lo) * (i + 0.5) / bins) for i in range(bins)]
    total = math.fsum(grid)
    return grid, [d / total for d in grid]


_GRID, _WEIGHTS = _diameter_grid_m()


def erosion_per_impactor_mass(density_kg_m3, scale_height_m, speed_ms, escape_ms):
    """`eta` averaged over the size distribution (weighted by impactor mass)."""
    return math.fsum(w * erosion_efficiency(density_kg_m3, scale_height_m, d, speed_ms, escape_ms) for d, w in zip(_GRID, _WEIGHTS))


# --- the balance ------------------------------------------------------------------------------------------------------------------------

def veneer_balance(mass_kg, radius_m, temperature_k, accreted_fraction, volatiles, orbital_speed_ms, *, caps_pa=None, air_kg=None, steps=STEPS):
    """
    What is in the air, kept at the surface and lost after the planet has accreted `accreted_fraction` of its mass in leftover planetesimals.

    Args:
        mass_kg (float): The planet's mass.
        radius_m (float): Its radius.
        temperature_k (float): The temperature of the air during the accretion (sets the scale height).
        accreted_fraction (float): The impactor mass over the planet's (`late_accreted_fraction`).
        volatiles (dict): Species (`SPECIES`) -> mass fraction in the impactors (`veneer_volatiles`).
        orbital_speed_ms (float): The planet's orbital speed.
        caps_pa (dict, optional): Species -> the largest partial pressure the surface lets the air keep; the excess is kept in the ground,
            the ocean or an ice cap. A species not in it is not capped.
        air_kg (dict, optional): Species -> the air the planet starts with.
        steps (int): Integration steps.

    Returns:
        dict: `air_kg`, `kept_kg` (condensed or dissolved), `lost_kg` (eroded or escaped with the impactor's vapour) and `delivered_kg`, each a
        dict of species; `eta` and `escape_share` of the last step, `speed_ms`.
    """
    caps_pa = caps_pa or {}
    gravity = G * mass_kg / radius_m ** 2
    escape = math.sqrt(2.0 * G * mass_kg / radius_m)
    speed = impact_speed_ms(mass_kg, radius_m, orbital_speed_ms)
    zeta = impactor_escape_share(speed, escape)
    air = {s: float((air_kg or {}).get(s, 0.0)) for s in SPECIES}
    kept = dict.fromkeys(SPECIES, 0.0)
    lost = dict.fromkeys(SPECIES, 0.0)
    delivered = {s: accreted_fraction * mass_kg * volatiles.get(s, 0.0) for s in SPECIES}
    step_mass = accreted_fraction * mass_kg / steps
    area = 4.0 * math.pi * radius_m ** 2
    eta = 0.0
    for _ in range(steps):
        total = math.fsum(air.values())
        if total > 0.0:
            moles = math.fsum(air[s] / MOLAR_MASS[s] for s in SPECIES)
            mean_mass = total / moles * AMU_KG   # kg per molecule (g/mol times the atomic mass unit)
            density = gravity * total / area * mean_mass / (K_BOLTZMANN * temperature_k)
            scale = K_BOLTZMANN * temperature_k / (mean_mass * gravity)
            eta = erosion_per_impactor_mass(density, scale, speed, escape)
        else:
            eta = 0.0
        for s in SPECIES:
            gain = (1.0 - zeta) * volatiles.get(s, 0.0) * step_mass
            lost[s] += zeta * volatiles.get(s, 0.0) * step_mass
            erode = min(eta * step_mass * (air[s] / total if total > 0.0 else 0.0), air[s])
            lost[s] += erode
            air[s] += gain - erode
            cap = caps_pa.get(s)
            if cap is not None:
                held = cap * area / gravity
                if air[s] > held:
                    kept[s] += air[s] - held
                    air[s] = held
    return {"air_kg": air, "kept_kg": kept, "lost_kg": lost, "delivered_kg": delivered, "eta": eta, "escape_share": zeta, "speed_ms": speed}
