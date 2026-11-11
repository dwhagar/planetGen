# planetgen/physics/disk_drift.py

"""
Pebble drift and evaporation at the ice lines of a protoplanetary disk (GEN.1473; docs/design/composition/disk-drift-evaporation.md)
====================================================================================================================================

`disk_condensation` holds a disk that stays as it was born: the gas at a place has the composition the temperature there gives. Pebbles do not stay: they grow, drift inward
faster than the gas, and evaporate where they cross each ice line, so the gas inside a line gets that line's ice at a rate set by the pebble flux over the gas flux
[R5208, R5209, R5210]. One number does the work, the **enrichment ratio** `E` of a species: the flux of its ice that pebbles carry in, over the flux the static budget
would carry with the gas. `E = 1` is the static budget (`disk_condensation`, GEN.1363) exactly; `E > 1` while the pebbles last; `E = 0` before the pebble front reaches the
line and after the pebbles are gone.

1. **The clocks.** The dust-growth front moves outward, `r_g = FRONT_AU_AT_1_MYR x (t / 1 Myr)^(2/3) x (M / Msun)^(1/3)` ([R5211: Eq. 10]); a line at `r_L` sees its pebbles from the time
   the front reaches it (`arrival_myr`) to the time the front leaves the disk's edge (`pebble_end_myr`: the drawn snow line times
   `tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER`). The gas needs time to carry the vapor from the line to a place inside it (`lag_myr`).
2. **The ratio.** `E(t) = Mdot_F(t) / (Mdot_g(t) x Z)` between the two times (`instant_ratio`): the pebble mass flux of [R5211: Eq. 14] over the viscous gas flux, in units of the
   disk's solids over gas `Z`. Both fluxes are proportional to the disk's mass (`tuning.DISK_MASS_STELLAR_MASS_EXPONENT`), so the ratio does not depend on it; the star's mass
   enters through the front and through the times only.
3. **The gas at a place.** The gas is the supply (the gas at the disk's edge) plus each ice whose line lies between the edge and the place, weighted by its `E` averaged over
   the lag (`ice_weights`). `disk_condensation._gas_and_solid_atoms(..., ice_weights=...)` takes the weights; weights of 1 give the static gas.
4. **A giant.** It gathers its gas when its core is built from the pebbles that pass its place (`gathering_age_myr`), which needs no new random draw: the star's mass and
   metals, its stored snow line and the giant's formation temperature give everything (`giant_ice_weights`).

Limits. The solids at a place stay the static ones; only the gas moves (DE8: the supply is the gas at the edge and each ice between the edge and the place is added with its weight; a
species is carried only by pebbles that began beyond its line, DE7). The pebble flux ends as a step at the disk's edge (DE3). Not modelled (DE10): a formed giant that blocks the
pebbles (the isolation mass), planets that eat pebbles on the way, the soot line, photoevaporation of the gas, and the cold finger
(the ice-to-rock ratio just outside a line; `COLD_FINGER_LEAK` and `cold_finger_enrichment` report it, nothing reads it yet). The nitrogen, sulfur and noble gases of
`disk_condensation.envelope_enrichments` keep their static values. Stand-ins are rows DE1 to DE10 of the registry (section 6.242 of `docs/design/research/gap-filling-plan.md`).
"""

import math

from planetgen import tuning
from planetgen.physics import constants as physical_constants, disk_condensation, scatter
from planetgen.util.checks import finite_domain

MYR_YR = 1.0e6
"""float: Years in a megayear."""

PEBBLE_FLUX_EARTH_MASSES_PER_YR = 9.5e-5
"""float: The mass flux of pebbles through a Sun-like disk at 1 Myr, Earth masses per year, for solids over gas `REFERENCE_SOLIDS_OVER_GAS`. The flux follows `(Z / 0.01)^(5/3)` and
`(t / 1 Myr)^(-1/3)` (`PEBBLE_FLUX_Z_EXPONENT`, `PEBBLE_FLUX_AGE_EXPONENT`), and the disk's mass (`disk_mass_scale`) in the place of the paper's star-mass term, so that the disk's
mass cancels in the enrichment ratio (design note 2.1). The form and number are [R5211: Eq. 14] (9.5e-5 Earth masses per year for beta 500 g/cm2, Z0 0.01, a Sun-like star at 10^6 yr). [Stand-in DE1, tier T1, tweak L0.]"""

PEBBLE_FLUX_Z_EXPONENT = 5.0 / 3.0
"""float: The power of the disk's solids over gas in the pebble flux ([R5211: Eq. 14]). [Stand-in DE1, tier T1, tweak L0.]"""

PEBBLE_FLUX_AGE_EXPONENT = -1.0 / 3.0
"""float: The power of the age in the pebble flux ([R5211: Eq. 14]). [Stand-in DE1, tier T1, tweak L0.]"""

REFERENCE_SOLIDS_OVER_GAS = 0.01
"""float: The solids over gas the flux of `PEBBLE_FLUX_EARTH_MASSES_PER_YR` is for ([R5211: Sec. 2.3]). [Stand-in DE1, tier T1, tweak L0.]"""

FRONT_AU_AT_1_MYR = 122.0
"""float: The distance the dust-growth front has reached at 1 Myr around a Sun-like star, AU: `r_g = 122 AU x (t / 1 Myr)^(2/3) x (M / Msun)^(1/3)` ([R5211: Eq. 10]; the number follows from a growth
efficiency of 0.05 (the sticking efficiency 0.5 over the ISM grain growth ratio of 10) and a solids over gas of 0.01). [Stand-in DE2, tier T1, tweak L0.]"""

FRONT_AGE_EXPONENT = 2.0 / 3.0
"""float: The power of the age in the front's distance ([R5211: Eq. 10]). [Stand-in DE2, tier T1, tweak L0.]"""

FRONT_MASS_EXPONENT = 1.0 / 3.0
"""float: The power of the star's mass in the front's distance ([R5211: Eq. 10]). [Stand-in DE2, tier T1, tweak L0.]"""

DISK_MASS_SOLAR = 0.05
"""float: The mass of a Sun-like star's disk at its birth, solar masses (scaled by `disk_mass_scale` for other stars). Standard disk theory; the number is chosen so the accretion rate at 1 to
2 Myr is the observed 1e-8 solar masses per year. [Stand-in DE4, tier T6, tweak L0.]"""

VISCOUS_TIME_MYR = 0.5
"""float: The viscous time of the disk, Myr: `Mdot_g = Mdot_0 x (1 + t / t_nu)^(-3/2)` with `Mdot_0 = Md_0 / (2 t_nu)`. [Stand-in DE4, tier T6, tweak L0.]"""

VISCOUS_GAS_EXPONENT = -1.5
"""float: The power of `(1 + t / t_nu)` in the gas flux (a viscosity linear in the radius). [Stand-in DE4, tier T6, tweak L0.]"""

ALPHA = 1.0e-3
"""float: The turbulence parameter of the disk, which sets the time the gas takes between two places. [Stand-in DE5, tier T6, tweak L0.]"""

LAG_MYR_PER_AU = 0.05
"""float: The time the gas takes to cross one AU inward at `REFERENCE_ALPHA`, Myr: `tau = (r_L - r) / (3 nu_1)` with `nu = nu_1 x (r / 1 AU)` (a steady disk). [Stand-in DE5, tier T6, tweak L0.]"""

REFERENCE_ALPHA = 1.0e-3
"""float: The `ALPHA` at which `LAG_MYR_PER_AU` holds; the lag scales as `REFERENCE_ALPHA / ALPHA`. [Stand-in DE5, tier T6, tweak L0.]"""

COLD_FINGER_LEAK = 0.9
"""float: The share of the vapor at a line that goes back out and freezes again just outside it; the ice there is enriched by `1 / (1 - leak)` and the net flux inward does not change
([R5209] abstract: the local vapor is enhanced by more than an order of magnitude). Reported by `cold_finger_enrichment`; nothing reads it yet. [Stand-in DE6, tier T6, tweak L0.]"""

CRITICAL_CORE_EARTH = 15.0
"""float: The core mass at which a giant starts to gather gas, Earth masses (the critical core, 10 to 20), the mass the pebbles must build before the gas is gathered. My derivation to avoid a new
draw. [Stand-in DE9, tier T6, tweak L0.]"""

INTEGRATION_INTERVALS = 32
"""int: The intervals of the Simpson rule that averages the ratio over a lag window (in the age to the 2/3, where the ratio is smooth)."""

_EARTH_MASSES_PER_SOLAR_MASS = physical_constants.SOLAR_MASS_TO_KG / physical_constants.EARTH_MASS_TO_KG

_SOLIDS_OVER_GAS_SUN = sum(disk_condensation.masses_percent(15.0).values()) / 100.0


def solids_over_gas(metal_ratio=1.0):
    """The disk's solids over gas (rock, grains and every ice, by mass): 1.43 percent for the Sun's composition, scaled by the star's metals over the Sun's (`metal_ratio`)."""
    if not metal_ratio > 0.0:
        raise ValueError("the metal ratio must be positive")
    return _SOLIDS_OVER_GAS_SUN * metal_ratio


def disk_mass_scale(mass_solar):
    """The disk's mass over a Sun-like star's, `(M / Msun)^DISK_MASS_STELLAR_MASS_EXPONENT` (it multiplies both fluxes and cancels in `instant_ratio`)."""
    return mass_solar ** tuning.DISK_MASS_STELLAR_MASS_EXPONENT


def edge_au(snow_line_au):
    """The disk's outer edge, AU: `tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER` snow lines."""
    return tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER * snow_line_au


def line_au_of(species, snow_line_au):
    """The distance of the line of `species` (`disk_condensation.LINES_K`), AU, from the star's drawn snow line (the water line): a disk's temperature falls as the square root of the distance."""
    if species not in disk_condensation.LINES_K:
        raise ValueError(f"no line for {species!r}")
    return snow_line_au * (disk_condensation.LINES_K["H2O"] / disk_condensation.LINES_K[species]) ** 2


def place_au(temperature_k, snow_line_au):
    """The distance, AU, at which the disk is at `temperature_k` (the inverse of the square-root law `line_au_of` uses)."""
    if not temperature_k > 0.0:
        raise ValueError("the temperature must be positive")
    return snow_line_au * (disk_condensation.LINES_K["H2O"] / temperature_k) ** 2


def edge_temperature_k():
    """The temperature at the disk's edge, K (the water line over the square root of the edge's multiple of the snow line): 42.9 K."""
    return disk_condensation.LINES_K["H2O"] / math.sqrt(tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER)


def front_au(age_myr, mass_solar=1.0):
    """The distance the dust-growth front has reached at `age_myr`, AU ([R5211: Eq. 10])."""
    if not age_myr >= 0.0 or not mass_solar > 0.0:
        raise ValueError("the age must be zero or more and the mass positive")
    return FRONT_AU_AT_1_MYR * age_myr ** FRONT_AGE_EXPONENT * mass_solar ** FRONT_MASS_EXPONENT


def arrival_myr(distance_au, mass_solar=1.0):
    """The age at which the dust-growth front reaches `distance_au` (the inverse of `front_au`), Myr."""
    if not distance_au >= 0.0 or not mass_solar > 0.0:
        raise ValueError("the distance must be zero or more and the mass positive")
    return (distance_au / (FRONT_AU_AT_1_MYR * mass_solar ** FRONT_MASS_EXPONENT)) ** (1.0 / FRONT_AGE_EXPONENT)


def pebble_end_myr(snow_line_au, mass_solar=1.0):
    """
    The age at which the front leaves the disk's edge and the pebbles are gone, Myr: 0.2 for the Sun (edge 41 AU). The flux ends as a step there; speculative, the simplest form that conserves
    the dust (a real disk loses its pebbles over a few Myr). [Stand-in DE3, tier T6, tweak L0.]
    """
    return arrival_myr(edge_au(snow_line_au), mass_solar)


def pebble_flux_earth_per_yr(age_myr, mass_solar=1.0, metal_ratio=1.0):
    """The pebble mass flux, Earth masses per year, at `age_myr` ([R5211: Eq. 14], with the disk-mass term of `disk_mass_scale`)."""
    if not age_myr > 0.0:
        raise ValueError("the age must be positive")
    z = solids_over_gas(metal_ratio)
    return (PEBBLE_FLUX_EARTH_MASSES_PER_YR * (z / REFERENCE_SOLIDS_OVER_GAS) ** PEBBLE_FLUX_Z_EXPONENT * age_myr ** PEBBLE_FLUX_AGE_EXPONENT
            * disk_mass_scale(mass_solar))


def gas_flux_solar_per_yr(age_myr, mass_solar=1.0):
    """The viscous gas flux, solar masses per year: `Mdot_0 (1 + t / t_nu)^(-3/2)` with `Mdot_0 = Md_0 / (2 t_nu)`, `Md_0` scaled by `disk_mass_scale`."""
    if not age_myr >= 0.0:
        raise ValueError("the age must be zero or more")
    initial = DISK_MASS_SOLAR * disk_mass_scale(mass_solar) / (2.0 * VISCOUS_TIME_MYR * MYR_YR)
    return initial * (1.0 + age_myr / VISCOUS_TIME_MYR) ** VISCOUS_GAS_EXPONENT


def instant_ratio(age_myr, mass_solar=1.0, metal_ratio=1.0):
    """
    The pebble flux over the gas flux in units of the solids over gas, at `age_myr`: `Mdot_F / (Mdot_g x Z)`. The disk's mass is in both fluxes and cancels. Not gated: see
    `enrichment_ratio` for the ratio between the arrival of a line's pebbles and their end. About 2.3 at 0.05 Myr for the Sun.
    """
    earth_per_solar = 1.0 / _EARTH_MASSES_PER_SOLAR_MASS
    flux = pebble_flux_earth_per_yr(age_myr, mass_solar, metal_ratio) * earth_per_solar
    return flux / (gas_flux_solar_per_yr(age_myr, mass_solar) * solids_over_gas(metal_ratio))


def lag_myr(from_au, to_au):
    """The time the gas takes to carry vapor from the line at `from_au` to a place at `to_au` inside it, Myr (DE5); zero for a place at or outside the line."""
    if not from_au >= 0.0 or not to_au >= 0.0:
        raise ValueError("the distances must be zero or more")
    return LAG_MYR_PER_AU * (REFERENCE_ALPHA / ALPHA) * max(from_au - to_au, 0.0)


def _mean_ratio(start_myr, end_myr, window_myr, mass_solar, metal_ratio):
    """The integral of `instant_ratio` from `start_myr` to `end_myr` over `window_myr` (the boxcar mean of the gated ratio over the window), by the Simpson rule in `u = t^(2/3)`."""
    if end_myr <= start_myr:
        return 0.0
    exponent = 1.0 / FRONT_AGE_EXPONENT
    u_start, u_end = start_myr ** FRONT_AGE_EXPONENT, end_myr ** FRONT_AGE_EXPONENT
    steps = INTEGRATION_INTERVALS

    def integrand(u):
        t = u ** exponent
        return instant_ratio(t, mass_solar, metal_ratio) * exponent * u ** (exponent - 1.0)

    width = (u_end - u_start) / steps
    total = integrand(u_start) + integrand(u_end)
    for index in range(1, steps):
        total += (4.0 if index % 2 else 2.0) * integrand(u_start + index * width)
    return total * width / 3.0 / window_myr


@finite_domain()
def enrichment_ratio(species, age_myr, snow_line_au, mass_solar=1.0, metal_ratio=1.0, to_au=None):
    """
    The enrichment ratio `E` of `species` at `age_myr`: the flux of its ice that the pebbles carry in over the flux the static budget carries with the gas.

    One, the static budget, for gas that has not crossed the line yet: with `to_au`, when `age_myr` is under the time the gas needs to carry vapor from the line to the place
    (`lag_myr`), the gas at the place began inside the line with its vapor (GEN.1517; the start state of the disk is the solar gas, not gas with the ices frozen out). Past that
    the gas crossed the line at `age_myr - lag`: zero when the line is not inside the disk's edge, before the front reaches the line (`arrival_myr`) and from the time the pebbles end (`pebble_end_myr`); between,
    `instant_ratio`. With `to_au` (a place inside the line) the value is the mean over the window `[age - lag_myr, age]` that the gas needs to carry the vapor from the line to the
    place, the window cut at the start of the disk (design note 2.1).

    Args:
        species (str): A key of `disk_condensation.LINES_K`.
        age_myr (float): The disk's age, Myr.
        snow_line_au (float): The star's drawn snow line, AU.
        mass_solar (float): The star's mass.
        metal_ratio (float): The star's metals over the Sun's (the disk's solids scale with it).
        to_au (float): The place the vapor is carried to; `None` for no lag.
    """
    if not age_myr >= 0.0:
        raise ValueError("the age must be zero or more")
    line = line_au_of(species, snow_line_au)
    if line >= edge_au(snow_line_au):
        return 0.0
    lag = 0.0 if to_au is None else lag_myr(line, to_au)
    if age_myr < lag:
        return 1.0                                  # the gas at the place has not crossed the line yet: it is the disk's own gas of the start, vapor and all (GEN.1517)
    first = arrival_myr(line, mass_solar)
    last = pebble_end_myr(snow_line_au, mass_solar)
    if first >= last:
        return 0.0
    if lag <= 0.0 or age_myr <= 0.0:
        return instant_ratio(age_myr, mass_solar, metal_ratio) if first <= age_myr < last else 0.0
    window_start = max(age_myr - lag, 0.0)
    return _mean_ratio(max(window_start, first), min(age_myr, last), age_myr - window_start, mass_solar, metal_ratio)


@finite_domain()
def ice_weights(temperature_k, age_myr, snow_line_au, mass_solar=1.0, metal_ratio=1.0):
    """
    The weight `E` of each ice whose line lies between the disk's edge and `temperature_k`, as `{species: E}`, for the gas at the place where the disk is at `temperature_k` at
    `age_myr` (the lag from each line to the place included). These are the ices that evaporated on the way in; the supply (the gas at the edge) is not weighted. Empty for a place
    at or beyond the edge. `disk_condensation.envelope_c_over_o(..., ice_weights=...)` takes the result.
    """
    place = place_au(temperature_k, snow_line_au)
    edge = edge_temperature_k()
    weights = {}
    for species, line in disk_condensation.LINES_K.items():
        if edge < line <= temperature_k:
            weights[species] = enrichment_ratio(species, age_myr, snow_line_au, mass_solar, metal_ratio, to_au=place)
    return weights


@finite_domain()
def gathering_age_myr(temperature_k, snow_line_au, mass_solar=1.0, metal_ratio=1.0):
    """
    The age at which a giant that formed where the disk was at `temperature_k` gathers its gas, Myr (DE9): when the pebbles that pass its place have built `CRITICAL_CORE_EARTH`
    from the front's arrival there, capped at the disk's lifetime (`scatter.DISK_LIFETIME_YR`) and at the end of the pebbles (a core the pebbles cannot build in time gathers
    its gas as the pebbles end and sees the pebble-depleted gas).
    """
    place = place_au(temperature_k, snow_line_au)
    arrived = arrival_myr(place, mass_solar)
    last = pebble_end_myr(snow_line_au, mass_solar)
    cap = min(scatter.DISK_LIFETIME_YR / MYR_YR, last)
    if arrived >= cap:
        return cap
    z = solids_over_gas(metal_ratio)
    # mass(t) = K x (t^(2/3) - t0^(2/3)) with the flux K' t^(-1/3): the integral of t^(-1/3) is 1.5 t^(2/3)
    coefficient = (PEBBLE_FLUX_EARTH_MASSES_PER_YR * (z / REFERENCE_SOLIDS_OVER_GAS) ** PEBBLE_FLUX_Z_EXPONENT * disk_mass_scale(mass_solar) * MYR_YR
                   / (PEBBLE_FLUX_AGE_EXPONENT + 1.0))
    needed = CRITICAL_CORE_EARTH / coefficient + arrived ** (PEBBLE_FLUX_AGE_EXPONENT + 1.0)
    age = needed ** (1.0 / (PEBBLE_FLUX_AGE_EXPONENT + 1.0))
    return min(age, cap)


def giant_ice_weights(temperature_k, snow_line_au, mass_solar=1.0, metal_ratio=1.0):
    """The `ice_weights` of the gas a giant gathers where the disk was at `temperature_k`, at the age it gathers it (`gathering_age_myr`); one function of the star's mass and metals and its stored snow line, no draw."""
    age = gathering_age_myr(temperature_k, snow_line_au, mass_solar, metal_ratio)
    return ice_weights(temperature_k, age, snow_line_au, mass_solar, metal_ratio)


def cold_finger_enrichment(share_of_flux_that_is_ice=1.0):
    """
    How many times the ice-to-rock ratio of the solids just outside a line rises over the static one when `COLD_FINGER_LEAK` of the vapor goes back out: `1 / (1 - leak)`, cut by
    the share of the pebble flux that is ice (`1 + (1 / (1 - leak) - 1) x share`). Reported, not read (design note 2.4).
    """
    if not 0.0 <= share_of_flux_that_is_ice <= 1.0:
        raise ValueError("the share must be between 0 and 1")
    return 1.0 + (1.0 / (1.0 - COLD_FINGER_LEAK) - 1.0) * share_of_flux_that_is_ice
