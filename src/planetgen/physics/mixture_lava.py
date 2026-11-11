# planetgen/physics/mixture_lava.py

"""
The vapor over a molten surface: lava worlds and the silicate exosphere of a white-dwarf or pulsar planet (GEN.499, GEN.456, GEN.479)
======================================================================================================================================

A surface hot enough to be molten (`MELT_MIN_K`) holds a vapor whose species are fixed by the melt's own chemistry, not by an element
budget: the melt is an infinite reservoir, and every gas that its atoms can make stands at the pressure the melt's chemical potentials
give it. The melt is an ideal solution of its oxides ([R3942], 2009; [R3059] treat the BSE the same way
to first order) whose liquid free energies are NASA Glenn's (`physics/species_thermo.py`, extended below their melting points as
supercooled liquids), and the oxygen fugacity of the vapor is held at the redox the recipe tables drew (a1, the magma-ocean
stage, or a4, the white-dwarf debris), in units of the iron-wuestite buffer. With each element's potential

    m_M = [ g(M_x O_y, liquid) + ln a(M_x O_y) - y m_O ] / x         m_O = (g(O2) + ln fO2) / 2
    m_Fe = g(Fe, liquid) + ln a(FeO) - (dIW / 2) ln 10               (the buffer fixes the FeO)

the partial pressure of any gas j of those elements follows with no iteration: `ln p_j = sum_el nu_el m_el - g(j)`, the standard
state being 1 bar. The mole fractions are the partial pressures over their sum; the pressure the vapor sums to is the vapor pressure
over the melt ([R2380]'s `p_eq_bar` of the recipe row `i6`, which this reproduces within a factor 10; it is the generator's
to use, not this module's). The oxides are all in the vapor's elements, aluminium included (GEN.818): its gases are Al, AlO, Al2O, AlO2 and Al2O2 of the
species table, from an Al2O3 activity that is its mole fraction (`ALUMINA_ACTIVITY_COEFFICIENT`).

The melt is the bulk silicate Earth's oxides (McDonough and Sun 1995, TEST composition, the only quantitative anchor of part 05
section 4) with each oxide's mass share scattered by `COMPOSITION_SCATTER` dex (design: no body's mantle is the Earth's; the bulk
composition of lane J replaces it when stored). A body whose recipe drew the volatile-bearing branch (`i6`) has a CO2 and CO gas of its
own over the vapor, at the surface pressure: carbon with the same held oxygen, by the equilibrium solver (the 55 Cnc e fit of Hu et al.
2024, TEST).
"""

import math

from planetgen.physics import elements, gas_equilibrium, oxygen_fugacity as fo2, species_thermo
from planetgen.physics.mixture_base import MixtureError, bracket, finish, require
from planetgen.util import draw

MELT_MIN_K = 1400.0
"""float: The least surface temperature of a molten surface ([R2380]'s range starts here); a hotter-zone body below it
outgasses and is not a lava world."""

MELT_MAX_K = 6000.0
"""float: The hottest temperature the vapor over a melt is worked out at (GEN.1294): the top of NASA Glenn's fits of the liquid oxides and of
47 of the 84 gases of the species table (`species_thermo`), above which `species_thermo` raises. A surface hotter than this (a dayside of a rock
next to a hot star, or a white-dwarf planet's exosphere) has its vapor worked out at this temperature, and the mixture reports the surface's own.
A stand-in: registry row AM95 (tier 6, tweak level 0)."""

SUPERCOOLED_BELOW_K = 1000.0
"""float: How far below its melting point a liquid's fit is used (a lava world's surface is at least `MELT_MIN_K`)."""

OXIDES = {
    "SiO2": (45.0, "Si", 1, 2), "MgO": (37.8, "Mg", 1, 1), "FeO": (8.05, "Fe", 1, 1), "Al2O3": (4.45, "Al", 2, 3),
    "CaO": (3.55, "Ca", 1, 1), "Na2O": (0.36, "Na", 2, 1), "TiO2": (0.201, "Ti", 1, 2), "K2O": (0.029, "K", 2, 1),
}
"""dict: Oxide -> (weight percent in the bulk silicate Earth, McDonough and Sun 1995; its metal; metal atoms per formula unit; oxygen
atoms per formula unit)."""

COMPOSITION_SCATTER = 0.15
"""float: The dex each oxide's share of the melt is scattered by (design)."""

LIQUID = {"SiO2": "SiO2(L)", "MgO": "MgO(L)", "Na2O": "Na2O(L)", "K2O": "K2O(L)", "CaO": "CaO(L)", "TiO2": "TiO2(L)",
          "Al2O3": "AL2O3(L)"}
"""dict: Oxide -> the liquid phase of the species table."""

ELEMENTS = ("O", "Na", "K", "Mg", "Si", "Fe", "Ca", "Ti", "Al")
"""tuple: The elements the vapor is made of (the melt's)."""

ALKALI_ACTIVITY_COEFFICIENT = (3e-9, 3e-8, "log")
"""tuple: The activity coefficient of Na2O and K2O in the silicate melt (an ideal solution has 1; the alkali oxides are held in the
silicate network and their activities are many orders of magnitude below their mole fractions). A stand-in: this bracket is the one
that reproduces [R2380]'s vapor pressure over a BSE melt (2e-4 bar at 2000 K, 8e-3 bar at 2500 K; READ) and the sequence
of [R3059] (sodium first, SiO, Mg and Fe from the 2000 to 2500 K range up; READ, tables not reached), which the ideal
solution misses by a factor 1e4 at 2000 K. The research item that replaces it is the melt activity model of GEN.499."""

ALUMINA_ACTIVITY_COEFFICIENT = 1.0
"""float: The activity coefficient of Al2O3 in the silicate melt, so that its activity is its mole fraction (the liquid at the melt's
temperature is the standard state). A stand-in (tier 3, tweak level 0): [R1341] compare activity-coefficient
models of the bulk silicate Earth melt with an ideal solution and find them agreeing except for the alkali elements, with aluminium
not named among the sensitive ones and its condensation temperature the same for the dry and the wet melt (read in the abstract only;
the paper and its Al2O3 value are not opened). The alkalis are the exception that `ALKALI_ACTIVITY_COEFFICIENT` carries. The research
item that replaces it is the melt activity model of GEN.499."""

OXYGEN_SHARE_MAX = 0.99
"""float: The largest share of the gas pressure over a melt that its held oxygen (O2 and O at the buffer's fugacity) may take in the carbon solve
of a volatile-bearing body (GEN.1296). The buffered oxygen fills the pressure when it reaches it: the solve then has no room for carbon and does not
converge (checked from 0.8 to 1.2 of the pressure, at six temperatures between 2600 and 6000 K and three redox states: it converges for every
share under 1 and fails for every share of 1 and over). A surface pressure lower than the oxygen holds is raised to what holds it, as the gas over a
melt cannot be thinner than the melt's own oxygen. A stand-in: registry row AM96 (tier 6, tweak level 0)."""

VOLATILE_BEARING_PER_CARBON = 1.0
"""float: Moles of carbon the volatile-bearing gas is computed for (only its CO2 to CO ratio matters)."""


def melt_composition(constants, oxides_wt_percent=None):
    """
    The oxides' mole fractions of the melt.

    Args:
        constants (dict): Receives `melt_mole_fractions`.
        oxides_wt_percent (dict, optional): The stored surface composition, oxide -> weight percent (any of `OXIDES`, the rest are
            zero); without it the bulk silicate Earth's weight shares scattered by `COMPOSITION_SCATTER` dex are drawn.
    """
    moles = {}
    for oxide, (weight_percent, metal, atoms, oxygen) in OXIDES.items():
        mass = atoms * elements.element(metal).atomic_weight + oxygen * elements.element("O").atomic_weight
        if oxides_wt_percent is not None:
            moles[oxide] = max(oxides_wt_percent.get(oxide, 0.0), 1e-12) / mass
            continue
        scatter = 10.0 ** draw.gauss(0.0, COMPOSITION_SCATTER)
        moles[oxide] = weight_percent * scatter / mass
    total = sum(moles.values())
    fractions = {oxide: value / total for oxide, value in moles.items()}
    constants["melt_mole_fractions"] = dict(fractions)
    return fractions


def element_potentials(activities, temperature_k, redox_diw):
    """
    The chemical potential over R T of each element of the vapor, from the melt's oxide activities (see the module text).

    Returns:
        dict: Element symbol -> m over R T, standard state 1 bar.
    """
    t = temperature_k
    m_oxygen = 0.5 * (species_thermo.g_over_rt("O2", t) + fo2.log_fo2(redox_diw, t) * math.log(10.0))
    potentials = {"O": m_oxygen}
    for oxide, (_, metal, atoms, oxygen) in OXIDES.items():
        if metal == "Fe":
            continue
        g_liquid = species_thermo.SPECIES[LIQUID[oxide]].supercooled_g_over_rt(t, SUPERCOOLED_BELOW_K)
        potentials[metal] = (g_liquid + math.log(activities[oxide]) - oxygen * m_oxygen) / atoms
    iron = species_thermo.SPECIES["Fe(L)"].supercooled_g_over_rt(t, SUPERCOOLED_BELOW_K)
    potentials["Fe"] = iron + math.log(activities["FeO"]) - 0.5 * redox_diw * math.log(10.0)
    return potentials


def vapor_pressures(potentials, temperature_k):
    """Species name -> partial pressure, bar, of every gas of the species table made of `ELEMENTS` whose fit covers the temperature."""
    pressures = {}
    for name in species_thermo.GASES:
        row = species_thermo.SPECIES[name]
        if not set(row.elements) <= set(ELEMENTS):
            continue
        low, high = row.range_k
        if not low <= temperature_k <= high:
            continue
        log_p = sum(count * potentials[symbol] for symbol, count in row.elements.items()) - row.g_over_rt(temperature_k)
        pressures[name] = math.exp(max(log_p, -700.0))
    return pressures


CARBON_GAS_PRESSURE_STEP = 1.25
"""float: The factor the gas pressure of the carbon solve rises by each time the solve does not converge (GEN.1513)."""
CARBON_GAS_PRESSURE_STEPS = 8
"""int: How many times the carbon solve is tried again at a higher pressure."""


def carbon_gas(chemistry_k, gas_pressure, log_fugacity):
    """
    `(Equilibrium, gas_pressure)` of carbon over a melt with its oxygen held at `log_fugacity`. When the buffered oxygen takes nearly all of the pressure the
    solve has no room for the carbon and does not converge (GEN.1296 keeps the oxygen share under `OXYGEN_SHARE_MAX` of the melt's own vapor; the buffer's
    fugacity at 6000 K can be higher, 1.3e5 bar against 1.4e5 bar of gas in the zeta Ophiuchi example): the gas is then taken at a higher pressure, by
    `CARBON_GAS_PRESSURE_STEP` at a time (it converged within three steps from 4500 to 6000 K and fugacities of 1e3 to 1e7 bar), as the share rule does.

    Raises:
        MixtureError: It does not converge at `CARBON_GAS_PRESSURE_STEPS` higher pressures either.
    """
    for _ in range(CARBON_GAS_PRESSURE_STEPS + 1):
        try:
            return gas_equilibrium.solve({"C": VOLATILE_BEARING_PER_CARBON}, chemistry_k, gas_pressure, fixed_o2_log_fugacity=log_fugacity), gas_pressure
        except gas_equilibrium.EquilibriumError as failure:
            if "no convergence" not in str(failure):
                raise
            gas_pressure *= CARBON_GAS_PRESSURE_STEP
    raise MixtureError(f"the carbon gas over the melt does not settle at {chemistry_k:g} K with the oxygen held at 10^{log_fugacity:.2f} bar up to {gas_pressure:.4g} bar")


def vapor_mixture(subject, recipe, *, temperature_name="tp_k", oxides_wt_percent=None):
    """
    The vapor over `subject`'s molten surface.

    Args:
        subject (Subject): With the surface temperature under `temperature_name` (`ts_k` when it lacks it) and, for the
            volatile-bearing branch, `pressure_bar`.
        recipe (RecipeDraw): With the redox (`a1` or `a4`) and the lava rows drawn.
        temperature_name (str): The quantity the temperature is read from.
        oxides_wt_percent (dict, optional): The melt's stored oxide composition (see `melt_composition`).

    Raises:
        MixtureError: The surface is colder than `MELT_MIN_K`, or the redox is missing.
    """
    constants = {}
    temperature = subject.get(temperature_name)
    if temperature is None:
        require(subject, ("ts_k",))
        temperature = subject["ts_k"]
    if temperature < MELT_MIN_K:
        raise MixtureError(f"a surface at {temperature:g} K is not molten (the vapor over a melt is taken from {MELT_MIN_K:g} K)")
    chemistry_k = min(temperature, MELT_MAX_K)
    redox = recipe.values.get("a4", recipe.values.get("a1"))
    if redox is None:
        raise MixtureError("no redox row applies to the melt")
    fractions = melt_composition(constants, oxides_wt_percent)
    gamma = bracket(constants, "alkali_activity_coefficient", ALKALI_ACTIVITY_COEFFICIENT)
    coefficient = {"Na2O": gamma, "K2O": gamma, "Al2O3": ALUMINA_ACTIVITY_COEFFICIENT}
    activities = {oxide: value * coefficient.get(oxide, 1.0) for oxide, value in fractions.items()}
    potentials = element_potentials(activities, chemistry_k, redox)
    pressures = vapor_pressures(potentials, chemistry_k)
    constants["vapor_pressure_bar"] = sum(pressures.values())
    volatile_bearing = recipe.binaries.get("i6", (False, None))[0]
    if volatile_bearing:
        require(subject, ("pressure_bar",))
        gas_pressure = max(subject["pressure_bar"], (pressures.get("O2", 0.0) + pressures.get("O", 0.0)) / OXYGEN_SHARE_MAX)
        volatile, gas_pressure = carbon_gas(chemistry_k, gas_pressure, fo2.log_fo2(redox, chemistry_k))
        for name, fraction in volatile.mole_fractions.items():
            pressures[name] = pressures.get(name, 0.0) + fraction * gas_pressure
        constants["volatile_bearing"] = True
    total = sum(pressures.values())
    equilibrium = gas_equilibrium.Equilibrium(
        temperature_k=temperature, pressure_bar=total, moles=pressures,
        mole_fractions={name: value / total for name, value in pressures.items()}, total_moles=total, iterations=0)
    return finish(equilibrium, family="lava", drawn={group: drawn.name for group, drawn in recipe.groups.items()},
                  redox_diw=redox, inventory=dict(fractions), constants=constants)
