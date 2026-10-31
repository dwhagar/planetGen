# planetgen/physics/mixture_giant.py

"""
The air of a giant planet, an ice giant, a brown dwarf or a hot Jupiter (GEN.454; docs/design/composition/mixture-draw.md)
==============================================================================================================

The giant family of `mixture_draw.draw_mixture`, after parts 05 sections 2, 7 and 8 of the atmosphere research. A giant's
air is the gas of the host's element mixture, enriched in the metals by a drawn factor and given a drawn carbon to oxygen
ratio, in equilibrium at the temperature and pressure of its 0.1 to 1 bar layer, after three things the research names:

1. *Clouds take their elements out.* Iron, the magnesium silicates, titanium, vanadium, calcium and chromium condense at the
   temperatures of Visscher, Lodders and Fegley (2006, 2010; part 05 section 7.2) and leave the gas above the cloud deck,
   taking about a fifth of the oxygen with the silicates; sodium takes its sulfur as Na2S and potassium goes with it; ammonia
   and hydrogen sulfide form NH4SH; and whatever ice of the gas table has a saturation curve (`physics/saturation.py`: water,
   ammonia, methane, hydrogen sulfide, hydrogen cyanide) is capped at its saturation pressure, the excess raining out. The
   ratio of hydrogen to helium is the host's and does not move.
2. *Quench.* Carbon monoxide and methane stop exchanging where the mixing time of the atmosphere (`H^2 / Kzz`) falls below
   the chemical time (Zahnle and Marley 2014: `t_q1 = 1.5e-6 p^-1 m^-0.7 exp(42000/T)`, `t_q2 = 40 p^-2 exp(25000/T)`, seconds,
   p in bar, m the metallicity in solar units). The mixing ratio at the quench level, on an adiabat down from the layer, sets the
   share of the carbon that is carbon monoxide up in the layer, so a cold giant keeps some CO (1e-9 to 1e-5) and a hot one has
   no more methane than equilibrium gives. The generator evaluates the formula with the object's own gravity, metallicity and
   mixing, as part 05 section 7.5 asks; it does not interpolate between 500 K and 1100 K.
3. *Below 300 K the chemistry is the hydride limit.* The fits of the thermodynamic table start at 300 K for several gases; a
   colder layer is equilibrated at 300 K (every heavy element is its hydride there: water, methane, ammonia and hydrogen
   sulfide, which is what part 05 section 2.2 computes and the equilibrium reproduces), and the clouds then work at the layer's
   own temperature.

The recipe rows (`f1` to `f8`, `g1` to `g6`, `h1` to `h5`) give the regime and the draws the chemistry reads: the metallicity
`f5` (Thorngren et al. 2016 mass-metallicity relation with 0.3 dex scatter), the carbon to oxygen ratio `f7`, whether titanium
and vanadium oxides stay in the gas `h4`, whether sodium and potassium are in the gas `m7` and `m8`. `f8`, `g6` are the
recipe's summary of the carbon carrier; the chemistry here decides it by the quench.
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.physics import constants, gas_equilibrium, saturation, species_thermo, unit_registry as units
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.physics.mixture_base import MixtureError, bracket, finish, require

_mearth_to_kg = units.converter("Mearth", "kg")
_mjup_to_kg = units.converter("Mjup", "kg")

CATEGORIES = "efgh"
"""str: The recipe tables a giant reads (`e` for a sub-Neptune envelope, `f` giants, `g` brown dwarfs, `h` hot Jupiters), and the
rows `m7` and `m8` (sodium and potassium gas) it draws itself after them."""

LAYER_PRESSURE_BAR = {"giant": 0.3, "hot-jupiter": 0.3, "envelope": 0.3, "brown-dwarf": 1.0}
"""dict: The pressure (bar) of the layer the mixture is quoted for. A planet's is the middle (geometric) of the research's 0.1 to
1 bar layer; a brown dwarf's is the photosphere of its higher gravity, 1 bar (design: with it the CH4 to CO crossover of the
quench formula falls at 850 to 1100 K for Kzz 1e2 to 1e5, Zahnle and Marley's 1100 +/- 100 K)."""

CHEMISTRY_FLOOR_K = 300.0
"""float: A colder layer is equilibrated at this temperature (the hydride limit; the fits of several gases start here)."""

ADIABAT_EXPONENT = 2.0 / 7.0
"""float: The dry adiabat of a diatomic gas, d ln T / d ln P, for the temperature on the way down to the quench level."""

QUENCH_DEEPEST_BAR = 1.0e5
"""float: A quench level is looked for down to this pressure; a layer that does not equilibrate by here is quenched here."""

KZZ_CM2_S = {"giant": (1.0e6, 1.0e11, "log"), "envelope": (1.0e6, 1.0e11, "log"), "hot-jupiter": (1.0e7, 1.0e11, "log"),
             "brown-dwarf": (1.0e2, 1.0e5, "log")}
"""dict: The vertical mixing coefficient by object, cm^2/s. The hot Jupiter's is the bracket of row h4 (Spiegel et al. 2009); the
other planets' is that widened down to 1e6 and a brown dwarf's, whose gravity is 30 times higher, is the lowest decades of the
literature's 1e2 to 1e6; design, to be replaced by the research item that tabulates Kzz by object."""

PHOTOSPHERE_LOG_EPS = {"He": 10.914, "P": 5.41, "Ti": 4.97, "Kr": 3.12, "Xe": 2.22, "V": 3.90, "Cr": 5.62}
"""dict: Asplund, Amarsi and Grevesse 2021 photospheric log eps (H = 12) of the elements `tuning.STELLAR_SOLAR_LOG_EPS` lacks. He, P,
Ti, Kr and Xe are the part 05 section 2.1 table (READ); V and Cr are RECALLED from the same paper and carry the stand-in label."""

METALS = ("C", "N", "O", "Ne", "Na", "Mg", "Si", "S", "Ar", "K", "Ca", "Fe", "P", "Ti", "V", "Cr", "Kr", "Xe")
"""tuple: The elements the budget carries beyond hydrogen and helium (every one the thermodynamic table has a gas of)."""

ALKALI_ONSET_K = {"Na": 962.0, "K": 775.0}
"""dict: Where the sodium and potassium gas appear (the centres of rows `m7` and `m8`), for a body whose recipe has no such row."""

MIN_PARTIAL_BAR = 1e-40
"""float: A species with a partial pressure below this is not tested against its saturation curve."""

MIN_CONDENSED = 1e-10
"""float: Ices that took fewer moles than this per hydrogen atom are not listed as condensed."""


def _log_eps(symbol):
    return PHOTOSPHERE_LOG_EPS.get(symbol, tuning.STELLAR_SOLAR_LOG_EPS.get(symbol))


SOLAR_PER_HYDROGEN = {symbol: 10.0 ** (_log_eps(symbol) - 12.0) for symbol in ("He",) + METALS}
"""dict: Atoms per hydrogen atom of the reference (photospheric) scale, the zero point of [Fe/H] and of the metallicity."""


def element_budget(log_enrichment, c_over_o=None):
    """
    The atoms per hydrogen atom of the envelope: helium as the reference, every metal at `10 ** log_enrichment` times its reference
    abundance, with the carbon and oxygen re-split to `c_over_o` at their total (when given).

    Args:
        log_enrichment (float): log10 of the metals over their reference abundance (the metallicity of the planet plus the
            host's [Fe/H]).
        c_over_o (float, optional): The carbon to oxygen number ratio; default the reference's.

    Returns:
        dict: Element symbol -> atoms per hydrogen atom, `H` included (1).
    """
    scale = 10.0 ** log_enrichment
    budget = {"H": 1.0, "He": SOLAR_PER_HYDROGEN["He"]}
    for symbol in METALS:
        budget[symbol] = SOLAR_PER_HYDROGEN[symbol] * scale
    if c_over_o is not None:
        total = budget["C"] + budget["O"]
        budget["C"] = total * c_over_o / (1.0 + c_over_o)
        budget["O"] = total / (1.0 + c_over_o)
    return budget


# --- the clouds -----------------------------------------------------------------------------------------------------------

def condensation_temperatures(pressure_bar, log_iron):
    """
    The condensation temperatures (K) of the iron and magnesium silicate clouds at `pressure_bar`, from the Visscher, Lodders and
    Fegley fits of part 05 section 7.2: `10^4 / T = a - b log P - c [Fe/H]`, with `log_iron` the log10 of the metals over the reference.

    Returns:
        dict: `iron`, `forsterite` (Mg2SiO4), `enstatite` (MgSiO3).
    """
    log_p = math.log10(pressure_bar)
    return {
        "iron": 1.0e4 / (5.44 - 0.48 * log_p - 0.48 * log_iron),
        "forsterite": 1.0e4 / (5.89 - 0.37 * log_p - 0.73 * log_iron),
        "enstatite": 1.0e4 / (6.26 - 0.35 * log_p - 0.70 * log_iron),
    }


OXIDE_OXYGEN = {"Mg": 1.0, "Si": 2.0, "Ca": 1.0, "Ti": 2.0, "V": 1.0, "Cr": 1.5}
"""dict: Oxygen atoms per atom of each rock-forming element when it condenses as an oxide or silicate (MgO + SiO2 for the
magnesium silicates)."""


def remove_clouds(budget, temperature_k, pressure_bar, log_enrichment, *, keep_titanium_vanadium=True, keep_sodium=True,
                  keep_potassium=True):
    """
    The budget of the gas above the cloud decks at `temperature_k` and `pressure_bar`.

    Args:
        budget (dict): Atoms per hydrogen atom (`element_budget`).
        temperature_k (float): The layer.
        pressure_bar (float): The layer.
        log_enrichment (float): log10 of the metals over reference (for the `[X/H]` terms of the fits).
        keep_titanium_vanadium (bool): False when the cold trap holds the titanium and vanadium oxides out of the gas even where the
            layer is hotter than their condensation (row `h4`: not leaky); where the layer is colder they condense anyway.
        keep_sodium, keep_potassium (bool): Whether the alkali metal is in the gas (rows `m7` and `m8`: the Na2S and KCl clouds of
            part 05 section 7.2 as the recipe tabulates their onset); sodium that goes takes half its atoms of sulfur as Na2S.

    Returns:
        tuple: `(budget, condensed)`: the new budget and a dict of the cloud names that formed (`iron`, `forsterite`, `enstatite`,
            `oxides`, `sodium_sulfide`, `potassium_chloride`) -> their condensation temperature (0 where the recipe decided).
    """
    left = dict(budget)
    condensed = {}
    cond = condensation_temperatures(pressure_bar, log_enrichment)

    def take(symbol, share=1.0):
        amount = left[symbol] * share
        left[symbol] -= amount
        return amount

    if temperature_k < cond["iron"]:
        take("Fe")
        condensed["iron"] = cond["iron"]
    if temperature_k < cond["forsterite"]:
        magnesium = take("Mg")
        silicon = min(left["Si"], magnesium / 2.0)
        left["Si"] -= silicon
        left["O"] -= OXIDE_OXYGEN["Mg"] * magnesium + OXIDE_OXYGEN["Si"] * silicon
        condensed["forsterite"] = cond["forsterite"]
        # the refractory oxides condense before the silicates
        for symbol in ("Ca", "Ti", "V", "Cr"):
            left["O"] -= OXIDE_OXYGEN[symbol] * take(symbol)
        condensed["oxides"] = cond["forsterite"]
    if temperature_k < cond["enstatite"]:
        silicon = take("Si")
        left["O"] -= OXIDE_OXYGEN["Si"] * silicon
        condensed["enstatite"] = cond["enstatite"]
    if not keep_titanium_vanadium:
        for symbol in ("Ti", "V"):
            left["O"] -= OXIDE_OXYGEN[symbol] * take(symbol)
    if not keep_sodium:
        sodium = take("Na")
        left["S"] -= min(left["S"], sodium / 2.0)
        condensed["sodium_sulfide"] = 0.0
    if not keep_potassium:
        take("K")
        condensed["potassium_chloride"] = 0.0
    return left, condensed


def remove_ammonium_hydrosulfide(budget, temperature_k, pressure_bar, log_enrichment):
    """
    NH4SH forms where `log P + 0.5 [N/H] + 0.5 [S/H] >= (48.91 - 10^4/T) / 4.15` (Visscher et al. 2006, part 05 section 7.2); it
    takes the nitrogen and sulfur out until the gas above the cloud is at that boundary, `(N - x)(S - x) = Kc`.

    Returns:
        tuple: `(budget, x)`: the new budget and the atoms per hydrogen atom of N and S taken (0 when it does not form).
    """
    nitrogen, sulfur = budget["N"], budget["S"]
    if nitrogen <= 0.0 or sulfur <= 0.0:
        return budget, 0.0
    reference = SOLAR_PER_HYDROGEN["N"] * SOLAR_PER_HYDROGEN["S"]
    boundary = (48.91 - 1.0e4 / temperature_k) / 4.15
    k_c = reference * 10.0 ** (2.0 * (boundary - math.log10(pressure_bar)))
    if nitrogen * sulfur <= k_c:
        return budget, 0.0
    x = 0.5 * ((nitrogen + sulfur) - math.sqrt((nitrogen - sulfur) ** 2 + 4.0 * k_c))
    left = dict(budget)
    left["N"] = max(nitrogen - x, 0.0)
    left["S"] = max(sulfur - x, 0.0)
    return left, x


def cap_ices(budget, temperature_k, pressure_bar, equilibrium_at, *, iterations=1000, tolerance=1e-4):
    """
    Caps every species with a saturation curve at its saturation pressure at `temperature_k`: the excess is taken out of the
    budget and the equilibrium solved again, until no species is over its cap by more than `tolerance` (relative).

    Args:
        budget (dict): Atoms per hydrogen atom.
        temperature_k (float): The layer's own temperature (the equilibrium may be solved at a warmer one).
        pressure_bar (float): The layer's pressure.
        equilibrium_at (callable): `budget -> Equilibrium` at the chemistry's temperature and pressure.

    Returns:
        tuple: `(Equilibrium, budget, condensed)` with `condensed` a dict of formula -> moles per hydrogen atom taken out.

    Raises:
        MixtureError: The caps do not settle in `iterations` rounds.
    """
    condensed = {}
    for _ in range(iterations):
        result = equilibrium_at(budget)
        takes = {}
        for name, fraction in result.mole_fractions.items():
            formula = species_thermo.SPECIES[name].formula
            if fraction * pressure_bar < MIN_PARTIAL_BAR or not saturation.has(formula):
                continue
            curve = saturation.saturation(formula, temperature_k)
            if curve.status == saturation.ABOVE_CRITICAL:
                continue
            allowed = curve.pressure_pa / saturation.BAR_TO_PA / pressure_bar
            if fraction > allowed * (1.0 + tolerance):
                takes[name] = (result.moles[name] - allowed * result.total_moles) / (1.0 - min(allowed, 0.999999))
        if not takes:
            return result, budget, condensed
        budget = dict(budget)
        for name, moles in takes.items():
            for symbol, count in species_thermo.SPECIES[name].elements.items():
                budget[symbol] = max(budget.get(symbol, 0.0) - count * moles, 0.0)
            formula = species_thermo.SPECIES[name].formula
            condensed[formula] = condensed.get(formula, 0.0) + moles
        condensed = {formula: moles for formula, moles in condensed.items() if moles > MIN_CONDENSED}
    raise MixtureError(f"the ice caps did not settle at {temperature_k} K and {pressure_bar} bar")


# --- the quench -------------------------------------------------------------------------------------------------------------

def log_carbon_monoxide_time(temperature_k, pressure_bar, metallicity):
    """ln of the chemical time of the CO to CH4 conversion in seconds (Zahnle and Marley 2014: `1/t = 1/t_q1 + 1/t_q2`), in logs so
    that a cold layer does not overflow."""
    log_t1 = math.log(1.5e-6) - math.log(pressure_bar) - 0.7 * math.log(metallicity) + 42000.0 / temperature_k
    log_t2 = math.log(40.0) - 2.0 * math.log(pressure_bar) + 25000.0 / temperature_k
    return -float(np.logaddexp(-log_t1, -log_t2))


def mixing_time_s(temperature_k, mean_molar_mass, gravity, kzz_cm2_s):
    """The time to mix one scale height, s: `H^2 / Kzz` with `H = R T / (M g)`."""
    scale_height = constants.R * temperature_k / (mean_molar_mass * gravity)
    return scale_height ** 2 / (kzz_cm2_s * 1e-4)


def quench_level(temperature_k, pressure_bar, metallicity, mean_molar_mass, gravity, kzz_cm2_s):
    """
    The `(T, P)` where the carbon monoxide chemistry freezes, on an adiabat down from the layer: the shallowest level below
    the layer where the chemical time is no longer longer than the mixing time. `None` when the layer itself is in
    equilibrium (the chemical time is already shorter than the mixing time there).
    """
    def excess(pressure):
        temperature = temperature_k * (pressure / pressure_bar) ** ADIABAT_EXPONENT
        return (log_carbon_monoxide_time(temperature, pressure, metallicity)
                - math.log(mixing_time_s(temperature, mean_molar_mass, gravity, kzz_cm2_s))), temperature

    if excess(pressure_bar)[0] <= 0.0:
        return None
    low = pressure_bar
    steps = 64
    ratio = (QUENCH_DEEPEST_BAR / pressure_bar) ** (1.0 / steps)
    high = None
    for step in range(1, steps + 1):
        trial = pressure_bar * ratio ** step
        if excess(trial)[0] <= 0.0:
            high = trial
            break
        low = trial
    if high is None:
        return excess(QUENCH_DEEPEST_BAR)[1], QUENCH_DEEPEST_BAR
    for _ in range(60):
        mid = math.sqrt(low * high)
        if excess(mid)[0] > 0.0:
            low = mid
        else:
            high = mid
    pressure = math.sqrt(low * high)
    return excess(pressure)[1], pressure


# --- the family ------------------------------------------------------------------------------------------------------------

def metal_mass_fraction(log_enrichment):
    """The mass fraction of the metals of `element_budget(log_enrichment)`: used to turn an envelope's drawn Z into a metallicity."""
    from planetgen.physics import elements
    budget = element_budget(log_enrichment)
    mass = {symbol: count * elements.element(symbol).atomic_weight for symbol, count in budget.items()}
    metals = sum(mass[symbol] for symbol in METALS)
    return metals / (metals + mass["H"] + mass["He"])


def enrichment_of_mass_fraction(z):
    """log10 of the metals over the reference that gives the metal mass fraction `z` (bisection; Z rises with the enrichment)."""
    low, high = -3.0, 4.0
    for _ in range(60):
        mid = 0.5 * (low + high)
        if metal_mass_fraction(mid) < z:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def _layer_temperature(subject, family):
    for name in {"hot-jupiter": ("t_layer_k", "tday_k", "teq_k"), "brown-dwarf": ("t_layer_k", "teff_k")}.get(family, ("t_layer_k", "teq_k")):
        if subject.get(name) is not None:
            return float(subject[name])
    raise MixtureError(f"the subject of class {subject.body_class!r} has no layer temperature")


def _gravity(subject):
    mass = subject.get("mass_jupiter")
    mass_kg = _mjup_to_kg(mass) if mass is not None else _mearth_to_kg(subject["mass_earth"])
    return constants.G * mass_kg / (subject["radius_km"] * 1e3) ** 2


def family_of(recipe):
    """Which giant regime the recipe drew: `brown-dwarf` (`g`), `hot-jupiter` (`f4` with a `h` row) or `giant`."""
    if "g" in recipe.groups:
        return "brown-dwarf"
    if "h" in recipe.groups:
        return "hot-jupiter"
    return "giant"


def _flag(recipe, row, default):
    return recipe.binaries[row][0] if row in recipe.binaries else default


def _alkali_in_gas(subject, row, temperature_k, onset_k):
    """Whether the alkali metal is in the gas of the layer: the seeded draw of the recipe row `m7` (Na) or `m8` (K) where the row
    applies to the subject, else the temperature test at the row's centre."""
    entry = RECIPES.entry(row)
    if entry.supports(subject):
        try:
            return RECIPES.draw_binary(row, subject)[0]
        except KeyError as error:
            raise MixtureError(str(error)) from error
    return temperature_k >= onset_k


def giant_mixture(subject, recipe, *, log_enrichment=None, family=None):
    """
    The mixture of a giant, brown dwarf or hot Jupiter for `subject` and its recipe draw (see the module text).

    Args:
        subject (Subject): With `radius_km`, `mass_jupiter` (or `mass_earth`), a layer temperature (`t_layer_k`, else `tday_k` for a hot
            Jupiter, `teff_k` for a brown dwarf, `teq_k`) and optionally `fe_h`, the host's [Fe/H] (0 where the host's is not known).
        recipe (RecipeDraw): The draw of `CATEGORIES`.
        log_enrichment (float, optional): log10 of the metals over the reference abundance, instead of the recipe's `f5` over the
            host's (the envelope family gives it from its drawn Z). The result's `constants["log_enrichment"]` is this number.
        family (str, optional): Overrides `family_of(recipe)`.

    Returns:
        MixtureDraw: The mixture, `redox_diw` None.
    """
    require(subject, ("radius_km",))
    family = family or family_of(recipe)
    constants_drawn = {}
    temperature = _layer_temperature(subject, family)
    pressure = LAYER_PRESSURE_BAR[family]
    fe_h = float(subject.get("fe_h") or 0.0)
    if log_enrichment is None:
        log_enrichment = math.log10(recipe.values.get("f5", 1.0)) + fe_h
    metallicity = 10.0 ** log_enrichment
    constants_drawn["log_enrichment"] = log_enrichment
    c_over_o = recipe.values.get("f7")
    if c_over_o is not None:
        constants_drawn["c_over_o"] = c_over_o
    budget = element_budget(log_enrichment, c_over_o)

    budget, condensed = remove_clouds(
        budget, temperature, pressure, log_enrichment,
        keep_titanium_vanadium=_flag(recipe, "h4", True),
        keep_sodium=_alkali_in_gas(subject, "m7", temperature, ALKALI_ONSET_K["Na"]),
        keep_potassium=_alkali_in_gas(subject, "m8", temperature, ALKALI_ONSET_K["K"]))
    budget, taken = remove_ammonium_hydrosulfide(budget, temperature, pressure, log_enrichment)
    if taken:
        condensed["ammonium_hydrosulfide"] = taken
    chemistry_k = max(temperature, CHEMISTRY_FLOOR_K)
    everything = tuple(species_thermo.GASES)
    without_co = tuple(name for name in everything if name != "CO")

    first = gas_equilibrium.solve(budget, chemistry_k, pressure)
    mean_molar_mass = sum(fraction * species_thermo.SPECIES[name].molar_mass for name, fraction in first.mole_fractions.items())
    kzz = bracket(constants_drawn, "kzz_cm2_s", KZZ_CM2_S[family])
    quench = quench_level(temperature, pressure, metallicity, mean_molar_mass, _gravity(subject), kzz)
    frozen = 0.0
    species = everything
    layer_budget = budget
    if quench is not None:
        quench_k, quench_bar = quench
        deep = gas_equilibrium.solve(budget, max(quench_k, CHEMISTRY_FLOOR_K), quench_bar).moles
        carriers = deep.get("CO", 0.0) + deep.get("CH4", 0.0)
        frozen = deep.get("CO", 0.0) / carriers if carriers > 0.0 else 0.0
        constants_drawn.update(quench_k=quench_k, quench_bar=quench_bar)
        carbon_monoxide = min(frozen * budget["C"], budget["O"])
        layer_budget = dict(budget, C=budget["C"] - carbon_monoxide, O=budget["O"] - carbon_monoxide)
        species = without_co
    else:
        carbon_monoxide = 0.0
    constants_drawn["carbon_monoxide_frozen"] = frozen

    def equilibrium_at(trial):
        return gas_equilibrium.solve(trial, chemistry_k, pressure, species=species)

    result, layer_budget, ices = cap_ices(layer_budget, temperature, pressure, equilibrium_at)
    for formula, moles in ices.items():
        condensed[formula] = moles
    if carbon_monoxide > 0.0:
        moles = dict(result.moles)
        moles["CO"] = moles.get("CO", 0.0) + carbon_monoxide
        total = sum(moles.values())
        result = gas_equilibrium.Equilibrium(
            temperature_k=result.temperature_k, pressure_bar=result.pressure_bar, moles=moles,
            mole_fractions={name: value / total for name, value in moles.items()}, total_moles=total,
            iterations=result.iterations, skipped=result.skipped)
    result = gas_equilibrium.Equilibrium(
        temperature_k=temperature, pressure_bar=pressure, moles=result.moles, mole_fractions=result.mole_fractions,
        total_moles=result.total_moles, iterations=result.iterations, skipped=result.skipped)
    return finish(result, family=family, drawn={group: drawn.name for group, drawn in recipe.groups.items()}, redox_diw=None,
                  inventory=element_budget(log_enrichment, c_over_o), constants=constants_drawn, condensed=condensed)
