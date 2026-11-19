# planetgen/physics/carbonate_buffer.py

"""
The dry-surface carbonate buffer: the CO2 a hot rock surface holds against the air once the ocean has gone (GEN.1767)
====================================================================================================================

`surface_air._weathers` switches the weathering of carbon off on a body whose water has all gone to vapor or whose surface is over the
critical temperature of water. The carbon its volcanoes brought up then stayed in the air whole, at pressures the rock does not allow:
a hot surface of calcium silicate and silica holds the CO2 of the air at the pressure of the Urey reaction

    CaCO3 (calcite) + SiO2 (quartz) = CaSiO3 (wollastonite) + CO2,

the buffer that Schaefer and Fegley take to hold the CO2 of Venus at 92.1 bar and 740 K ([R5313]: "Our model explicitly assumes that this
buffer (C1) controls the CO2 pressure in the lower Venusian atmosphere"). Above the pressure of the reaction the rock carbonates and takes
CO2 out of the air; below it the carbonate breaks down and gives it back.

* *The pressure of the reaction.* `equilibrium_bar(T)` is the equilibrium constant of the reaction from the standard-state data of
  [R1981: Table 2a] (enthalpy of formation, entropy and heat capacity `a + b T + c / T^2 + d / T^(1/2)` of calcite, quartz, wollastonite
  and CO2), the CO2 an ideal gas and the three solids pure at one bar. Against the sources read: 88 bar at 740 K (Venus: 92.1 bar, [R5313]); 257 bar
  at 800 K (the experiments of [R5314], read through the pressure and enthalpy a later paper prints, reach pure CO2 at about 275 bar near 800 K). The pressure term of the three solids (their volume
  falls by 1.97 cm3 per mol of reaction) raises the pressure by 3 percent at 100 bar and is left out.
* *The limit.* The carbonate is made of the crust's calcium: `capacity_bar` is the CO2 pressure that all the calcium of the crust column that
  the volcanism built (`sulfur_balance.CRUST_CALCIUM_WT`), less the calcium the sulfate sink took, would hold. A body with more carbon than that
  keeps the rest in the air.
* *Where it applies.* To a body with a surface whose carbon is not weathered (`surface_air._shares`: not a giant, weathering off: a body that never had an
  ocean, or whose water is all vapor, or whose surface is over the critical temperature of water) and whose surface is at least `DRY_CARBONATION_MIN_K`. Venus is the
  case the equilibrium reading is made for ([R5313] take it dry); [R5315] call the CO2-silicate reactions of Venus kinetically inhibited, and the newer source wins. Below that
  temperature a dry rock carbonates too slowly to count ([R5316], citing the earlier measurements it names: without liquid water the CO2 "may not be controlled by silicate-carbonate
  equilibria because of the inefficient carbonation of dry near-surface rocks"), so a cool dry body keeps its air as it is. No upper limit is needed: above about 800 K the pressure of the reaction is over
  what any body holds, and the buffer takes nothing.

The kept pressure is `min(supply, max(equilibrium, supply - capacity))`: the air holds the equilibrium pressure while the crust has calcium to
carbonate and the whole surplus when it does not.

Every number that is not read from a source is a design constant with its tier and tweak level in `docs/design/research/gap-filling-plan.md`
(GEN.725 rule, section 6.310).
"""

import math
from functools import lru_cache

from planetgen.physics import elements, gases, sulfur_balance

GAS_CONSTANT_KJ = 8.314462618e-3
"""float: The gas constant, kJ per mole per kelvin ([R4361], measured)."""

REFERENCE_K = 298.15
"""float: The temperature of the standard state of the data set, K [R1981]."""

# Table 2a of [R1981] ([R1981], J. Metamorphic Geol. 29, 333): enthalpy of formation (kJ/mol), entropy (J/K/mol) and the heat capacity a, b, c, d of
# Cp = a + b T + c T^-2 + d T^-1/2 in kJ/K/mol (b read as the 1e-5 kJ/K^2 scale that reproduces the 298 K heat capacities to 3 percent: [Stand-in 6.310 CB4, tier T6 educated guess from physics, tweak L0]).
PHASES = {
    "calcite": (-1207.88, 92.50, 0.1409, 0.5029e-5, -950.7, -0.8584),
    "wollastonite": (-1633.75, 82.50, 0.1593, 0.0, -967.3, -1.0754),
    "quartz": (-910.70, 41.43, 0.0929, -0.0642e-5, -714.9, -0.7161),
    "carbon dioxide": (-393.51, 213.70, 0.0878, -0.2644e-5, 706.4, -0.9989),
}
"""dict: The four phases of the reaction: `(enthalpy of formation, entropy, a, b, c, d)` ([R1981: Table 2a]). [Stand-in 6.310 CB1, tier T1 verified, tweak L0.]"""

DRY_CARBONATION_MIN_K = 500.0
"""float: The coolest surface at which a dry rock buffers the CO2 of its air, K. [Stand-in 6.310 CB2, tier T6 educated guess from physics, tweak L3: below it dry rock
carbonates too slowly to count [R5316]; the pressure of the reaction is 0.08 bar there, so a body just over it holds under a tenth of a bar.]"""

CALCIUM_SHARE = 1.0
"""float: The share of the crust's calcium (less the sulfate's) that can carbonate. [Stand-in 6.310 CB3, tier T6 educated guess from physics, tweak L3: the item sets the limit at the crust's
calcium; the Earth's own crust holds about 17 percent of its calcium as carbonate, so a lower share is the first thing to tune.]"""


def _enthalpy(phase, temperature_k):
    """kJ/mol: the enthalpy of formation of `phase` at `temperature_k`."""
    h, _, a, b, c, d = PHASES[phase]
    t0 = REFERENCE_K
    return h + a * (temperature_k - t0) + 0.5 * b * (temperature_k ** 2 - t0 ** 2) - c * (1.0 / temperature_k - 1.0 / t0) + 2.0 * d * (math.sqrt(temperature_k) - math.sqrt(t0))


def _entropy(phase, temperature_k):
    """kJ/K/mol: the entropy of `phase` at `temperature_k`."""
    _, s, a, b, c, d = PHASES[phase]
    t0 = REFERENCE_K
    return (s / 1000.0 + a * math.log(temperature_k / t0) + b * (temperature_k - t0) - 0.5 * c * (1.0 / temperature_k ** 2 - 1.0 / t0 ** 2)
            - 2.0 * d * (1.0 / math.sqrt(temperature_k) - 1.0 / math.sqrt(t0)))


@lru_cache(maxsize=4096)
def equilibrium_bar(temperature_k):
    """The CO2 pressure of the Urey reaction at `temperature_k`, bar: `exp(-dG / R T)` with `dG` the Gibbs energy of wollastonite and CO2 less that of calcite and quartz; the
    CO2 an ideal gas, the solids pure at one bar. 0.08 bar at 500 K, 3 bar at 600 K, 88 bar at 740 K, 1,100 bar at 900 K."""
    delta_g = sum(sign * (_enthalpy(phase, temperature_k) - temperature_k * _entropy(phase, temperature_k))
                  for phase, sign in (("wollastonite", 1.0), ("carbon dioxide", 1.0), ("calcite", -1.0), ("quartz", -1.0)))
    exponent = -delta_g / (GAS_CONSTANT_KJ * temperature_k)
    return math.exp(min(exponent, 700.0))


def capacity_bar(crust_kg_m2, sulfate_mol_m2, gravity_m_s2):
    """The CO2 pressure, bar, the crust's calcium can carbonate: the calcium of `crust_kg_m2` of crust (`sulfur_balance.CRUST_CALCIUM_WT`) less the calcium the sulfate sink took
    (`sulfate_mol_m2`, one per sulfur), times `CALCIUM_SHARE`, one CO2 per calcium, over the surface."""
    calcium_mol_m2 = crust_kg_m2 * sulfur_balance.CRUST_CALCIUM_WT / elements.molar_mass("Ca") - sulfate_mol_m2
    return max(calcium_mol_m2, 0.0) * CALCIUM_SHARE * gases.gas("CO2").molar_mass * gravity_m_s2 / 1.0e5


def kept_bar(supply_bar, temperature_k, capacity_bar_):
    """The CO2 pressure of the air of a dry body at `temperature_k`, bar, when its volcanoes delivered `supply_bar` and its crust can carbonate `capacity_bar_`: the pressure of the
    reaction while the calcium lasts, the surplus over the calcium after, never over the supply, and the whole supply under `DRY_CARBONATION_MIN_K`."""
    if temperature_k < DRY_CARBONATION_MIN_K:
        return supply_bar
    return min(supply_bar, max(equilibrium_bar(temperature_k), supply_bar - capacity_bar_))
