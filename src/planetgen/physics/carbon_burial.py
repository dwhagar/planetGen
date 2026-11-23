# planetgen/physics/carbon_burial.py

"""
The carbon a biosphere buries, and what it does to the CO2 of its world (GEN.1821)
==================================================================================

The Boss chose, at 3:25 PM EST 11-22, the honest way to say that life draws CO2 down (report `reports/co2-carbon-sink-gen1821-2026-11-21.md`, census `reports/co2-sink-census-gen1821-2026-11-21/`): the weathering law is
levelled so that the Earth balances **with** its biosphere (`surface_pressure.EARTH_NET_BURIAL_SHARE`, 3 Tmol of carbon a year net over 6.6 of volcanic outgassing), a world with no life weathers 1 / (1 - 0.455) times the
modern rate to balance its volcanoes, and a world with life gets the CO2 that its own burial leaves. This module is the burial, `b`, in units of Earth's outgassing per area:

* **Oxygenic classes** (those that split water). The carbon they bury is the source `S` of the oxygen law (`biogenic_oxygen.source_of`); what the oxygen takes back is the oxidative weathering `W(p)` of the kerogen
  (GEN.1733), so the net sink is `S - W(p)`: all of `S` while no oxygen stands, `R` at the balance of the oxygen. [Stand-in 6.322 CB3, tier T1 (the oxygen law's own numbers), tweak L0.]
* **Classes 2 to 4**, which make no oxygen (the burial rows of `life_ledger.json`). The carbon fixed is the electrons the world's donors or hydrogen give (`life_supply`, a share of Earth's, times the area), a quarter of them
  to carbon at the ceiling (the light-fed Class 2, the figure of [R3393]: 80 Tmol of Fe2+ a year fix 20 Tmol of carbon) or at the growth yield of the chemical classes (`life_supply.json` `sealed_fuel_efficiency`, GEN.1825, [R5348]);
  of that `anoxic_burial_share` (2 percent, 10 times the modern ocean's, the anoxic bottom water of the Archean) is buried and the rest returns as carbon dioxide and methane. Classes that eat from one fuel
  are counted once (the larger), as the methane of the ledger is. Classes 6 to 9 have no row (the carbonate-silicate cycle is not the sink of a dense carbon dioxide fluid, an ammonia sea, an acid or a methane sea).
  [Stand-in 6.322 CB4, tier T5 (the share) and T3 (the electrons), tweak L1.]
* **The floor.** The burial falls in a line to nothing as the CO2 falls to the 0.5 Pa floor of the producers (`surface_pressure.burial_rate`), so a biosphere that buries more than the volcanoes give drains the air to the
  floor and not under it. [Stand-in 6.322 CB2, tier T7, tweak L1.]

A body with no carbon sink (no carbon to weather, a sealed ocean's chemistry, a dry carbonate buffer), a giant, or no life takes no burial. The burial is worked out once from the air the life has made: its temperature moves the
oxygenic source a little (`biogenic_oxygen.temperature_factor`), and the CO2 that results moves the temperature in turn; one pass is the law here (the census of the item measured the same pass).
"""

from planetgen.physics import atmosphere, biogenic_oxygen, habitability_world, life_air, life_supply

OUTGASSING_EARTH_MOL_YR = 6.6e12
"""float: Earth's modern volcanic outgassing of CO2, mol of carbon a year (6.6 Tmol, [R2489], registry AM60; `waterworld_carbon.SEAFLOOR_SHARE_OF_OUTGASSING`'s reference)."""

ELECTRONS_PER_CARBON = 4.0
"""float: The electrons that reduce one carbon dioxide to carbohydrate (derived: the stoichiometry of CO2 + 4 H+ + 4 e- to CH2O + H2O), so a share `yield` of a mole of oxygen-equivalent of electrons (4) is a mole of carbon at the ceiling."""

TMOL = life_air.TMOL


def oxygenic_net_mol_yr(body, living):
    """The net carbon, mol a year, that the oxygenic classes of `living` bury on `body`: `S - min(W(p), S)` (see the module's docstring); 0 for a body with none."""
    pairs = life_air._oxygenic(body, living, forever=False)
    if not pairs:
        return 0.0
    area = biogenic_oxygen.area_m2(body)
    source = biogenic_oxygen.source_of(body, area, [(share, band) for _, share, band in pairs])
    pascal = atmosphere.gas_partials_pa(body).get("O2", 0.0)
    weathering = biogenic_oxygen.EARTH_SURPLUS_MOL_YR * (area / biogenic_oxygen.EARTH_AREA_M2) * (pascal / biogenic_oxygen.EARTH_O2_PA) ** biogenic_oxygen.OXIDATIVE_WEATHERING_EXPONENT
    return max(source - min(weathering, source), 0.0)


def fixed_carbon_mol_yr(body, row, supplies):
    """The carbon, mol a year, the class of ledger row `row` fixes on `body` (0 for a row with no burial): the electrons of its fuel (Earth's, times the supply's share of Earth's and the area's) over four, at its yield."""
    burial = row.get("burial")
    if row["status"] != "built" or not burial:
        return 0.0
    area_ratio = biogenic_oxygen.area_m2(body) / biogenic_oxygen.EARTH_AREA_M2
    electrons = life_air.number(burial["electrons_tmol_yr"]) * burial["electrons_per_unit"] * TMOL * supplies.get(burial["fuel"], 0.0) * area_ratio
    share = 1.0 if burial["yield"] == "ceiling" else life_supply.number("sealed_fuel_efficiency")
    return electrons / ELECTRONS_PER_CARBON * share


def anoxic_net_mol_yr(body, living, supplies=None):
    """The carbon, mol a year, that the classes of `living` with a burial row bury on `body`: each fuel's largest fixed carbon (classes on one fuel are not added), times `anoxic_burial_share`."""
    if supplies is None:
        supplies = life_supply.supplies(habitability_world.values_from_body(body))
    star_age = getattr(getattr(body, "star", None), "age", None)
    by_fuel = {}
    for one in living:
        started = one.started(body)
        row = life_air.row_of(one.life_class)
        if not row.get("burial") or started is None or star_age is None or started > star_age:
            continue
        fuel = row["burial"]["fuel"]
        by_fuel[fuel] = max(by_fuel.get(fuel, 0.0), fixed_carbon_mol_yr(body, row, supplies))
    return life_air.number("anoxic_burial_share") * sum(by_fuel.values())


def rate_of(body, living):
    """`b`, the net organic sink of the life on `body` over Earth's outgassing per area (0 for a giant, a body with no radius or no life)."""
    if getattr(body, "body_type", "t") == "g" or not getattr(body, "radius", None) or not living:
        return 0.0
    burial = oxygenic_net_mol_yr(body, living) + anoxic_net_mol_yr(body, living)
    area = biogenic_oxygen.area_m2(body)
    return (burial / area) / (OUTGASSING_EARTH_MOL_YR / biogenic_oxygen.EARTH_AREA_M2)


def apply(body, living):
    """Solves the carbon sink of `body` again with the burial of its life (`biogenic_oxygen.rebalance_co2`): the CO2 its own biosphere leaves. Returns whether the air changed."""
    burial = rate_of(body, living)
    return burial > 0.0 and biogenic_oxygen.rebalance_co2(body, burial)
