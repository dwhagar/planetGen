# planetgen/physics/ocean_oxygen.py

"""
The oxygen of a sealed ocean (GEN.1824; stand-ins section 6.324)
================================================================

An ocean under an ice lid (hydrosphere `ice-covered ocean`) has no air to breathe from. Its dissolved oxygen is a stock that builds while oxidant arrives faster than the floor's reductants take it, as in the air law of
`physics/biogenic_oxygen.py` (a source `S`, a reductant flux `R`, nothing held while `S` is not above `R`) but read the way [R5346] read it for Europa: the oxidant is made in the surface ice all the time and reaches the ocean in pulses,
once for each resurfacing of the lid, so between two pulses the ocean holds what the surplus `S - R` has put into its water, and the reductants, not the oxidant, are the limit. All in mol of O2 per m2 of the ocean per year:

* **The lid** (`lid_delivery`): oxidant made in the surface ice by the particles that strike it (the radiolysis of ice), `SURFACE_O2_MOL_M2_YR` at the dose of Europa's surface `REFERENCE_DOSE_MSV_YR`, going with the energy the ice takes in, so with the body's stored
  surface dose (`surface_dose_msv_yr`: cosmic rays, solar particles and, for a moon, the belts of its giant); `LID_DELIVERY_SHARE` of it reaches the ocean [R5336]; at most `LID_DELIVERY_MAX_MOL_M2_YR` (the top of the published range, [R5347]).
* **The water** (`radiolysis`): the potassium-40, uranium and thorium in the water split it, 0.2 oxygen for each 100 eV ([R2691] for the hydrogen yield, the oxygen half of it by redox balance), `RADIOLYSIS_O2_MOL_KG_YR` for each kilogram of
  water with Earth's seawater potassium ([R5337] for an icy moon's). It goes with the depth of the ocean, to `RADIOLYSIS_DEPTH_MAX_KM`.
* **The faults** (`fracture`): the radicals of fresh silicate surfaces broken in water [R5338], `FRACTURE_MAX_MOL_M2_YR` at the heat flow of the Earth's, less where there is less heat to crack the floor.
* **The sink** (`sink`): the reductants of the floor, the hydrogen of the water's reaction with the rock and the volcanic reduced gases and iron, which are two channels of one pool of hydrothermal reductants (Earth: 1e12 mol of H2 a year at its ridges, [R5339]):
  `H2_O2_EQUIVALENT_MOL_M2_YR` for each unit of the mean of `life_supply`'s `h2_supply` and `donor_supply`, half an oxygen for each hydrogen.
  The water's radiolysis makes hydrogen with its oxygen, in equal electron equivalents, so its hydrogen is in the sink (`RADIOLYSIS_H2_ESCAPE` of it leaves): the radiolysis is a net oxidant only for the hydrogen that goes.
* **The bank** (`held_um`, `held_kpa`): the surplus times the interval between pulses `PULSE_INTERVAL_YR` (30 to 70 Myr, the resurfacing of Europa's lid, [R5346]) over the water the ocean holds, as micromolar of dissolved oxygen, at most the pressure the air law lets the air hold
  (`mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA`), read as the partial pressure `HENRY_UM_PER_KPA` gives (0.1 kPa over water is about 2 micromolar). With a surplus below zero the ocean holds none. A share `INHERITED_O2_RETAINED` of the oxygen
  the ice accreted with (a comet's 3.8 percent of its water) is a second, older stock, drained by the reductants the surplus does not cover over the age of the body; it is zero until a source constrains it.

`ocean_o2_kpa` is the partial pressure that stock holds; `None` for a body that is not a sealed ocean. `life_complexity.assess` reads it as the oxygen of the animals of every class that lives there: a chemical class breathes the ocean's oxygen for its grade and not
the 100 K surface's air, and a class that breathes oxygen reads the larger of the air's and the ocean's.

What it leaves out: the oxygen that an ice-covered ocean with light under thin ice gets from photosynthesis (the oxygenic classes, `biogenic_oxygen`, whose oxygen goes to the air); the nodule "dark oxygen" of the research report (disputed; not built); the mixing of the ocean (one
well-mixed value, taken at the end of a pulse interval); the delivery of oxidant by subduction and clathrate (they change when it arrives, not how much); the fuel cap on the production of life (a separate card); the oxygen that an air above the lid gives (the larger of the two is read, not the sum).

Every number that is not read from a source is a design constant with its tier and tweak level in `docs/design/research/gap-filling-plan.md` (section 6.324), and the research is
`reports/dark-oxygen-production-2026-11-21.md` and `reports/dark-ocean-untested-theories-2026-11-21.md` (section 1).
"""

import dataclasses

from planetgen.physics import life_supply, mixture_oxygen

SEALED = ("ice-covered ocean",)
"""tuple: The hydrospheres whose ocean has no air and no light of its own: the ocean lies under an ice lid."""

SURFACE_O2_MOL_M2_YR = 3.9e-4
"""float: The oxygen the particles make in the surface ice of Europa, mol of O2 per m2 of its surface per year: 12 kg a second (the Juno measurement of Szalay and others 2024, second-hand through the research report) is 1.2e10 mol a year over its area of
3.1e13 m2. [Stand-in 6.324 SO1, tier T4, tweak L0.]"""

REFERENCE_DOSE_MSV_YR = 1.97e6
"""float: The dose of the surface of Europa, mSv a year (5.4 Sv a day in the belts of Jupiter, recalled), the dose at which `SURFACE_O2_MOL_M2_YR` is made. The production of oxidant in the ice goes with the energy the ice takes in (a constant yield for each eV),
so it is taken to go with the stored surface dose. [Stand-in 6.324 SO2, tier T5, tweak L1.]"""

LID_DELIVERY_SHARE = 0.85
"""float: The share of the surface oxidant that reaches the ocean through the lid (the model of Hesse and others 2022, [R5336], whose chaos-terrain limit on the rate is the lower end of the delivery the pulses of the bank stand for). Read for Europa; taken for every lid. [Stand-in 6.324 SO3, tier T3, tweak L0.]"""

LID_DELIVERY_MAX_MOL_M2_YR = 9.7e-3
"""float: The most the lid delivers, mol of O2 per m2 per year: the top of the published range for Europa, 3e11 mol a year ([R5347]; the abstract was read, the first research report quoted 1e11 second-hand), over its area of 3.1e13 m2. [Stand-in 6.324 SO4, tier T3, tweak L0.]"""

RADIOLYSIS_O2_MOL_KG_YR = 8.0e-13
"""float: The oxygen the radioactive decay in seawater makes, mol of O2 for each kg of water a year: Earth's seawater potassium (0.388 g a kg, 40K 1.17e-4 of it, 12 Bq a kg, 0.66 MeV a decay, 1.3e-12 W a kg) at 0.2 O2 for each 100 eV, half the hydrogen
yield of 0.4 for each 100 eV of gamma and beta [R2691]. The icy moons' potassium is of the same order [R5337]. [Stand-in 6.324 SO5, tier T5, tweak L0.]"""

RADIOLYSIS_DEPTH_MAX_KM = 100.0
"""float: The depth of ocean that the radiolysis counts, km: a column of 100 km gives 8e-5 mol per m2 per year; a deeper ocean has ice on its floor and its water below the depth is out of the circulation. [Stand-in 6.324 SO6, tier T6, tweak L1.]"""

FRACTURE_MAX_MOL_M2_YR = 1.0e-5
"""float: The oxidant of the radicals of a sealed seafloor's faults, mol of O2 per m2 per year, at the heat flow of Earth's: 1e-6 to 1e-5, the top of the range, which is the glacial abrasion of [R5338] with a tenth of the
sediment (no rivers or shallows in a sealed ocean) and a fault's fresh surface for the rest. It goes with the heat flow the body has over Earth's, to 1. [Stand-in 6.324 SO7, tier T6, tweak L1.]"""

H2_O2_EQUIVALENT_MOL_M2_YR = 0.5 * 1.0e12 / 3.61e14
"""float: The oxygen Earth's hydrothermal reductants take, mol of O2 for each m2 of the floor a year: 1e12 mol of H2 a year at the ridges (Worman and others 2016, [R5339]; the abstract was read, the ranges in the full text were not) over Earth's
ocean of 3.61e14 m2, half an oxygen for each hydrogen. The sink of a body is this times the mean of its `h2_supply` and `donor_supply` (Earth's are 1). Cassini's Enceladus (1e9 to 5e9 mol of H2 a year) is 1 to 4 times it per area, the estimate of Europa (1.6e10, second-hand) 0.2 times. [Stand-in 6.324 SO8, tier T3, tweak L0.]"""

PULSE_INTERVAL_YR = 50.0e6
"""float: The years between two deliveries of the surface oxidant to the ocean, the resurfacing of the lid: 30 to 70 Myr for Europa, the middle ([R5346]: at that interval the ocean holds the oxygen of terrestrial surface waters even with the hydrothermal reductants taking most of it). Taken for every lid. [Stand-in 6.324 SO9, tier T3, tweak L1.]"""

HENRY_UM_PER_KPA = 22.0
"""float: The dissolved oxygen, micromolar, in water at 0 C over a partial pressure of 1 kPa of oxygen (2.2 millimolar over 1 bar: 0.1 kPa is about 2 micromolar and 6 kPa about 130). [Stand-in 6.324 SO10, tier T1, tweak L0.]"""

WATER_MOL_PER_L = 55.5
"""float: The moles of water in a litre (derived: 1,000 g over 18.015 g a mole)."""

COMET_O2_PER_WATER = 0.038
"""float: The oxygen in the ice of comet 67P, mol of O2 for each mol of water (3.80 plus or minus 0.85 percent, the Rosetta measurement; two papers argue it is primordial, trapped in the grains). Only a share of it is kept. [Stand-in 6.324 SO11, tier T1, tweak L0.]"""

INHERITED_O2_RETAINED = 0.0
"""float: The share of the oxygen the ice accreted with (`COMET_O2_PER_WATER`) that an ocean keeps through the heating, hydration and oxidation of the rock of its first Gyr. No source constrains it: the research report suggests 0.01 (range 0.0001 to 0.1), enough to keep an
oxygen-rich ocean alive against 4.5 Gyr of reductants, and it is off here (0) until Boss turns it on. It would put every sealed ocean over the oxygen floors. [Stand-in 6.324 SO12, tier T4, tweak L2.]"""

RADIOLYSIS_H2_ESCAPE = 0.0
"""float: The share of the hydrogen the radiolysis makes with its oxygen that leaves the ocean (into the rock or the ice). The water splits into both in equal electron equivalents (0.4 H2 and 0.2 O2 for each 100 eV, `RADIOLYSIS_O2_MOL_KG_YR`), so
the oxygen is a net oxidant only for the part of the hydrogen that goes: in a sealed ocean the rest takes the oxygen back, as a reductant (the research report, section 3.2: "unless it leaves, the water does not become net oxidizing"). No source gives the share:
0 keeps the water's radiolysis net zero (the oxygen counted and its hydrogen counted in the sink, `sink`). [Stand-in 6.324 SO13, tier T6, tweak L1.]"""

BODY_AGE_YR = 4.5e9
"""float: The age the inherited stock is drained over, years: the Solar System's 4.5 Gyr (a body's own age is not stored where the grade is read, and the grade must read the same numbers in the build and in a rescore). [Stand-in 6.324 SO14, tier T6, tweak L1.]"""


@dataclasses.dataclass(frozen=True)
class Balance:
    """The oxygen of a sealed ocean: `lid`, `water` and `faults` the sources and `sink` the reductants, mol of O2 per m2 per year; `um` the dissolved oxygen that holds, micromolar, and `kpa` the partial pressure it stands for."""
    lid: float
    water: float
    faults: float
    sink: float
    um: float
    kpa: float

    @property
    def source(self):
        return self.lid + self.water + self.faults


def is_sealed(values):
    """Whether the body's ocean lies under an ice lid with no air of its own."""
    return values.get("hydrosphere") in SEALED


def lid_delivery(values):
    """The oxidant the ice lid delivers to the ocean, mol of O2 per m2 per year, on average over a pulse interval: Europa's surface production at its dose, to the dose of the body, the share that gets through, at most `LID_DELIVERY_MAX_MOL_M2_YR`."""
    dose = max(0.0, values.get("surface_dose_msv_yr") or 0.0)
    return min(LID_DELIVERY_MAX_MOL_M2_YR, LID_DELIVERY_SHARE * SURFACE_O2_MOL_M2_YR * dose / REFERENCE_DOSE_MSV_YR)


def radiolysis(values):
    """The oxidant the decay in the ocean's water makes, mol of O2 per m2 of the ocean per year: the water over a square metre (1,000 kg a m3 of the depth, at most `RADIOLYSIS_DEPTH_MAX_KM`) times the rate for each kg."""
    depth_km = min(RADIOLYSIS_DEPTH_MAX_KM, max(0.0, values.get("ocean_depth_km") or 0.0))
    return RADIOLYSIS_O2_MOL_KG_YR * depth_km * 1000.0 * 1000.0


def fracture(values):
    """The oxidant of the faults' radicals, mol of O2 per m2 per year: `FRACTURE_MAX_MOL_M2_YR` at Earth's heat flow, less with less heat (to 1)."""
    return FRACTURE_MAX_MOL_M2_YR * min(1.0, life_supply.heat_share(values))


def sink(values, supplies=None):
    """The reductants that take the oxygen, mol of O2-equivalent per m2 per year: Earth's hydrothermal flux times the mean of the body's `h2_supply` and `donor_supply` (`life_supply`; a share stored in `values` replaces the one worked out, as it does for the classes' scores),
    and the hydrogen of the water's own radiolysis that does not leave the ocean (`RADIOLYSIS_H2_ESCAPE`)."""
    shares = supplies if supplies is not None else {**life_supply.supplies(values), **{name: values[name] for name in life_supply.QUANTITIES if values.get(name) is not None}}
    return H2_O2_EQUIVALENT_MOL_M2_YR * 0.5 * (shares["h2_supply"] + shares["donor_supply"]) + (1.0 - RADIOLYSIS_H2_ESCAPE) * radiolysis(values)


def held_um(source, reductants, values):
    """
    The dissolved oxygen, micromolar, that a source and a sink (mol of O2 per m2 per year, `Balance.source` and `Balance.sink`) hold in the ocean of `values` at the end of a pulse interval: the surplus over the interval in the water the ocean holds
    (`ocean_depth_km`, 1e6 litres under a m2 for each km), plus the stock the ice accreted with that the reductants the source does not cover have not drained. 0 for an ocean with no depth.
    """
    litres_per_m2 = max(0.0, values.get("ocean_depth_km") or 0.0) * 1.0e6
    if litres_per_m2 <= 0.0:
        return 0.0
    banked = max(0.0, source - reductants) * PULSE_INTERVAL_YR / litres_per_m2 * 1.0e6
    inherited = INHERITED_O2_RETAINED * COMET_O2_PER_WATER * WATER_MOL_PER_L * 1.0e6
    drained = max(0.0, reductants - source) * BODY_AGE_YR / litres_per_m2 * 1.0e6
    return banked + max(0.0, inherited - drained)


def held_kpa(source, reductants, values):
    """The partial pressure, kPa, of the oxygen `held_um` stands for (Henry's law at 0 C), at most the pressure the air law lets the air hold."""
    return min(mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA / 1000.0, held_um(source, reductants, values) / HENRY_UM_PER_KPA)


def balance(values):
    """The `Balance` of the sealed ocean of the body whose column-named `values` are given; `None` for a body that is not a sealed ocean."""
    if not is_sealed(values):
        return None
    lid, water, faults, reductants = lid_delivery(values), radiolysis(values), fracture(values), sink(values)
    source = lid + water + faults
    um = held_um(source, reductants, values)
    return Balance(lid, water, faults, reductants, um, held_kpa(source, reductants, values))


def ocean_o2_kpa(values):
    """The oxygen partial pressure, kPa, of a sealed ocean (see `balance`), or `None` for a body whose ocean is not under a lid."""
    held = balance(values)
    return None if held is None else held.kpa
