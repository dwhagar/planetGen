# planetgen/physics/hydrothermal_deposits.py

"""
The hydrothermal and sedimentary-hosted deposits of a rocky body (GEN.1381, the second slice of GEN.1379's design)
==================================================================================================================

`docs/design/composition/ore-processes.md` rows 1 to 5 and 15: porphyry Cu-Mo-Au, volcanogenic massive sulfide (VMS), SEDEX and MVT Zn-Pb, orogenic gold, epithermal Au-Ag and roll-front uranium.
They are rows of the table the magmatic types fill (`ore_deposits`): `planet_deposit_types`, `moon_deposit_types` and `rogue_planet_deposit_types`, so this slice stores no new column and
needs no migration. They are assessed **after the life pass**, next to the sedimentary slice, because the uranium type reads the oxygen the air ended with and the reductant that life leaves.

The record, the count law (a Poisson count of deposits scaled by the opportunity and the window, each deposit a lognormal ore mass) and the status flag are the sedimentary slice's: this module
supplies its own types and settings and calls `sedimentary_deposits.assess`.

**The settings** (rule 2 of the design: liquid water in the crust and a heat engine; a mobile lid for the arc and collision types). Each returns `None` where the gate fails:

* **Porphyry** (`porphyry_cu_mo_au`): a mobile lid, liquid water, an oxidized mantle (a sulfate-rich arc magma, as the magmatic sulfide gate has it the other way) and a crust thick enough to hold a
  1 to 4 km intrusion (`PORPHYRY_MIN_CRUST_KM`). The grade is Cu, capped by the star's copper.
* **VMS** (`vms_massive_sulfide`): a mobile lid and a liquid ocean (spreading under water) with heat flow of at least `VMS_MIN_HEAT_RATIO` of the Earth's. The grade is Zn.
* **SEDEX and MVT** (`sedex_mvt_zn_pb`): liquid water, a basin (land and an ocean to fill it), a crust sulfur index of `SEDEX_MIN_SULFUR` or more (sulfate reduction) and sediment old enough to
  hold brines (`SEDEX_MIN_AGE_GY`). The grade is Zn.
* **Orogenic gold** (`orogenic_gold`): a mobile lid, liquid water, a continental crust (a granitoid share and a crust of `OROGENIC_MIN_CRUST_KM`) to metamorphose. The grade is Au.
* **Epithermal** (`epithermal_au_ag`): a mobile lid, liquid water and volcanism (a processed mantle of `EPITHERMAL_MIN_MAGMA_RATIO`), kept only in the youngest part of the body's age (these
  deposits erode away). The grade is Au.
* **Roll-front uranium** (`roll_front_uranium`): oxidizing groundwater (oxygen in the air over `ROLL_FRONT_MIN_O2_PA`), land, a granitoid source rock and a reductant (the crust's sulfur or
  life). The grade is U3O8, capped by the uranium the body's rock stores (`c_u0_ppm`) times `U_CRUST_OVER_BULK`.

**Endowment and grade.** The grade is drawn inside the type's range and held to `concentration * crustal mean`: the crustal mean of copper and zinc is the Earth's, scaled by the star's Cu/Fe and
Zn/Fe over the reference mix's (`ore_deposits.over_iron`); gold has no star ratio and follows the iron. The opportunity is the surface area over the Earth's times the heat flow, the ocean or
the land the type needs, over the Earth's. Every constant is a stand-in in section 6.224 of the gap-filling plan (HT1 to HT7). The grade ranges, the concentration factors and the tonnage
range of each type are the design file's, recalled.
"""

from typing import NamedTuple, Optional

from planetgen.physics import constants, interior_heat, mixture, ore_deposits, sedimentary_deposits as sd

PORPHYRY, VMS, SEDEX, OROGENIC, EPITHERMAL, ROLL_FRONT = ("porphyry_cu_mo_au", "vms_massive_sulfide", "sedex_mvt_zn_pb", "orogenic_gold", "epithermal_au_ag", "roll_front_uranium")
HYDROTHERMAL_IDS = (PORPHYRY, VMS, SEDEX, OROGENIC, EPITHERMAL, ROLL_FRONT)
"""tuple: The ids this module stores in `ore_deposits.ATTRIBUTE` (the magmatic and sedimentary ids stay as their passes made them)."""

EARTH_HEAT_FLUX_W_M2 = 0.087
EARTH_OCEAN_FRACTION = 0.71
EARTH_LAND_FRACTION = 0.29
EARTH_O2_PA = 21278.0
"""float: The Earth's mean heat flow, ocean share, land share and oxygen partial pressure: the references the opportunities scale from. Row HT1, tier T3, L0."""

PORPHYRY_MIN_CRUST_KM = 20.0
OROGENIC_MIN_CRUST_KM = 25.0
VMS_MIN_HEAT_RATIO = 0.3
SEDEX_MIN_SULFUR = 0.3
SEDEX_MIN_AGE_GY = 0.5
EPITHERMAL_MIN_MAGMA_RATIO = 0.2
ROLL_FRONT_MIN_O2_PA = 1000.0
"""Stand-in HT2, tier T7, L3: the thinnest crust that holds a porphyry stock and a collision belt, the coldest ridge that drives a hydrothermal cell, the least sulfur for a brine to reduce sulfate,
the youngest age with basin sediment, the least volcanism for a shallow system and the least oxygen for an oxidizing groundwater."""
HEAT_RATIO_CAP = 3.0
MAGMA_RATIO_CAP = 3.0
"""float: The ceilings on the heat-flow and volcanism factors of the opportunity, so a very hot body does not multiply the count without bound. Stand-in HT3, tier T7, L3."""

CRUST_PPM = {"Cu": 60.0e-6, "Zn": 70.0e-6, "Au": 4.0e-9}
"""dict: The Earth's continental crust mean (mass fraction) of copper, zinc and gold, the references the star's ratios scale. Stand-in HT4, tier T4 (recalled), L0."""
CONCENTRATION = {"Cu": 140.0, "Zn": 750.0, "Au": 2300.0, "U": 740.0}
"""dict: The factor by which the process concentrates the crustal mean into ore (the design file's, recalled). Stand-in HT4."""
U_CRUST_OVER_BULK = 80.0
U3O8_PER_U = 1.1792
"""float: The crust's uranium over the bulk silicate's (an incompatible element, the crust holds about 80 times the mantle's mean) and the U3O8 mass per unit mass of uranium. Stand-in HT5, tier T6, L2 (the Earth's 1.7 ppm over 0.02 ppm)."""

WINDOWS = {PORPHYRY: (0.5, 1.0), VMS: (0.3, 1.0), SEDEX: (0.3, 1.0), OROGENIC: (0.3, 1.0), EPITHERMAL: (0.9, 1.0), ROLL_FRONT: (0.6, 1.0)}
"""dict: The share of the body's age in which each type formed and survived. Stand-in HT6, tier T7, L3."""

PORPHYRY_DEPOSITS, VMS_DEPOSITS, SEDEX_DEPOSITS, OROGENIC_DEPOSITS, EPITHERMAL_DEPOSITS, ROLL_FRONT_DEPOSITS = 30.0, 100.0, 80.0, 100.0, 150.0, 120.0
"""Stand-in HT7, tier T7, L3: the deposits of an Earth-like body over a window of the reference length (the count the Earth's known deposits suggest, recalled)."""


class State(NamedTuple):
    """What the hydrothermal types read of a finished rocky body (all stored values; no draw)."""

    age_gy: float
    lid: str
    redox: Optional[str]
    wet: bool
    ocean: bool
    ocean_fraction: float
    land_fraction: float
    radius_km: float
    heat_flux_w_m2: float
    crust_km: float
    sulfur_index: float
    felsic_share: float
    magma_ratio: float
    o2_pa: float
    reductant: bool
    crust_ppm: dict
    u_ppm: float


def _heat_ratio(state):
    return min(HEAT_RATIO_CAP, state.heat_flux_w_m2 / EARTH_HEAT_FLUX_W_M2)


def _area_ratio(state):
    return (state.radius_km / constants.EARTH_RADIUS_KM) ** 2


def _cap(state, element):
    """The highest grade the host supplies: the concentration factor times the crustal mean of the element."""
    return min(1.0, CONCENTRATION[element] * state.crust_ppm[element])


def _setting(state, kind, opportunity, cap):
    start, end = WINDOWS[kind]
    return sd.Setting(opportunity, start * state.age_gy, end * state.age_gy, cap)


def _mobile_wet(state):
    return state.lid == interior_heat.MOBILE and state.wet


def _porphyry(state):
    if not _mobile_wet(state) or state.redox != "oxidized" or state.crust_km < PORPHYRY_MIN_CRUST_KM:
        return None
    return _setting(state, PORPHYRY, _area_ratio(state) * _heat_ratio(state), _cap(state, "Cu"))


def _vms(state):
    if not _mobile_wet(state) or not state.ocean or state.heat_flux_w_m2 < VMS_MIN_HEAT_RATIO * EARTH_HEAT_FLUX_W_M2:
        return None
    return _setting(state, VMS, _area_ratio(state) * _heat_ratio(state) * state.ocean_fraction / EARTH_OCEAN_FRACTION, _cap(state, "Zn"))


def _sedex(state):
    if not state.wet or not state.ocean or state.land_fraction <= 0.0 or state.sulfur_index < SEDEX_MIN_SULFUR or state.age_gy < SEDEX_MIN_AGE_GY:
        return None
    return _setting(state, SEDEX, _area_ratio(state) * state.ocean_fraction / EARTH_OCEAN_FRACTION * min(1.5, state.sulfur_index / 0.8), _cap(state, "Zn"))


def _orogenic(state):
    if not _mobile_wet(state) or state.felsic_share <= 0.0 or state.crust_km < OROGENIC_MIN_CRUST_KM:
        return None
    return _setting(state, OROGENIC, _area_ratio(state) * _heat_ratio(state), _cap(state, "Au"))


def _epithermal(state):
    if not _mobile_wet(state) or state.magma_ratio < EPITHERMAL_MIN_MAGMA_RATIO:
        return None
    return _setting(state, EPITHERMAL, _area_ratio(state) * min(MAGMA_RATIO_CAP, state.magma_ratio), _cap(state, "Au"))


def _roll_front(state):
    if state.o2_pa < ROLL_FRONT_MIN_O2_PA or state.land_fraction <= 0.0 or not state.wet or state.felsic_share <= 0.0 or not state.reductant:
        return None
    cap = min(1.0, CONCENTRATION["U"] * state.u_ppm * 1.0e-6 * U_CRUST_OVER_BULK * U3O8_PER_U)
    return _setting(state, ROLL_FRONT, _area_ratio(state) * state.land_fraction / EARTH_LAND_FRACTION * min(1.0, state.o2_pa / EARTH_O2_PA), cap)


TYPES = (
    sd.SedimentaryType(PORPHYRY, "Cu", (0.003, 0.010), PORPHYRY_DEPOSITS, 1.0e12, 1.5, _porphyry),
    sd.SedimentaryType(VMS, "Zn", (0.03, 0.10), VMS_DEPOSITS, 1.0e10, 1.5, _vms),
    sd.SedimentaryType(SEDEX, "Zn", (0.05, 0.12), SEDEX_DEPOSITS, 1.0e10, 1.5, _sedex),
    sd.SedimentaryType(OROGENIC, "Au", (3.0e-6, 1.0e-5), OROGENIC_DEPOSITS, 1.0e10, 1.5, _orogenic),
    sd.SedimentaryType(EPITHERMAL, "Au", (1.0e-6, 1.0e-5), EPITHERMAL_DEPOSITS, 1.0e9, 1.5, _epithermal),
    sd.SedimentaryType(ROLL_FRONT, "U3O8", (0.0005, 0.003), ROLL_FRONT_DEPOSITS, 1.0e10, 1.5, _roll_front),
)
"""tuple: The six types. Grade ranges, concentration factors and the tonnage range of one deposit (1e8 to 1e10 t, 1e6 to 1e8 t, 1e6 to 1e8 t, 1e6 to 1e8 t, 1e5 to 1e7 t, 1e3 to 1e5 t of U) are the design
file's, recalled; the medians are inside those ranges (stand-ins HT7). Tier T6 to T7, L3."""

def state_of(body):
    """The `State` of a finished rocky planet or moon, or `None` for a body the ore pass did not assess (a giant, a belt, a row from before GEN.699). No draws."""
    ore = ore_deposits.state_of(body) if getattr(body, ore_deposits.ATTRIBUTE, None) is not None else None
    if ore is None or getattr(body, "internal_heat_flux_w_m2", None) is None:
        return None
    star = getattr(body, "star", None)
    crust_ppm = {element: ppm * ore_deposits.over_iron(star, element) for element, ppm in CRUST_PPM.items()}
    hydrosphere = getattr(body, "hydrosphere", None)
    return State(age_gy=ore.age_gy, lid=ore.lid, redox=ore.redox, wet=ore.wet, ocean=hydrosphere in sd.WET_OCEAN, ocean_fraction=getattr(body, "ocean_fraction", None) or 0.0,
                 land_fraction=getattr(body, "land_fraction", None) or 0.0, radius_km=body.radius, heat_flux_w_m2=body.internal_heat_flux_w_m2, crust_km=body.crust_thickness_km or 0.0,
                 sulfur_index=body.crust_sulfur_index or 0.0, felsic_share=ore.felsic_share, magma_ratio=ore.magma_ratio, o2_pa=mixture.partial_pressures_pa(body).get("O2", 0.0),
                 reductant=(body.crust_sulfur_index or 0.0) >= SEDEX_MIN_SULFUR or getattr(body, "life_age_photosynthesis_gyr", None) is not None, crust_ppm=crust_ppm,
                 u_ppm=getattr(body, "c_u0_ppm", None) or 0.0)


def stream_of(body, kind):
    """The body's own stream for one type (`ore_deposits.body_stream`): seeded from what is final about the body and the type, so no draw of the run's stream moves and adding a type moves no other."""
    return ore_deposits.body_stream("hydrothermal", body, kind)


def generate(body):
    """
    Adds the six rows to `ore_deposits.ATTRIBUTE` of a finished rocky planet or moon (call after the life pass; no draw of the run's stream), replacing any it already holds so a body made
    again or edited has each once. A body the ore pass did not assess keeps `None`.
    """
    state = state_of(body)
    if state is None:
        return
    kept = [row for row in getattr(body, ore_deposits.ATTRIBUTE) if row["deposit_type"] not in HYDROTHERMAL_IDS]
    setattr(body, ore_deposits.ATTRIBUTE, kept + [sd.assess(kind, state, stream_of(body, kind.id)) for kind in TYPES])
