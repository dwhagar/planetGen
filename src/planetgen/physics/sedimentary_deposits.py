# planetgen/physics/sedimentary_deposits.py

"""
The sedimentary and biogenic deposits of a rocky body (GEN.1378, the slice of GEN.1379's design that reads the life timeline)
===========================================================================================================================

`docs/design/composition/ore-processes.md` rows 9, 16 and 17 (banded iron formation, phosphorite, evaporite) and the four biogenic rocks the item names (limestone, chert, coal, kerogen) are
rows of the same table the magmatic types fill (`ore_deposits`): `planet_deposit_types`, `moon_deposit_types` and `rogue_planet_deposit_types`. They are assessed **after the life pass**
(`generation.life.apply_life_data`), because four of them read the timeline of life (`life_age_*_gyr`) and the rest are cheap to keep in one place. A body whose life never reached a stage the
type needs has the type `absent_by_state`: a deposit that needs the life timeline is not made without it.

**The settings.** Each type has a setting, a function of the state that returns `None` where the gate fails, else the opportunity (a count scale), the window in which the setting stood (Gyr
of the body's age) and the cap the host puts on the grade:

* **Banded iron formation** (`banded_iron_formation`): a neutral ocean under starlight with iron to give (the rock's FeO), from the first oxygenic photosynthesis (`life_age_photosynthesis_gyr`,
  the oxygen source that oxidises dissolved Fe(II)) for `BIF_SPAN_GY` Gyr, the time the deep ocean takes to turn oxic and shut the iron off. An acid sulfate, soda, brine or sealed ocean has none.
  The abiotic route (ultraviolet oxidation of Fe(II), hydrogen escape) is not modelled: it oxidises too little iron to lay a formation (open point 1 of the design, decided here).
* **Phosphorite** (`phosphorite`): a surface ocean with life (from `BIF`'s same first photosynthesis plus `PHOSPHORITE_LAG_GY`) and phosphate to upwell: the ocean's phosphate of
  `phosphorus_budget` over `PHOSPHORITE_MIN_UMOL_L`, scaling the count; the grade is capped by the crust's P2O5 times `PHOSPHORITE_CONCENTRATION`.
* **Evaporite** (`evaporite`): a liquid surface ocean (open, not under ice) with land basins at a temperature that dries them; the salts are the sequence of the body's own brine
  (`surface_composition.dried_salts`, GEN.1366), the grade is the share of the dominant salt in the bed. It needs no life.
* **Limestone** (`limestone`): a neutral ocean over a shelf with carbon dioxide in the air; abiotic from the start of the ocean's age window, richer once complex cells (calcifiers) exist.
* **Chert** (`chert`): biogenic silica, from the first complex cells (`life_age_complex_cells_gyr`); the abiotic Precambrian chert rides with the iron formation and is not counted separately.
* **Coal** (`coal`): a wet land with flora (`life_flora` at least `COAL_MIN_FLORA`) from the time of multicellular life (`life_age_multicellularity_gyr`).
* **Kerogen** (`kerogen`): organic-rich marine mud, from the first photosynthesis plus `KEROGEN_LAG_GY` in a surface or ice-covered ocean.

**Endowment** is as the magmatic types': a Poisson count of deposits (`deposits` times the opportunity times the window's share of the type's reference window) of lognormal size.
Every constant is a stand-in in section 6.223 of the gap-filling plan (SD1 to SD9); the grade ranges and the per-deposit tonnage are the design file's (recalled) where it has them.
"""

import math
from dataclasses import dataclass
from typing import Callable, NamedTuple, Optional

from planetgen.physics import mixture, ore_deposits, phosphorus_budget, surface_composition
from planetgen.util import draw

BIF, PHOSPHORITE, EVAPORITE, LIMESTONE, CHERT, COAL, KEROGEN = ("banded_iron_formation", "phosphorite", "evaporite", "limestone", "chert", "coal", "kerogen")

OPEN_OCEAN = ("surface ocean", "hycean")
"""tuple: The hydrospheres whose ocean has a free surface in the starlight and an air above it."""
WET_OCEAN = OPEN_OCEAN + ("ice-covered ocean",)
"""tuple: The hydrospheres with a liquid ocean the life reaches (an ice cover lets the light through, thinly)."""

EARTH_FEO = 0.081
EARTH_OCEAN_FRACTION = 0.71
EARTH_LAND_FRACTION = 0.29
"""float: The Earth's mantle FeO (mass fraction of the rock), ocean share and land share: the references the opportunities scale from."""

BIF_SPAN_GY = 1.2
"""float: How long the iron formations are laid after the first oxygenic photosynthesis, Gyr, before the deep ocean is oxic (the Earth: 3.0 to 1.85 Ga). Stand-in SD1, tier T6, L3."""
BIF_DEPOSITS = 12.0
BIF_FE_PER_FEO = 0.7773
BIF_CRUST_OVER_BULK = 0.9
BIF_CONCENTRATION = 11.0
"""Stand-in SD2, tier T6, L3: the count of formations on an Earth-like ocean, and the design's Fe concentration factor 11 on the crust's Fe (0.9 of the bulk silicate's) to cap the grade."""
PHOSPHORITE_LAG_GY = 0.5
PHOSPHORITE_MIN_UMOL_L = 1.0
PHOSPHORITE_DEPOSITS = 25.0
PHOSPHORITE_CONCENTRATION = 200.0
"""Stand-in SD3, tier T6, L3: the delay of productive upwelling after the first photosynthesis, the least ocean phosphate, the deposits of an Earth-like ocean, and the factor from the crust's P2O5 to the ore's."""
EVAPORITE_DEPOSITS = 40.0
EVAPORITE_MIN_TEMPERATURE_K = 270.0
EVAPORITE_WINDOW = (0.2, 1.0)
"""Stand-in SD4, tier T6, L3: the basins of an Earth-like land, the coldest mean temperature at which a basin dries, and the share of the body's age in which evaporites formed."""
LIMESTONE_DEPOSITS = 120.0
LIMESTONE_MIN_CO2_PA = 1.0
LIMESTONE_BIOTIC_BOOST = 2.0
"""Stand-in SD5, tier T7, L3: the shelf deposits of an Earth-like ocean, the least carbon dioxide in the air, and the boost once complex cells (calcifiers) exist."""
CHERT_DEPOSITS = 60.0
"""Stand-in SD6, tier T7, L3: the chert beds of an Earth-like ocean once siliceous life exists."""
COAL_MIN_FLORA = 3
COAL_DEPOSITS = 80.0
"""Stand-in SD7, tier T7, L3: the least `life_flora` score (0 to 9) that makes peat swamps, and the coal fields of an Earth-like land."""
KEROGEN_LAG_GY = 0.3
KEROGEN_DEPOSITS = 60.0
"""Stand-in SD8, tier T7, L3: the delay of organic-rich mud after the first photosynthesis and the source-rock basins of an Earth-like ocean."""
REFERENCE_WINDOW_GY = 2.0
"""float: The window length (Gyr) at which a type's `deposits` count applies in full; a shorter window scales it down, a longer one is capped at it. Stand-in SD9, tier T7, L3."""

SEDIMENTARY_IDS = (BIF, PHOSPHORITE, EVAPORITE, LIMESTONE, CHERT, COAL, KEROGEN)
"""tuple: The ids this module stores in `ore_deposits.ATTRIBUTE` (the magmatic ids stay as the ore pass made them)."""


class State(NamedTuple):
    """What the sedimentary types read of a finished rocky body (all stored values; no draw)."""

    age_gy: float
    hydrosphere: Optional[str]
    ocean_class: Optional[str]
    ocean_fraction: float
    land_fraction: float
    temperature_k: float
    co2_pa: float
    feo: float
    crust_p2o5_wt: float
    phosphate_umol_l: float
    photosynthesis_gy: Optional[float]
    complex_cells_gy: Optional[float]
    multicellular_gy: Optional[float]
    flora: int
    salts: dict


class Setting(NamedTuple):
    """The setting of a type on a body: a count scale, the window it stood in (Gyr of the body's age), the host's cap on the grade and, where the type fixes them, its key and its grade."""

    opportunity: float
    start_gy: float
    end_gy: float
    cap: float
    key: Optional[str] = None
    grade: Optional[float] = None


@dataclass(frozen=True)
class SedimentaryType:
    """One type: its id, the key it grades (element or oxide), the grade range, the deposits of an Earth-like body, the median ore of one and its lognormal spread, and its setting."""

    id: str
    key_element: str
    grade_range: tuple
    deposits: float
    median_ore_kg: float
    sigma: float
    setting: Callable


def _after(start_gy, state, span_gy=None):
    """The window from `start_gy` to the body's age (or `span_gy` later, if sooner); `None` where it has not begun."""
    if start_gy is None or start_gy >= state.age_gy:
        return None
    end = state.age_gy if span_gy is None else min(state.age_gy, start_gy + span_gy)
    return start_gy, end


def _bif(state):
    window = _after(state.photosynthesis_gy, state, BIF_SPAN_GY)
    if window is None or state.hydrosphere not in WET_OCEAN or state.ocean_class != "neutral":
        return None
    iron = state.feo / EARTH_FEO
    cap = min(1.0, BIF_CONCENTRATION * BIF_CRUST_OVER_BULK * BIF_FE_PER_FEO * state.feo)
    return Setting(state.ocean_fraction / EARTH_OCEAN_FRACTION * iron, *window, cap)


def _phosphorite(state):
    window = _after(None if state.photosynthesis_gy is None else state.photosynthesis_gy + PHOSPHORITE_LAG_GY, state)
    if window is None or state.hydrosphere not in OPEN_OCEAN or state.phosphate_umol_l < PHOSPHORITE_MIN_UMOL_L:
        return None
    supply = min(5.0, state.phosphate_umol_l / phosphorus_budget.EARTH_OCEAN_UMOL_L)
    return Setting(state.ocean_fraction / EARTH_OCEAN_FRACTION * supply, *window, min(1.0, PHOSPHORITE_CONCENTRATION * state.crust_p2o5_wt / 100.0))


def _evaporite(state):
    if state.hydrosphere not in OPEN_OCEAN or state.land_fraction <= 0.0 or state.temperature_k < EVAPORITE_MIN_TEMPERATURE_K or not state.salts:
        return None
    key = max(state.salts, key=lambda name: (state.salts[name], name))
    return Setting(state.land_fraction / EARTH_LAND_FRACTION, EVAPORITE_WINDOW[0] * state.age_gy, EVAPORITE_WINDOW[1] * state.age_gy, 1.0, "salt", state.salts[key])


def _limestone(state):
    if state.hydrosphere not in WET_OCEAN or state.ocean_class != "neutral" or state.co2_pa < LIMESTONE_MIN_CO2_PA:
        return None
    biotic = LIMESTONE_BIOTIC_BOOST if state.complex_cells_gy is not None and state.complex_cells_gy < state.age_gy else 1.0
    return Setting(state.ocean_fraction / EARTH_OCEAN_FRACTION * biotic, 0.1 * state.age_gy, state.age_gy, 1.0)


def _chert(state):
    window = _after(state.complex_cells_gy, state)
    if window is None or state.hydrosphere not in WET_OCEAN:
        return None
    return Setting(state.ocean_fraction / EARTH_OCEAN_FRACTION, *window, 1.0)


def _coal(state):
    window = _after(state.multicellular_gy, state)
    if window is None or state.flora < COAL_MIN_FLORA or state.hydrosphere not in OPEN_OCEAN or state.land_fraction <= 0.0:
        return None
    return Setting(state.land_fraction / EARTH_LAND_FRACTION, *window, 1.0)


def _kerogen(state):
    window = _after(None if state.photosynthesis_gy is None else state.photosynthesis_gy + KEROGEN_LAG_GY, state)
    if window is None or state.hydrosphere not in WET_OCEAN:
        return None
    return Setting(state.ocean_fraction / EARTH_OCEAN_FRACTION, *window, 1.0)


TYPES = (
    SedimentaryType(BIF, "Fe", (0.25, 0.35), BIF_DEPOSITS, 1.0e13, 1.0, _bif),
    SedimentaryType(PHOSPHORITE, "P2O5", (0.18, 0.35), PHOSPHORITE_DEPOSITS, 1.0e12, 1.5, _phosphorite),
    SedimentaryType(EVAPORITE, "salt", (0.5, 1.0), EVAPORITE_DEPOSITS, 1.0e12, 2.0, _evaporite),
    SedimentaryType(LIMESTONE, "CaCO3", (0.85, 0.98), LIMESTONE_DEPOSITS, 1.0e13, 1.5, _limestone),
    SedimentaryType(CHERT, "SiO2", (0.90, 0.99), CHERT_DEPOSITS, 1.0e12, 1.5, _chert),
    SedimentaryType(COAL, "C", (0.60, 0.90), COAL_DEPOSITS, 3.0e12, 1.5, _coal),
    SedimentaryType(KEROGEN, "C", (0.02, 0.20), KEROGEN_DEPOSITS, 1.0e13, 1.5, _kerogen),
)
"""tuple: The seven types. Grade ranges and the iron formation's, phosphorite's and evaporite's tonnage (1e9 to 1e11 t, 1e8 to 1e10 t, 1e8 to 1e12 t per deposit) are the design file's, recalled;
the medians are inside those ranges (the others are stand-ins SD5 to SD8). Stand-in tier T6 to T7, L3."""
TYPE_BY_ID = {kind.id: kind for kind in TYPES}


def state_of(body):
    """The `State` of a finished rocky planet or moon, or `None` for a body the ore pass did not assess (a giant, a belt, a row from before GEN.699). No draws."""
    if getattr(body, ore_deposits.ATTRIBUTE, None) is None or ore_deposits.state_of(body) is None:
        return None
    partials = mixture.partial_pressures_pa(body)
    inputs = surface_composition.inputs_of(body)
    shares = surface_composition.dried_salts(inputs) if inputs.rock_oxides is not None else {"halite": 1.0}
    return State(age_gy=getattr(body.star, "age", None) or 0.0, hydrosphere=getattr(body, "hydrosphere", None), ocean_class=getattr(body, "ocean_class", None),
                 ocean_fraction=getattr(body, "ocean_fraction", None) or 0.0, land_fraction=getattr(body, "land_fraction", None) or 0.0,
                 temperature_k=getattr(body, "surface_temperature", None) or 0.0, co2_pa=partials.get("CO2", 0.0),
                 feo=getattr(body, "rock_feo_fraction", None) or 0.0, crust_p2o5_wt=getattr(body, "crust_p2o5_wt", None) or 0.0,
                 phosphate_umol_l=getattr(body, "ocean_phosphate_umol_l", None) or 0.0,
                 photosynthesis_gy=getattr(body, "life_age_photosynthesis_gyr", None), complex_cells_gy=getattr(body, "life_age_complex_cells_gyr", None),
                 multicellular_gy=getattr(body, "life_age_multicellularity_gyr", None), flora=getattr(body, "life_flora", None) or 0, salts=shares)


def stream_of(body, kind):
    """The body's own stream for one type: seeded from what is final about the body and the type, so no draw of the run's stream moves and adding a type moves no other."""
    star = getattr(body, "star", None)
    return draw.body_stream_of(f"sediment:{kind}", body, slot=str(bool(getattr(body, 'is_moon', False))))      # GEN.1547: keyed to the body's birth, not its radius, orbit or the star's age


def assess(kind, state, stream):
    """
    The stored row of one type for a body (`ore_deposits.ROW_FIELDS`): `absent_by_state` where the setting fails or the host's cap is under the bottom of the grade range, `absent_by_chance`
    where the count is zero, else `present` with the grade (drawn in range, held to the cap), the endowment and the window.
    """
    row = {"deposit_type": kind.id, "status": ore_deposits.ABSENT_BY_STATE, "endowment_kg": None, "grade": None, "key_element": kind.key_element,
           "age_window_gyr_start": None, "age_window_gyr_end": None, "capped_by_host": 0}
    setting = kind.setting(state)
    if setting is None:
        return row
    low, high = kind.grade_range
    if setting.cap < low or setting.end_gy <= setting.start_gy:
        return row
    window_share = min(1.0, (setting.end_gy - setting.start_gy) / REFERENCE_WINDOW_GY)
    count = ore_deposits.poisson(stream, kind.deposits * setting.opportunity * window_share)
    drawn = stream.uniform(low, high)
    if count == 0:
        row["status"] = ore_deposits.ABSENT_BY_CHANCE
        return row
    grade = setting.grade if setting.grade is not None else min(drawn, setting.cap)
    row.update(status=ore_deposits.PRESENT, grade=grade, endowment_kg=ore_deposits.lognormal_sum(kind.median_ore_kg, kind.sigma, count, stream),
               capped_by_host=int(setting.grade is None and drawn > setting.cap), age_window_gyr_start=setting.start_gy, age_window_gyr_end=setting.end_gy)
    if setting.key is not None:
        row["key_element"] = setting.key
    return row


def generate(body):
    """
    Adds the seven rows to `ore_deposits.ATTRIBUTE` of a finished rocky planet or moon (call after the life pass; no draw of the run's stream), replacing any it already holds so a body made
    again or edited has each once. A body the ore pass did not assess keeps `None`.
    """
    state = state_of(body)
    if state is None:
        return
    kept = [row for row in getattr(body, ore_deposits.ATTRIBUTE) if row["deposit_type"] not in SEDIMENTARY_IDS]
    setattr(body, ore_deposits.ATTRIBUTE, kept + [assess(kind, state, stream_of(body, kind.id)) for kind in TYPES])
