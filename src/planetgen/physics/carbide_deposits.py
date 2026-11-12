# planetgen/physics/carbide_deposits.py

"""
The carbide deposits of a carbon-rich rocky body (GEN.1423)
===========================================================

A rocky body whose star has C/O at or above `bulk_draw.CARBON_OVER_OXYGEN` makes carbon rock: silicon carbide and graphite in place of silicate (`carbon_rock`, [R0455]). The ground of such a body is already
the `carbon` type of `surface_composition`; this module adds the two minerals the item names as rows of the same deposit table the other slices fill (`ore_deposits.ATTRIBUTE`: `planet_deposit_types`,
`moon_deposit_types`, `rogue_planet_deposit_types`), so no column is stored and no migration is needed:

* **Moissanite crust** (`moissanite_crust`): the silicon carbide of the carbon rock. The crust is the rock the planet condensed, so the grade is the SiC share of the SiC-and-graphite rock,
  `1 / (1 + g)` with `g` the graphite against the rest of the rock, which the bulk draw takes in `bulk_draw.RANGES["graphite_over_rock"]` (0.05 to 0.40): 71 to 95 percent. The count scales with the
  surface area and the crust thickness the ore pass stored (`crust_thickness_km`).
* **Cohenite metal** (`cohenite_metal`): iron carbide, Fe3C, in the metal of a carbon-rich body. Iron that meets carbon at the iron-carbon eutectic holds up to 4.3 percent carbon, which is
  `CARBON_SATURATED_C_WT / COHENITE_C_WT` of the metal as Fe3C (0.64); the grade is drawn between `COHENITE_GRADE_FLOOR` and that. The count scales with the body's metal mass over the Earth's core share.

Both are `absent_by_state` on a body whose rock is silicate (nothing here is drawn for an oxygen-rich world) and on a body the ore pass did not assess. They use the count law and the record of
`sedimentary_deposits` (a Poisson count of lognormal sizes, on a stream of the body's own), so adding them moves no other draw. Every constant is a stand-in in section 6.235 of the gap-filling plan (CH1 to CH3).
"""

from typing import NamedTuple

from planetgen.physics import bulk_draw, constants, ore_deposits, sedimentary_deposits as sd
from planetgen.util import draw

MOISSANITE, COHENITE = "moissanite_crust", "cohenite_metal"
CARBIDE_IDS = (MOISSANITE, COHENITE)
"""tuple: The ids this module stores in `ore_deposits.ATTRIBUTE`."""

CRUST_REFERENCE_KM = 35.0
EARTH_CORE_MASS_FRACTION = 0.325
"""float: The Earth's continental crust thickness (km) and core mass fraction, the references the counts scale from. Stand-in CH1, tier T4 (recalled), L0."""

CARBON_SATURATED_C_WT = 4.3
COHENITE_C_WT = 6.69
COHENITE_GRADE_FLOOR = 0.30
"""float: Carbon (percent by mass) in iron saturated at the iron-carbon eutectic, carbon in Fe3C (percent, 12.011 / 179.55), and the lowest Fe3C share of the metal drawn. The first is the textbook
eutectic (recalled), the second is stoichiometry, the floor a design edge. Stand-in CH2, tier T4, L2."""

MOISSANITE_DEPOSITS = 40.0
COHENITE_DEPOSITS = 20.0
MOISSANITE_MEDIAN_KG = 1.0e12
COHENITE_MEDIAN_KG = 1.0e11
SPREAD = 1.5
"""Stand-in CH3, tier T7, L3: the deposits over a window of the reference length on a carbon-rich body of the Earth's size, their median ore (kg) and the lognormal spread. No carbon-rich body has a
surveyed endowment; the numbers are set beside the other types' (sedimentary_deposits.TYPES)."""

WINDOWS = {MOISSANITE: (0.0, 1.0), COHENITE: (0.2, 1.0)}
"""dict: The share of the body's age in which each formed: the carbon rock condensed with the body; the metal separated and was carbon-saturated later. Stand-in CH3."""

SIC_GRADE = (1.0 / (1.0 + bulk_draw.RANGES["graphite_over_rock"][1]), 1.0 / (1.0 + bulk_draw.RANGES["graphite_over_rock"][0]))
"""tuple: The SiC share of the carbon rock at the two ends of the graphite range the bulk draw takes (0.714 to 0.952)."""
COHENITE_GRADE = (COHENITE_GRADE_FLOOR, CARBON_SATURATED_C_WT / COHENITE_C_WT)
"""tuple: The Fe3C share of the metal: from the floor to the carbon-saturated share (0.64)."""


class State(NamedTuple):
    """What the carbide types read of a finished rocky body (all stored values; no draw)."""

    age_gy: float
    carbon_rock: bool
    radius_km: float
    crust_km: float
    metal_mass_fraction: float


def _setting(state, kind, opportunity):
    start, end = WINDOWS[kind]
    return sd.Setting(opportunity, start * state.age_gy, end * state.age_gy, 1.0)


def _moissanite(state):
    if not state.carbon_rock or state.crust_km <= 0.0:
        return None
    return _setting(state, MOISSANITE, (state.radius_km / constants.EARTH_RADIUS_KM) ** 2 * state.crust_km / CRUST_REFERENCE_KM)


def _cohenite(state):
    if not state.carbon_rock or state.metal_mass_fraction <= 0.0:
        return None
    return _setting(state, COHENITE, (state.radius_km / constants.EARTH_RADIUS_KM) ** 2 * state.metal_mass_fraction / EARTH_CORE_MASS_FRACTION)


TYPES = (
    sd.SedimentaryType(MOISSANITE, "SiC", SIC_GRADE, MOISSANITE_DEPOSITS, MOISSANITE_MEDIAN_KG, SPREAD, _moissanite),
    sd.SedimentaryType(COHENITE, "Fe3C", COHENITE_GRADE, COHENITE_DEPOSITS, COHENITE_MEDIAN_KG, SPREAD, _cohenite),
)
"""tuple: The two types."""


def state_of(body):
    """The `State` of a finished rocky planet or moon, or `None` for a body the ore pass did not assess (a giant, a belt, a row from before GEN.699). No draws."""
    ore = ore_deposits.state_of(body) if getattr(body, ore_deposits.ATTRIBUTE, None) is not None else None
    if ore is None:
        return None
    return State(age_gy=ore.age_gy, carbon_rock=bool(getattr(body, "carbon_rock", False)), radius_km=body.radius,
                 crust_km=getattr(body, "crust_thickness_km", None) or 0.0, metal_mass_fraction=getattr(body, "metal_mass_fraction", None) or 0.0)


def stream_of(body, kind):
    """The body's own stream for one type: seeded from what is final about the body and the type, so no draw of the run's stream moves and adding a type moves no other."""
    return draw.body_stream_of(f"carbide:{kind}", body, slot=str(bool(getattr(body, 'is_moon', False))))      # GEN.1547: keyed to the body's birth, not its radius, orbit or the star's age


def generate(body):
    """
    Adds the two rows to `ore_deposits.ATTRIBUTE` of a finished rocky planet or moon (call after the ore pass; no draw of the run's stream), replacing any it already holds so a body made again
    or edited has each once. A body the ore pass did not assess keeps `None`.
    """
    state = state_of(body)
    if state is None:
        return
    kept = [row for row in getattr(body, ore_deposits.ATTRIBUTE) if row["deposit_type"] not in CARBIDE_IDS]
    setattr(body, ore_deposits.ATTRIBUTE, kept + [sd.assess(kind, state, stream_of(body, kind.id)) for kind in TYPES])
