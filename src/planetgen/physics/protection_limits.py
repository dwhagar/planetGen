# planetgen/physics/protection_limits.py

"""
Every limit of the human visitor's habitability bands and protective gear in one
table (GEN.208), so that a later per-species limit (the canon has races that live in
conditions humans cannot) replaces a row here and nothing else.

`habitability` builds its band edges from these rows, and `protection` (the gear a
planet needs, and the sentence that says why) reads the same rows, so the score, the
kit and the description cannot disagree.

Each row has the value in use and, where there are two, the value the built table
had before the research notes arrived (`built`) and the one in the notes (`notes`:
Environmental.Extremes.and.Protection.Research.md and
Thermal.Respiratory.Protection.Research.md). Where they differ, Boss's rule
(2026-10-14 13:58Z) is that the more permissive one wins, the one that lets a person
go with lighter gear or over a wider range, because the limits will be widened later
anyway. `CHOSEN_PAIRS` lists those choices. Units: pressures kPa, temperatures deg C
(dry bulb unless named), doses Sv/yr.

This module imports nothing from the package, so anything may read it.
"""

from collections import namedtuple

Limit = namedtuple("Limit", "name value unit meaning built notes looser", defaults=(None,))
"""A limit: `name` (the key), `value` in use, `unit`, `meaning`, the `built` and `notes` values (`None`
where that source has none) and `looser`, which direction is the more permissive when the two differ
(`"lower"` or `"higher"`)."""

LIMIT_ROWS = (
    # --- Pressure ---
    Limit("armstrong_kpa", 6.3, "kPa", "water boils at body temperature below this: a full gas-pressurised suit", 6.3, 6.3),
    Limit("mask_min_kpa", 14.3, "kPa", "least ambient pressure at which a pure-oxygen mask still gives 8 kPa of inspired O2; "
          "below it a counterpressure or pressurised suit is needed", 14.3, 20.0, "lower"),
    Limit("pressure_no_gear_kpa", 50.0, "kPa", "least pressure with no gear at all (Blue)", 50.0, 70.0, "lower"),
    Limit("pressure_blue_max_kpa", 250.0, "kPa", "greatest pressure with no gear at all (Blue)", 250.0, None),
    Limit("pressure_narcosis_kpa", 400.0, "kPa", "nitrogen narcosis starts: the greatest pressure of modular gear (Green)", 400.0, 400.0),
    Limit("pressure_hyperbaric_max_kpa", 10_000.0, "kPa", "greatest pressure of a hyperbaric pressure suit (Yellow); "
          "past it a rigid armoured suit", 10_000.0, 5_000.0, "higher"),
    Limit("pressure_outer_kpa", 100_000.0, "kPa", "the score reaches 0 here", 100_000.0, None),
    Limit("pressure_low_outer_kpa", 0.63, "kPa", "the score reaches 0 here at the low end", 0.63, None),
    # --- Oxygen ---
    Limit("inspired_o2_min_kpa", 8.0, "kPa", "least inspired O2 for consciousness; below it a breathing mask supplies it", 8.0, 19.5, "lower"),
    Limit("inspired_o2_chronic_kpa", 50.0, "kPa", "most inspired O2 for months", 50.0, None),
    Limit("inspired_o2_acute_kpa", 160.0, "kPa", "acute CNS oxygen toxicity (1.6 bar)", 160.0, None),
    # --- Temperature ---
    Limit("cold_outer_c", -120.0, "C", "the score reaches 0 here", -120.0, None),
    Limit("cold_suit_c", -50.0, "C", "colder than this a sealed suit with active heating (Red); down to it an insulated suit (Yellow)", -50.0, -40.0, "lower"),
    Limit("cold_clothing_c", -20.0, "C", "colder than this cold-weather clothing is not enough (Green below Blue)", -20.0, -20.0),
    Limit("comfort_low_c", 0.0, "C", "least temperature with no gear (Blue)", 0.0, None),
    Limit("comfort_high_c", 31.0, "C", "greatest temperature with no gear (Blue)", 31.0, None),
    Limit("heat_garment_c", 45.0, "C", "hotter than this an actively cooled sealed suit (Yellow); below it cooling garments (Green)", 45.0, 45.0),
    Limit("heat_max_c", 122.0, "C", "the limit of carbon-based life (Takai et al., 2008): the score's Yellow edge", 122.0, None),
    Limit("heat_shielded_c", 140.0, "C", "hotter than this a mask-only visitor's skin burns at the seal: a heat-shielded suit with a liquid cooling garment", None, 140.0),
    Limit("hot_outer_c", 200.0, "C", "the score reaches 0 here", 200.0, None),
    Limit("wet_bulb_c", 31.0, "C", "uncompensable wet-bulb temperature: cooling garments needed at or above it", 31.0, 30.0, "higher"),
    # --- Radiation ---
    Limit("dose_blue_sv_yr", 0.02, "Sv/yr", "the ICRP occupational limit averaged over five years (Blue)", 0.02, None),
    Limit("dose_green_sv_yr", 0.1, "Sv/yr", "a visit is still allowed (Green)", 0.1, None),
    Limit("dose_yellow_sv_yr", 1.0, "Sv/yr", "the threshold of radiation sickness: shelter or shielding needed (Yellow); "
          "past it a hardened habitat (Red)", 1.0, 1.0),
    Limit("dose_outer_sv_yr", 1000.0, "Sv/yr", "the score reaches 0 here", 1000.0, None),
    Limit("dose_single_year_sv", 0.05, "Sv", "ICRP's limit for any one year", 0.05, None),
)

LIMITS = {row.name: row for row in LIMIT_ROWS}
"""dict: The limits by name."""


def limit(name):
    """The value in use of the limit `name`."""
    return LIMITS[name].value


def chosen_pairs():
    """
    The limits where the built table and the research notes gave two different values, as
    `(name, built, notes, value in use)`. The value in use is always the more permissive of the two
    (`Limit.looser`); the merge note lists these.
    """
    out = []
    for row in LIMIT_ROWS:
        if row.built is None or row.notes is None or row.built == row.notes:
            continue
        pick = min if row.looser == "lower" else max
        assert row.value == pick(row.built, row.notes), row.name
        out.append((row.name, row.built, row.notes, row.value))
    return out
