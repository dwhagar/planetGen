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
    Limit("cold_suit_c", -50.0, "C", "colder than this an environment suit with active heating (Red); down to it an insulated suit (Yellow)", -50.0, -40.0, "lower"),
    Limit("cold_clothing_c", -20.0, "C", "colder than this cold-weather clothing is not enough (Green below Blue)", -20.0, -20.0),
    Limit("comfort_low_c", 0.0, "C", "least temperature with no gear (Blue)", 0.0, None),
    Limit("comfort_high_c", 31.0, "C", "greatest temperature with no gear (Blue)", 31.0, None),
    Limit("heat_garment_c", 45.0, "C", "hotter than this an actively cooled environment suit (Yellow); below it cooling garments (Green)", 45.0, 45.0),
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
    # --- Survival equipment (GEN.655): where gear is required and where it is only recommended. The gear follows the
    # physics even where the score's bands are looser (Boss, 2026-10-23, decision B) ---
    Limit("suit_required_kpa", 16.0, "kPa", "below this a mask on pure oxygen no longer works: a pressure suit is required", None, None),
    Limit("suit_recommended_kpa", 20.0, "kPa", "from here down a mask on pure oxygen is marginal: a pressure suit is recommended", None, 20.0),
    Limit("inspired_o2_rec_kpa", 16.0, "kPa", "below this inspired O2 (down to the consciousness floor) a gas-mix mask is recommended", None, None),
    Limit("cold_heated_rec_c", -40.0, "C", "from here down a heated sealed suit is recommended (the notes' edge; the table's required edge is colder)", None, -40.0),
    Limit("wet_bulb_rec_c", 30.0, "C", "from here a wet bulb makes cooling clothing recommended (Vecellio's 30.55 less its spread)", None, 30.0),
    Limit("density_rec_g_l", 6.0, "g/L", "best-mix breathing-gas density where work capacity halves: an armoured suit is recommended", None, 6.0),
    Limit("density_max_g_l", 10.0, "g/L", "best-mix breathing-gas density where ventilation fails: an armoured suit is required", None, 10.0),
    Limit("armoured_heat_c", 400.0, "C", "hotter than this an armoured, thermal-capacitor suit (hours of sortie)", None, None),
    Limit("vehicle_heat_c", 500.0, "C", "hotter than this no suit exists: a vehicle", None, None),
    Limit("vehicle_pressure_kpa", 100_000.0, "kPa", "higher than this no suit exists: a vehicle", None, None),
    Limit("frame_gravity_rec_g", 1.25, "g", "from here a powered frame is recommended", None, None),
    Limit("frame_gravity_g", 1.5, "g", "from here a powered frame is required", None, None),
    Limit("vehicle_gravity_g", 3.0, "g", "greater than this no suit or frame helps: a vehicle", None, None),
    Limit("uv_cover_index", 3.0, "x Earth", "a DNA-weighted UV index from here: sun-cover clothing", None, None),
    Limit("uv_suit_index", 30.0, "x Earth", "a DNA-weighted UV index from here: a full-cover suit", None, None),
    Limit("corrosive_ph_low", 2.0, "pH", "a surface liquid below this pH needs a chemical-barrier suit", None, None),
    Limit("corrosive_ph_high", 12.0, "pH", "a surface liquid above this pH needs a chemical-barrier suit", None, None),
    Limit("dust_pressure_kpa", 10.0, "kPa", "interim rule: a rocky body under this pressure has regolith dust that wrecks open seals", None, None),
    Limit("work_hours_outside_yr", 2000.0, "h", "a working year outside, for the suit worn to go out from a shelter", None, None),
)

Toxic = namedtuple("Toxic", "gas chronic_kpa idlh_kpa skin source")
"""A gas a filter can take out: its long-term limit and its IDLH (immediately dangerous to life or health) as partial
pressures, kPa; `skin` is true when it is corrosive to skin and eyes (a chemical-barrier suit, not a filter, past its
IDLH); `source` says where the numbers come from."""

PPM_KPA = 101.325e-6
"""float: One ppm of a gas at 1 atm, as a partial pressure in kPa."""

TOXIC_GASES = (
    Toxic("CO", 0.005, 0.12, False, "built long-term limit; NIOSH IDLH 1,200 ppm"),
    Toxic("H2S", 0.001, 100 * PPM_KPA, False, "built long-term limit; NIOSH IDLH 100 ppm (the built acute limit, 0.05, is looser)"),
    Toxic("SO2", 0.0005, 100 * PPM_KPA, False, "built long-term limit; NIOSH IDLH 100 ppm"),
    Toxic("NH3", 25 * PPM_KPA, 300 * PPM_KPA, True, "recalled: NIOSH REL 25 ppm, IDLH 300 ppm (check at build time)"),
    Toxic("HF", 3 * PPM_KPA, 30 * PPM_KPA, True, "recalled: ACGIH ceiling 3 ppm, IDLH 30 ppm"),
    Toxic("HCl", 5 * PPM_KPA, 45 * PPM_KPA, True, "NIOSH ceiling 5 ppm; IDLH 45 ppm (NIOSH 2025)"),
    Toxic("Cl2", 1 * PPM_KPA, 10 * PPM_KPA, True, "recalled: NIOSH ceiling 0.5 to 1 ppm, IDLH 10 ppm"),
    Toxic("NO2", 1 * PPM_KPA, 20 * PPM_KPA, False, "recalled: NIOSH STEL 1 ppm, IDLH 20 ppm"),
    Toxic("HCN", 4.7 * PPM_KPA, 50 * PPM_KPA, False, "recalled: NIOSH STEL 4.7 ppm, IDLH 50 ppm"),
    Toxic("O3", 0.1 * PPM_KPA, 5 * PPM_KPA, False, "recalled: OSHA 0.1 ppm, IDLH 5 ppm"),
)
"""tuple: The `Toxic` rows. CO, H2S and SO2 are the species the generator stores today; the others are read once the
atmosphere stores every gas (GEN.424, GEN.657). A world's gas that is not a row here is read by its own rule."""

TOXIC = {row.gas: row for row in TOXIC_GASES}
"""dict: The `TOXIC_GASES` rows by gas symbol."""

RadiationStep = namedtuple("RadiationStep", "key name layers")
"""A step of the radiation ladder (decision F): `key`, `name` and `layers`, `((share of the body, g/cm2 of
hydrogen-rich shielding), ...)`: the body is the sum of its layers (a vest covers the torso, 60 percent of the
effective dose, with more)."""

RADIATION_STEPS = (
    RadiationStep("garment", "shielded garment", ((1.0, 1.0),)),
    RadiationStep("vest", "storm vest", ((0.4, 1.0), (0.6, 6.0))),
    RadiationStep("hardened", "hardened suit", ((0.4, 5.0), (0.6, 10.0))),
    RadiationStep("shelter", "shielded shelter", ((1.0, 50.0),)),
    RadiationStep("habitat", "hardened habitat", ((1.0, 300.0),)),
)
"""tuple: The five steps, lightest first (the real gear each stands for is in
docs/design/survival-equipment-fields.md 4.5; the thicknesses are planetGen defaults)."""

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
