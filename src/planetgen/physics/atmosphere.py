"""
Mantle redox and the stored air of planets and moons (GEN.85, GEN.411, GEN.422)
==============================================================================

Each body with air stores the partial pressure, Pa, of every gas in its air (`GAS_FIELD`, a dict of gas-table formula to Pa, kept in
the child tables `planet_gases` and `moon_gases`; DB.71) and its mantle's redox state; the free-text `atmosphere` is written from
them. The air itself is made by `surface_air` (docs/design/composition/surface-pressure-supply-law.md): the mixture draw of the
body's regimes, the pressure the interior has degassed and the sinks have left, the freeze-out of what its temperature condenses and
the greenhouse that gives the temperature, solved together. This module keeps the result (`store_air`), reads it (the accessors
below) and says it in words (`describe`).

1. **Mantle redox.** The oxygen fugacity of the mantle that outgassed the air, in log units above the iron-wustite buffer (`delta_iw`,
   on the one scale of `physics/oxygen_fugacity.py`, the [R1951] buffer), is the mixture draw's (`MixtureDraw.redox_diw`:
   the recipe rows a1 to a4, mass- and age-dependent, [R4154]). Reduced at or below IW, oxidized from
   `REDOX_OXIDIZED_DELTA_IW`. Gas giants have no mantle redox (`None`).
2. **Cold side.** No species can exceed its saturation pressure at the surface temperature (`physics/saturation.py`: the Fray and
   Schmitt solid curves below the triple point, the Wagner liquid curves above it, no limit above the critical point): the excess is
   condensate, and the total surface pressure drops by it (`surface_air.freeze_out`).
"""

from collections.abc import Mapping

from planetgen.physics import gases, saturation

SPECIES = ("o2", "co2", "co", "n2", "ar", "h2", "h2o", "ch4", "h2s", "so2")
"""tuple: The stored species, in column order (`p_<species>_kpa`)."""

GAS_FIELD = "gas_partials_pa"
"""str: The body attribute holding the air: gas-table formula to partial pressure in Pa, only the gases above zero. Stored in
the child tables `planet_gases` and `moon_gases` (DB.71), not in columns."""

ATMOSPHERE_FIELDS = ("mantle_redox", "mantle_delta_iw")
"""tuple: The attributes (and columns) this module sets besides `GAS_FIELD`."""

REDOX_OXIDIZED_DELTA_IW = 2.5
"""float: Above this (about QFM - 1) a mantle counts as oxidized [D]."""

THIN_KPA = 1.0
DENSE_KPA = 500.0
MAJOR_FRACTION = 0.01
"""float: The share of the air from which a gas is named as one of the mix (`describe`; at most four)."""
TRACE_FRACTION = 1e-4
"""float: The share from which a smaller gas is named as a trace (at most two)."""
HUMID_H2O_FRACTION = 0.02

PENDING = "an air not yet drawn"
"""str: What a body whose class has air says until its air is made (`calculate_atmospheric_conditions` replaces it by `describe` of the
solved mixture, GEN.415): a class names no gas, only whether it has air."""


def redox_class(delta_iw):
    """'reduced' (at or below IW), 'oxidized' (from
    `REDOX_OXIDIZED_DELTA_IW`) or 'intermediate'."""
    if delta_iw <= 0.0:
        return "reduced"
    if delta_iw >= REDOX_OXIDIZED_DELTA_IW:
        return "oxidized"
    return "intermediate"


def vapour_pressure_kpa(gas, temperature_k):
    """The most of `gas` (a stored species, `SPECIES`) that can be vapour at `temperature_k`, kPa: its saturation pressure
    (`saturation.saturation`), infinite above its critical temperature."""
    return saturation.pressure_pa(gases.stored(gas).formula, temperature_k) / 1000.0


def gas_name(formula):
    """The name of a gas of the table in words, as the gas table spells it (`gases.gas(formula).name`, American: "water vapor"): the
    one source of every gas name in the text the generator writes (GEN.415)."""
    return gases.gas(formula).name


def species_name(species):
    """The name of one of the ten stored species (`SPECIES`, lower case: `"h2o"`) in words, as `gas_name` gives it."""
    return gases.stored(species).name


def describe(partials_pa):
    """The free-text atmosphere, e.g. "a mix of nitrogen, oxygen, and argon, with traces of water vapor and carbon dioxide", from the
    partial pressures (gas-table formula to Pa) of the air; "None" for an empty one."""
    total_pa = sum(partials_pa.values())
    if total_pa <= 0.0:
        return "None"
    total_kpa = total_pa / KPA_TO_PA
    parts = sorted(((value / total_pa, gas_name(formula)) for formula, value in partials_pa.items() if value > 0.0),
                   key=lambda part: -part[0])
    majors = [name for fraction, name in parts if fraction >= MAJOR_FRACTION][:4]
    traces = [name for fraction, name in parts if TRACE_FRACTION <= fraction < MAJOR_FRACTION][:2]

    ten = {gas: partials_pa.get(gases.stored(gas).formula, 0.0) / total_pa for gas in SPECIES}
    qualifiers = []
    if total_kpa < THIN_KPA:
        qualifiers.append("thin")
    elif total_kpa > DENSE_KPA:
        qualifiers.append("dense")
    if ten["h2o"] >= HUMID_H2O_FRACTION:
        qualifiers.append("humid")
    reduced = sum(ten[gas] for gas in ("h2", "co", "ch4", "h2s"))
    oxidized = sum(ten[gas] for gas in ("o2", "co2", "so2"))
    if ten["h2"] < 0.5 and reduced > oxidized:
        qualifiers.append("reducing")
    lead = "a " + ", ".join(qualifiers) + " " if qualifiers else "a "

    if len(majors) == 1:
        text = f"{lead}mix of mostly {majors[0]}"
    else:
        text = f"{lead}mix of {_join(majors)}"
    if traces:
        text += f", with traces of {_join(traces)}"
    return text


def _join(names):
    if len(names) <= 2:
        return " and ".join(names)
    return ", ".join(names[:-1]) + ", and " + names[-1]


def store_air(planet, partials_pa, redox_diw):
    """
    Sets `planet`'s `ATMOSPHERE_FIELDS`, its `GAS_FIELD` and its `atmosphere` text from the solved air's partial pressures (gas-table
    formula to Pa; empty for an airless body) and the mantle's redox (`None` for a giant).
    """
    planet.mantle_delta_iw = redox_diw
    planet.mantle_redox = None if redox_diw is None else redox_class(redox_diw)
    setattr(planet, GAS_FIELD, dict(partials_pa))
    planet.atmosphere = describe(partials_pa)


KPA_TO_PA = 1000.0


def gas_partials_pa(body):
    """`body`'s air as a dict of gas-table formula to partial pressure in Pa: every gas it stores, the balance gases included. `body`
    is a `Planet` or a column-named dict (a row with its gases, `db.store.with_gases`). Empty for an airless body or one saved
    before GEN.85."""
    air = body.get(GAS_FIELD) if isinstance(body, Mapping) else getattr(body, GAS_FIELD, None)
    return dict(air or {})


def partial_kpa(body, gas):
    """`body`'s partial pressure of one of the ten species (`SPECIES`, lower case), kPa; 0 for a gas it lacks."""
    return gas_partials_pa(body).get(gases.stored(gas).formula, 0.0) / KPA_TO_PA


def partial_pressures_kpa(body):
    """`body`'s partial pressures of the ten species as a dict gas -> kPa (zeros for a gas it lacks, and for a body saved before
    GEN.85)."""
    air = gas_partials_pa(body)
    return {gas: air.get(gases.stored(gas).formula, 0.0) / KPA_TO_PA for gas in SPECIES}


