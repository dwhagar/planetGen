# planetgen/physics/stored_air.py

"""
The stored air history as the life-air API gives it (API.51; GEN.1745 slice 4b stores it, `physics/air_record.py` packs it)
==========================================================================================================================

`GET /api/planets/<id>/life-air` answered from the air as the body holds it now, worked backwards on request (`life_air.air_history`: oxygen, carbon dioxide and methane only). A body that has the
history of its air stored (the `air_history` column, a packed `air_record.Record`) answers from that instead: every gas the record keeps (any gas of `air_record.PARTIAL_FLOOR_PA` or more at some age),
the surface temperature and the total pressure at each recorded age, the mass extinctions that start a step and the life ages between them. The answer keeps the shape of the one it replaces (`age_gyr`,
`producers_gyr`, `rose_gyr`, `classes`, `gases`, `stages`, `speculative`), so the chart and every client that read it go on reading it, and adds `source`, `temperature_k`, `pressure_kpa`, `events` and
`life_ages`. What life changed in a gas, class by class (`by_class`), is not in the record: it comes from the rebuilt history where that has the gas, else it is empty.
"""

from planetgen.physics import air_record, gases

SOURCE_STORED = "stored"
"""str: The `source` of an answer made from the stored history (the rebuilt one says `rebuilt`)."""
SOURCE_REBUILT = "rebuilt"
"""str: The `source` of an answer worked backwards from the air now."""

PASCAL_PER_KPA = 1000.0
"""float: Pa in a kPa (the API gives kPa)."""


def _name(formula):
    """The gas's name as the chart gives it, with a capital (the table spells them in lower case), or the formula for a gas the table does not have."""
    if not gases.has(formula):
        return formula
    name = gases.gas(formula).name
    return name[:1].upper() + name[1:]


def history_of(record, star_age_gyr, stage_ages, rebuilt=None):
    """
    The life-air answer of a body from its stored `record` (`air_record.Record`), for a star `star_age_gyr` Gyr old. `stage_ages` are the milestone ages by key (`evolution.MILESTONE_KEYS`; `None` for one never
    reached), and `rebuilt` the answer of `life_air.air_history` for the same body or `None`: it supplies the classes, the producers' and the oxygen rise's ages, what each class made of a gas (`by_class`)
    and the less-sure note, and is the same body's life as the record knew it.

    Returns:
        dict: As `life_air.air_history` (see there) with `gases` holding every gas of the record, the biggest first (each `steps` is `[[age_gyr, kPa], ...]`, 0 where the gas is under the floor), and
            `source`, `temperature_k` (`[[age_gyr, K], ...]`), `pressure_kpa` (`[[age_gyr, kPa], ...]`), `events` (`{"age_gyr", "cause", "kill_share"}`) and `life_ages`
            (`{"number", "start_gyr", "end_gyr", "cause", "kill_share"}`, `end_gyr` `None` for the age still going).
    """
    ages = [step.age_gyr for step in record.steps]
    formulas = sorted({formula for step in record.steps for formula in step.partials_pa})
    made_by = {gas["formula"]: gas.get("by_class") or [] for gas in (rebuilt or {}).get("gases", ())}
    series = [{"formula": formula, "name": _name(formula),
               "steps": [[age, round(step.partials_pa.get(formula, 0.0) / PASCAL_PER_KPA, 6)] for age, step in zip(ages, record.steps)],
               "by_class": made_by.get(formula, [])} for formula in formulas]
    series.sort(key=lambda gas: (-max(value for _age, value in gas["steps"]), gas["formula"]))
    base = rebuilt or {"producers_gyr": None, "rose_gyr": None, "classes": [], "speculative": None}
    return {"age_gyr": star_age_gyr, "producers_gyr": base["producers_gyr"], "rose_gyr": base["rose_gyr"], "classes": base["classes"], "gases": series,
            "stages": sorted(({"key": key, "age_gyr": age} for key, age in stage_ages.items() if age is not None and age <= star_age_gyr), key=lambda stage: stage["age_gyr"]),
            "speculative": base["speculative"], "source": SOURCE_STORED,
            "temperature_k": [[step.age_gyr, step.temperature_k] for step in record.steps if step.temperature_k is not None],
            "pressure_kpa": [[step.age_gyr, round(step.pressure_pa / PASCAL_PER_KPA, 6)] for step in record.steps],
            "events": [{"age_gyr": step.age_gyr, "cause": step.event[0], "kill_share": step.event[1]} for step in record.steps if step.event is not None],
            "life_ages": [{"number": age.number, "start_gyr": age.start_gyr, "end_gyr": age.end_gyr, "cause": age.cause, "kill_share": age.kill_share} for age in record.ages]}


def from_text(text, star_age_gyr, stage_ages, rebuilt=None):
    """`history_of` the record packed in `text` (the `air_history` column), or `None` for no text."""
    if not text:
        return None
    return history_of(air_record.unpack(text), star_age_gyr, stage_ages, rebuilt)
