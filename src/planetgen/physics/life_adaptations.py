# planetgen/physics/life_adaptations.py

"""
The adaptations of life to cold, heat, dryness, salt and dim light (GEN.1725, docs/design/habitability/life-adaptations.md)
==========================================================================================================================

The data file `physics/data/life_adaptations.json` holds one row for each way life on Earth has widened where it can live (from the research of GEN.1724,
`reports/life-adaptations-climate-2026-11-14.md`): the climate it serves, the kind of life that carries it, the axis whose limit it widens, the limit at the tolerant
(grade 2) and extreme (grade 3) grades, the stage the life must have reached before it can have it, and the habitat it opens. This module loads and checks the file and
answers the questions the generator asks of it:

* `ordinary(axis)`: the limit of a life with no adaptation (one number for each axis, so no two rows disagree about it);
* `limit(axis, grade, ...)`: how far the limit goes at a grade, for the kinds of life and the class a body has and the age its life has reached (the widest value of every
  row that applies; the ordinary limit when none does);
* `habitat(id)`: a refuge kind the gate of GEN.1723 counts, with the rows that open it.

The axes are `t_cold` and `t_hot` (the local temperature of the water the life works in, K), `light_floor` (the least light for net growth, umol photons m-2 s-1),
`water_hours` (the hours of liquid water a year it needs), `water_activity`, `ground_water_depth` (how deep a water table its roots reach, m) and `insulation` (how much
warmer the life keeps the ground than the air, K). A row that moves no number (the snow algae's sunscreen) has no axis and only opens a habitat.

What reads it: `adaptation_traits.limit` (the traits cold and heat of GEN.876's draw, which store which grade a body's life has reached) and, once GEN.1723 has merged,
the refuge gate that widens a class's windows. Nothing in the generation reads a limit yet. A dormant row (`activity` is `dormant`: freeze tolerance, desiccation
tolerance, spores) lets life survive a limit but not work past it, so `limit` leaves it out unless asked. Every number is a stand-in of section 6.294 of
`docs/design/research/gap-filling-plan.md` with its tier and tweak level, carried by the row itself.
"""

import json
import os
from dataclasses import dataclass

_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "life_adaptations.json")

LIFE_TYPES = ("microbe", "oxygenic_microbe", "alga", "lichen_crust", "moss", "plant", "animal")
"""tuple: The kinds of life a row can be carried by, from the simplest up: any microbe, an oxygen-making microbe, an alga, a lichen or crust, a moss, a vascular plant, an animal."""
CLIMATES = ("IC", "FM", "SN", "CT", "WT", "HW", "HD", "XA", "SB", "DL", "any_cold", "any_hot")
"""tuple: The climate codes of the research (IC ice-covered water, FM frozen land with melt, SN snow and glacier, CT cold temperate, WT warm wet, HW hot wet, HD hot dry, XA hyperarid, SB brine, DL dim light)."""
ACTIVITIES = ("active", "dormant")
TIERS = tuple(f"T{n}" for n in range(1, 10))
TWEAKS = tuple(f"L{n}" for n in range(0, 5))
OXYGEN_CLASSES = ("1a", "1b", "1c")
"""tuple: The life classes that make oxygen, the ones the refuge term of GEN.1723 is built on."""


class AdaptationError(ValueError):
    """The data file is refused."""


@dataclass(frozen=True)
class Adaptation:
    """One row of the table."""

    id: str
    name: str
    climate: tuple
    life_types: tuple
    classes: tuple
    axis: object
    grade2: object
    grade3: object
    activity: str
    production_factor: object
    stage_needed: str
    opens_habitat: object
    effect_on_habitat: object
    source_rows: str
    tier: str
    tweak: str
    note: object

    def values(self, grade):
        """The limits this row gives at `grade` (2 or 3; grade 3 holds the grade 2 limit too), nearest the ordinary first: a list of 0 to 2 numbers."""
        found = []
        if grade >= 2 and self.grade2 is not None:
            found.append(self.grade2)
        if grade >= 3 and self.grade3 is not None:
            found.append(self.grade3)
        return found


def _load(path=_PATH):
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    axes, stages, habitats = data["axes"], data["stages"], {row["id"]: row for row in data["habitats"]}
    rows, seen = [], set()
    for raw in data["adaptations"]:
        where = f"adaptation {raw.get('id')!r}"
        if raw["id"] in seen:
            raise AdaptationError(f"{where} appears twice")
        seen.add(raw["id"])
        if raw["axis"] is not None and raw["axis"] not in axes:
            raise AdaptationError(f"{where}: the axis {raw['axis']!r} is not one of {sorted(axes)}")
        if raw["axis"] is None and (raw["grade2"] is not None or raw["grade3"] is not None):
            raise AdaptationError(f"{where}: a row with no axis cannot hold a limit")
        if raw["axis"] is not None and not any(raw[g] is not None for g in ("grade2", "grade3")) and raw["opens_habitat"] is None and not raw["effect_on_habitat"] and not raw["note"]:
            raise AdaptationError(f"{where}: an axis with no limit, habitat or effect says nothing")
        for kind in raw["life_types"]:
            if kind not in LIFE_TYPES:
                raise AdaptationError(f"{where}: {kind!r} is not a kind of life")
        for code in raw["climate"]:
            if code not in CLIMATES:
                raise AdaptationError(f"{where}: {code!r} is not a climate")
        if raw["activity"] not in ACTIVITIES:
            raise AdaptationError(f"{where}: the activity is 'active' or 'dormant'")
        if raw["stage_needed"] not in stages:
            raise AdaptationError(f"{where}: the stage {raw['stage_needed']!r} is not one of {sorted(stages)}")
        if raw["opens_habitat"] is not None and raw["opens_habitat"] not in habitats:
            raise AdaptationError(f"{where}: the habitat {raw['opens_habitat']!r} is not one of {sorted(habitats)}")
        if raw["tier"] not in TIERS or raw["tweak"] not in TWEAKS:
            raise AdaptationError(f"{where}: a tier T1 to T9 and a tweak L0 to L4 are needed")
        if raw["classes"] != ["all"] and any(c not in OXYGEN_CLASSES for c in raw["classes"]):
            raise AdaptationError(f"{where}: the classes are 'all' or oxygen-making ones {OXYGEN_CLASSES}")
        if raw["axis"] is not None:
            ordinary, wider = axes[raw["axis"]]["ordinary"], axes[raw["axis"]]["wider"]
            sign = -1.0 if wider == "lower" else 1.0
            previous = ordinary
            for grade in ("grade2", "grade3"):
                if raw[grade] is None:
                    continue
                if sign * (raw[grade] - previous) < 0.0:
                    raise AdaptationError(f"{where}: {grade} {raw[grade]} is not wider than {previous} (the axis widens {wider})")
                previous = raw[grade]
        row = {key: raw[key] for key in ("id", "name", "axis", "grade2", "grade3", "activity", "production_factor", "stage_needed", "opens_habitat",
                                         "effect_on_habitat", "source_rows", "tier", "tweak", "note")}
        rows.append(Adaptation(climate=tuple(raw["climate"]), life_types=tuple(raw["life_types"]), classes=tuple(raw["classes"]), **row))
    for axis, row in axes.items():
        if row["wider"] not in ("lower", "higher") or row["ordinary_tier"] not in TIERS or row["ordinary_tweak"] not in TWEAKS:
            raise AdaptationError(f"axis {axis!r} needs a direction (lower or higher), a tier and a tweak")
    return axes, stages, habitats, tuple(rows)


AXES, STAGES, HABITATS, ADAPTATIONS = _load()
"""The axes, the stages, the habitats by id and the rows of the data file, checked."""
BY_ID = {row.id: row for row in ADAPTATIONS}


def ordinary(axis):
    """The limit of a life with no adaptation on `axis`."""
    return AXES[axis]["ordinary"]


def stress_axis(stress):
    """The axis a trait stress moves (`adaptation_traits`: `cold` moves `t_cold`, `heat` moves `t_hot`), or `None` for a stress that moves no axis here."""
    for axis, row in AXES.items():
        if row.get("stress") == stress:
            return axis
    return None


def stage_delay_gyr(stage):
    """The earliest a life that began at time 0 has reached `stage`, Gyr (0 for the first life); `None` when the stage has no delay of its own (the animal stage follows the life timeline)."""
    delay = STAGES[stage]["delay_gyr"]
    return None if delay is None else delay[0]


def available(adaptation, life_age_gyr=None, animals=False):
    """Whether a life of age `life_age_gyr` (Gyr since it began; `None` for no age limit) can have `adaptation`: its stage reached (an animal stage needs `animals`)."""
    delay = stage_delay_gyr(adaptation.stage_needed)
    if delay is None:
        return bool(animals)
    return life_age_gyr is None or life_age_gyr >= delay


def rows_for(axis, life_types=None, class_id=None, life_age_gyr=None, animals=False, dormant=False):
    """The rows that move `axis` for a body: the kinds of life it has (any, when `life_types` is `None`), the class that lives on it (any, when `class_id` is `None`), the age its life
    has reached and whether it has animals; dormant rows only when `dormant`."""
    found = []
    for row in ADAPTATIONS:
        if row.axis != axis or (row.activity == "dormant" and not dormant):
            continue
        if life_types is not None and not set(life_types) & set(row.life_types):
            continue
        if class_id is not None and row.classes != ("all",) and class_id not in row.classes:
            continue
        if not available(row, life_age_gyr, animals):
            continue
        found.append(row)
    return found


def limit(axis, grade, life_types=None, class_id=None, life_age_gyr=None, animals=False, dormant=False):
    """
    How far the limit of `axis` goes for a body whose life has reached `grade` (1 ordinary, 2 tolerant, 3 extreme): the widest value any row that applies gives at that grade
    or below (`rows_for`), or the ordinary limit when none does or the grade is 1.
    """
    if grade not in (1, 2, 3):
        raise ValueError(f"the grade is 1, 2 or 3, not {grade!r}")
    base, wider = ordinary(axis), AXES[axis]["wider"]
    if grade == 1:
        return base
    values = [value for row in rows_for(axis, life_types, class_id, life_age_gyr, animals, dormant) for value in row.values(grade)]
    if not values:
        return base
    return min(values) if wider == "lower" else max(values)


def widest(axis, **where):
    """The widest limit any row on `axis` reaches (grade 3) for a body described as in `limit`; the end the draw's full stress is set at."""
    return limit(axis, 3, **where)


def habitat(habitat_id):
    """The habitat row `habitat_id` of the file (the refuge kinds the gate of GEN.1723 counts)."""
    return HABITATS[habitat_id]

