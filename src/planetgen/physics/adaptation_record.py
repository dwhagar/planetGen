# planetgen/physics/adaptation_record.py

"""
What the life of a finished body had to become to live where it does (GEN.1808; Boss, 5:54 PM EST 2026-11-20: "apply to ALL life form types but adding information to each about the adaptations that were necessary to make it work")
=====================================================================================================================================================================================================================

The refuge gate of GEN.1723 (`surface_refuge`) lets life live in a melt film, a brine seep, a spring or under a lid of ice, in a narrow range of temperature, light and water; the table of GEN.1725
(`life_adaptations`) says which adaptations of the life on Earth widen each limit and how far. This module puts the two together for **one finished body**: it reads the conditions of the places the
refuge counts (the coldest and the hottest wet lit step, the least light under the ice, the hours of liquid water (the lit hours of the refuge over the share of the time the star is up), the depth of the water table the springs rise from, the water activity of a brine), asks
for each what the life must be able to do (the demand), and finds the lowest grade of the table that reaches it (`life_adaptations.limit`: 1 is the ordinary life, 2 the tolerant, 3 the extreme grade),
for the kinds of life and the class the body has and the age its life has reached, and the rows of the table that give it.

* *A need* is one axis whose ordinary limit the life must pass. `grade` is the lowest grade that reaches the demand, `None` when even the extreme grade of every row that applies falls short (the
  refuge gate counts the place; the record says the life there is past the table). `rows` are the ids of the rows whose own limit at that grade reaches the demand.
* *A habitat row* has no axis (a pigment that works in no light, a sunscreen, a crust that wakes on water): where the refuge counted a place of the habitat it opens, the rows that open it and apply to
  the class are listed, with no demand and no grade.
* The record reads **only** the finished body; it moves no number of it (nothing in the generation reads it yet; GEN.1809 stores it). Classes without a refuge (`1a`, `1b`, `1c` and `2` are the
  photosynthetic ones the gate of GEN.1723 serves, the rest follow in GEN.1810 and GEN.1811) get no record.
* The light of the table is in micromoles of photons a square metre a second and the gate's in watts a square metre: they are tied by the gate's own floor, 2e-3 W/m2 being 0.01 umol (RW2), so one watt is
  five micromoles here (`UMOL_PER_W`).

Every number that is not read from the table is a stand-in of section 6.318 of `docs/design/research/gap-filling-plan.md` (GEN.725 rule).
"""

from dataclasses import dataclass

from planetgen.physics import life_adaptations, surface_refuge

CLASSES = ("1a", "1b", "1c", "2")
"""tuple: The life classes a record is made for: the photosynthetic ones the refuge gate of GEN.1723 serves (GEN.1808, first slice)."""
UMOL_PER_W = 5.0
"""float: Micromoles of photons a square metre a second in one watt a square metre of usable light: the gate's floor, 2e-3 W/m2, is 0.01 umol (`surface_refuge.LIGHT_FLOOR_W_M2`, RW2). [Stand-in 6.318 AE1, tier T6, tweak L0.]"""
FIELD = "life_adaptations"
"""str: The attribute of a planet or moon that holds its record, a list of the dicts `entries_of` makes (`None` for a body with none)."""


@dataclass(frozen=True)
class Need:
    """
    One thing the life of a body must be able to do.

    Attributes:
        axis (str or None): The axis of `life_adaptations.AXES` the demand is on; `None` for a habitat row with no axis.
        habitat (str): The habitat of `life_adaptations.HABITATS` the demand comes from.
        demand (float or None): What the places demand on the axis, in the unit of the axis (a temperature K, a light in umol m-2 s-1, the hours of liquid water a year, a depth m, a water activity).
        grade (int or None): The lowest grade (2 or 3) of the table that reaches `demand`; `None` when none does (or for a habitat row).
        rows (tuple): The ids of the rows that give it (or, for a habitat row, the row).
    """

    axis: str | None
    habitat: str
    demand: float | None
    grade: int | None
    rows: tuple


def demands(refuge, values):
    """`[(axis, habitat, demand)]` the places of `refuge` make on the life, from the conditions of the body (`values`, the column-named values `habitability_world.values_from_body` gives): the cold and the hot
    limit at the temperatures of the wet lit steps, the light floor under the ice, the hours of liquid water of a place that is not an open sea, the water table's depth for springs and the water activity of
    a brine. An open sea (no place of its own) makes none."""
    found = []
    habitats = refuge.habitats
    if refuge.coldest_k is not None:
        found.append(("t_cold", _first(habitats, ("brine", "melt_film", "seep")), refuge.coldest_k))
    if refuge.hottest_k is not None:
        found.append(("t_hot", _first(habitats, ("brine", "melt_film", "seep")), refuge.hottest_k))
    if "water_under_ice" in habitats and refuge.under_ice_w_m2 is not None:
        found.append(("light_floor", "water_under_ice", refuge.under_ice_w_m2 * UMOL_PER_W))
    if refuge.hours and set(habitats) & {"melt_film", "brine", "seep"}:
        lit = values.get("day_share_mean")
        lit = surface_refuge.DEFAULT_LIT_SHARE if lit is None else min(max(float(lit), 1.0e-3), 1.0)
        found.append(("water_hours", _first(habitats, ("melt_film", "brine", "seep")), min(refuge.hours / lit, surface_refuge.YEAR_HOURS)))        # the hours the water is liquid, lit or not
    if "seep" in habitats and values.get("aquifer_depth_km") is not None:
        found.append(("ground_water_depth", "seep", values["aquifer_depth_km"] * 1000.0))
    if "brine" in habitats and values.get("water_activity"):
        found.append(("water_activity", "brine", values["water_activity"]))
    return found


def _first(habitats, wanted):
    """The first of `wanted` among `habitats` (the habitat a temperature or hours demand is put to)."""
    return next((habitat for habitat in wanted if habitat in habitats), wanted[0])


def _reaches(axis, limit, demand):
    """Whether `limit` goes as far as `demand` on `axis` (the lower is wider on some axes, the higher on others)."""
    wider = life_adaptations.AXES[axis]["wider"]
    return limit <= demand if wider == "lower" else limit >= demand


def needs(refuge, values, class_id, life_age_gyr=None, animals=False):
    """
    The needs of the life of class `class_id` on a body with `refuge` (`surface_refuge.refuge_of`) and `values`, at the life age `life_age_gyr` Gyr (`None` for no limit): one `Need` for each demand
    that the ordinary limit of its axis does not reach, then one for each habitat row that opens a habitat the refuge counted. `[]` for a class with no record, an open sea and a body with no refuge.
    """
    if class_id not in CLASSES or not refuge.habitats:
        return []
    found = []
    for axis, habitat, demand in demands(refuge, values):
        if _reaches(axis, life_adaptations.ordinary(axis), demand):
            continue
        grade = next((g for g in (2, 3) if _reaches(axis, life_adaptations.limit(axis, g, class_id=class_id, life_age_gyr=life_age_gyr, animals=animals), demand)), None)
        rows = ()
        if grade is not None:
            rows = tuple(row.id for row in life_adaptations.rows_for(axis, class_id=class_id, life_age_gyr=life_age_gyr, animals=animals)
                         if any(_reaches(axis, value, demand) for value in row.values(grade)))
        found.append(Need(axis, habitat, demand, grade, rows))
    seen = {id_ for need in found for id_ in need.rows}
    for row in life_adaptations.ADAPTATIONS:
        if row.axis is None and row.opens_habitat in refuge.habitats and row.activity == "active" and row.id not in seen:
            if row.classes != ("all",) and class_id not in row.classes:
                continue
            if not life_adaptations.available(row, life_age_gyr, animals):
                continue
            found.append(Need(None, row.opens_habitat, None, None, (row.id,)))
    return found


def entries_of(found, class_id):
    """The record of `needs` of the class `class_id` as a list of plain dicts (`life_class`, `axis`, `habitat`, `demand`, `grade`, `rows`), or `None` for none: the form the body keeps."""
    if not found:
        return None
    return [{"life_class": class_id, "axis": need.axis, "habitat": need.habitat, "demand": need.demand, "grade": need.grade, "rows": list(need.rows)} for need in found]


def record_of(body, class_ids, life_age_gyr=None, animals=False):
    """
    The record of the finished `body`: the needs of each of `class_ids` (the classes that live on it, the dominant one first) that has a record, as `entries_of` gives them, one class after the other (Boss, 5:54 PM EST
    2026-11-20: all the life form types). The refuge and the values are read from the body (`habitability_world.values_from_body`, `surface_refuge.refuge_of`), once. `None` when no class has a record or the body has no refuge.
    """
    from planetgen.physics import habitability_world

    wanted = [class_id for class_id in class_ids if class_id in CLASSES]
    if not wanted or not (getattr(body, "surface_temperature", None) or 0.0):
        return None
    values = habitability_world.values_from_body(body)
    refuge = surface_refuge.refuge_of(values)
    entries = []
    for class_id in wanted:
        entries.extend(entries_of(needs(refuge, values, class_id, life_age_gyr, animals), class_id) or ())
    return entries or None
