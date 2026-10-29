# planetgen/db/population_slot.py

"""
Population as a metric slot (MAP.318, regions chain step 3d; docs/design/habitability/population-and-politics.md).

The map ranks worlds for its inhabited label slice (MAP.317) and its charted-regions list (MAP.316) by head count, and a head count
is stored by nothing yet: the population pass stores species, civilizations, polities and system ownership, and the world
population simulation that stores the head counts is POP.30 (8.3). So the slot is built first and reads what it is given:

* `rank_key(world)`: a world with a population ranks by it, before every world without one; the worlds without one rank by their
  habitability score (so until head counts exist the inhabited slice is the habitable worlds, best first). Ties go to the name,
  then the id, so the order never depends on the order the worlds came in.
* `region_totals(conn)`: the head count of each region, `{(shell, band, wedge): people}`, holding only regions that have one. It is
  empty while nothing stores head counts, and POP.30 fills it by summing its stored head counts per 9-block.
* `population_of(world)`: the head count a world states, or `None` (a missing, zero or negative value is no population, which keeps
  the hidden-when-absent rule).

`region_stats.charted_regions` puts `population` on a region entry that has a total and leaves it off the others, so a list that
shows a population column has nothing to show it for until `region_totals` returns something.
"""


def population_of(world):
    """The head count `world` (a mapping) states, or `None` when it has none."""
    value = world.get("population")
    if value is None or not value > 0:                        # also false for NaN
        return None
    return float(value)


def rank_key(world):
    """
    A sort key that puts the first-ranked world first: worlds with a population by population (largest first), then the worlds
    without one by `habitability` (the overall score, largest first; a missing score is 0), ties by `name` and `id`.
    """
    people = population_of(world)
    score = world.get("habitability")
    return (0 if people is not None else 1, -(people or 0.0), -(float(score) if score is not None else 0.0),
            str(world.get("name") or ""), str(world.get("id") or ""))


def rank(worlds):
    """`worlds` (mappings with optional `population`, `habitability`, `name`, `id`) best-ranked first."""
    return sorted(worlds, key=rank_key)


def region_totals(conn):
    """
    The head count of each region that has one, `{(shell, band, wedge): people}` of the 9-block cell. Empty: no table stores head
    counts yet (POP.30 stores them, and this sums them per cell then).
    """
    return {}
