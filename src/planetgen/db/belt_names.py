"""
GEN.1508: the number a belt is shown with. A belt has no stored name; it is numbered among the belts of its own star by its place out from the star, closest first, from 1
(only the belts that exist count, so there are no gaps, and a lone belt is `Asteroid belt 1`). `orbital_index` is the belt's slot among all the star's bodies (planets and belts), so it is
used only to put the belts of a star in order, never as the number. A wide pair's belt belongs to the star it orbits (`star_id`); every other belt of a system shares the one count.
"""

from planetgen.names import bodies

label = bodies.belt_label
number_of = bodies.belt_number


def numbers_for_belts(conn, belt_ids):
    """`{belt id: number}` for the belts `belt_ids`, counted among all the belts of each one's system."""
    belt_ids = list(belt_ids)
    if not belt_ids:
        return {}
    marks = ", ".join("?" * len(belt_ids))
    systems = sorted({row["star_system_id"] for row in conn.execute(f"SELECT star_system_id FROM asteroid_belts WHERE id IN ({marks})", tuple(belt_ids)).fetchall()})
    return {belt_id: number for belt_id, number in numbers_for_systems(conn, systems).items() if belt_id in set(belt_ids)}


def numbers_for_systems(conn, system_ids):
    """`{belt id: number}` for every belt of the systems `system_ids`."""
    system_ids = list(system_ids)
    if not system_ids:
        return {}
    marks = ", ".join("?" * len(system_ids))
    rows = [dict(row) for row in conn.execute(
        f"SELECT id, star_system_id, star_id, orbital_index FROM asteroid_belts WHERE star_system_id IN ({marks}) ORDER BY star_system_id, orbital_index, id", tuple(system_ids)).fetchall()]
    count, found = {}, {}
    for row in rows:
        key = (row["star_system_id"], row["star_id"])
        count[key] = count.get(key, 0) + 1
        found[row["id"]] = count[key]
    return found
