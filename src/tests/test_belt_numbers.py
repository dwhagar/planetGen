"""GEN.1508: a belt is numbered among the belts of its own star, closest first, from 1 (never by its slot among the star's planets and belts)."""

from planetgen.db import belt_names
from planetgen.names import bodies


def _belt(row_id, index, star_id=None):
    return {"id": row_id, "orbital_index": index, "star_id": star_id}


def test_a_lone_belt_behind_four_planets_is_belt_1():
    belt = _belt(9, 4)                       # the fifth slot out
    assert bodies.belt_label(bodies.belt_number(belt, [belt])) == "Asteroid belt 1"


def test_belts_count_out_from_the_star_whatever_the_order_of_the_list():
    inner, middle, outer = _belt(7, 1), _belt(3, 4), _belt(5, 9)
    belts = [outer, inner, middle]
    assert [bodies.belt_number(b, belts) for b in (inner, middle, outer)] == [1, 2, 3]


def test_each_star_of_a_wide_pair_counts_its_own_belts():
    a1, a2, b1 = _belt(1, 2, star_id=10), _belt(2, 5, star_id=10), _belt(3, 3, star_id=11)
    belts = [a1, a2, b1]
    assert [bodies.belt_number(b, belts) for b in belts] == [1, 2, 1]


class _Rows:
    def __init__(self, rows):
        self.rows = rows

    def fetchall(self):
        return self.rows


class _Conn:
    def __init__(self, rows):
        self.rows = rows

    def execute(self, sql, params=()):
        return _Rows(self.rows)


def test_the_database_numbers_agree_with_the_pure_count():
    rows = [{"id": 21, "star_system_id": 1, "star_id": None, "orbital_index": 3},
            {"id": 22, "star_system_id": 1, "star_id": None, "orbital_index": 7},
            {"id": 30, "star_system_id": 2, "star_id": None, "orbital_index": 5}]
    assert belt_names.numbers_for_systems(_Conn(rows), [1, 2]) == {21: 1, 22: 2, 30: 1}
