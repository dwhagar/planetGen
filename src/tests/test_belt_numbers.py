"""GEN.1508: a belt is numbered among the belts of its own star, closest first, from 1 (never by its slot among the star's planets and belts)."""

from planetgen.db import belt_names, near
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
    rows = [{"id": 21, "star_system_id": 1, "star_id": None, "orbital_index": 3, "lagrange_point": None},
            {"id": 22, "star_system_id": 1, "star_id": None, "orbital_index": 7, "lagrange_point": None},
            {"id": 30, "star_system_id": 2, "star_id": None, "orbital_index": 5, "lagrange_point": None},
            {"id": 23, "star_system_id": 1, "star_id": None, "orbital_index": 4, "lagrange_point": "L4"}]       # a trojan region (GEN.1150): a name of its own, no number, and no gap in the belts' counting
    assert belt_names.numbers_for_systems(_Conn(rows), [1, 2]) == {21: 1, 22: 2, 30: 1}
    system_one = [row for row in rows if row["star_system_id"] == 1]
    assert [bodies.belt_number(row, system_one) for row in system_one if not row["lagrange_point"]] == [1, 2]


class _SystemConn:
    """Answers the four statements that name the belts of one system, from `belts` (each `id`, `orbital_index`, `lagrange_point`, `host`)."""

    def __init__(self, belts):
        self.belts = belts

    def execute(self, sql, params=()):
        if "p.name AS host_name" in sql:
            return _Rows([{"id": b["id"], "lagrange_point": b["lagrange_point"], "host_name": b["host"]} for b in self.belts if b["lagrange_point"] and b["id"] in params])
        if "orbital_index" in sql:
            return _Rows([{"id": b["id"], "star_system_id": 1, "star_id": None, "orbital_index": b["orbital_index"], "lagrange_point": b["lagrange_point"]} for b in self.belts])
        if "SELECT star_system_id" in sql:
            return _Rows([{"star_system_id": 1}])
        return _Rows([{"id": b["id"]} for b in self.belts])


def test_a_system_with_a_trojan_region_lists_its_belts_by_name_and_the_region_by_its_hosts_name():
    """NAV.79: `near` read a number for every belt of a system, and a trojan region takes none (GEN.1150), so it raised KeyError."""
    belts = [{"id": 21, "orbital_index": 3, "lagrange_point": None, "host": None},
             {"id": 23, "orbital_index": 4, "lagrange_point": "L4", "host": "Zorax IV"},
             {"id": 22, "orbital_index": 7, "lagrange_point": None, "host": None}]
    rows = near._bodies_of(_SystemConn(belts), 1, "Zorax", {"belt"})
    assert [(kind, row_id, name) for kind, row_id, name, _parent in rows] == [
        ("belt", 21, "Asteroid belt 1"), ("belt", 23, "Zorax IV L4 trojan region"), ("belt", 22, "Asteroid belt 2")]
