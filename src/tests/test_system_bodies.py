"""GEN.1066: the system ledger and the barycentre that includes the belts."""

import math

import pytest

from planetgen.physics import system_bodies as sb
from planetgen.physics.system_bodies import Body, SystemBodiesError

SUN = 2.0e30
JUPITERISH = 1.9e27


def _star(mass=SUN):
    return Body("s", None, "star", mass, 7e5, 1400.0, status="star", tier="typed")


def _system():
    rows = [
        _star(),
        Body("p1", "s", "planet", JUPITERISH, 7e4, 1300.0, 5.0, 0.05, 1.0, "planet", "typed", (5.0, 0.0, 0.0)),
        Body("m1", "p1", "moon", 1e23, 2e3, 3000.0, 0.003, 0.0, 0.0, "moon", "typed", (0.0, 0.003, 0.0)),
        Body("b", "s", "belt", 0.0, status="belt", tier="typed"),
    ]
    return rows


def _ring_of_members(count, mass_each, radius_au, remainder):
    rows = []
    for k in range(count):
        angle = 2 * math.pi * k / count
        rows.append(Body(f"bm{k}", "b", "belt_member", mass_each, 100.0, 2000.0, radius_au, 0.0, 0.0, "minor body", "L1",
                         (radius_au * math.cos(angle), radius_au * math.sin(angle), 0.0)))
    rows.append(Body("bx", "b", "extended", remainder, status="extended", tier="L3"))
    return rows


def test_a_valid_system_orders_every_parent_before_its_children():
    order = sb.check(_system())
    assert order.index("s") < order.index("p1") < order.index("m1")
    assert set(order) == {"s", "p1", "m1", "b"}


@pytest.mark.parametrize("change,message", [
    (lambda r: r + [r[1]], "duplicate"),
    (lambda r: r + [Body("x", "s", "wibble", 1.0)], "unknown kind"),
    (lambda r: r + [Body("x", "s", "planet", -1.0)], "mass"),
    (lambda r: r + [Body("x", "s", "planet", float("nan"))], "mass"),
    (lambda r: r + [Body("x", "nowhere", "planet", 1.0)], "not in the ledger"),
    (lambda r: r + [Body("x", "m1", "planet", 1.0)], "cannot orbit"),
    (lambda r: r + [Body("x", "s", "planet", 1.0, position_au=(1.0, 2.0))], "three finite"),
    (lambda r: r + [Body("x", "s", "planet", 1.0, status="mystery")], "unknown status"),
    (lambda r: r + [Body("x", "s", "planet", 1.0, tier="L9")], "unknown tier"),
    (lambda r: r + [Body("x", None, "star", 1.0)], "exactly one root"),
    (lambda r: [Body("b0", None, "belt", 0.0)], "must be one of"),
    (lambda r: [Body("b", "s", "belt", 5.0)] + r[:1], "container"),
])
def test_a_bad_ledger_is_refused_with_the_reason(change, message):
    rows = _system()
    if message == "container":
        rows = [rows[0], Body("b", "s", "belt", 5.0)]
    else:
        rows = change(rows)
    with pytest.raises(SystemBodiesError, match=message):
        sb.check(rows)


def test_a_cycle_is_refused():
    rows = [_star(), Body("p", "q", "planet", 1.0), Body("q", "p", "planet", 1.0)]
    with pytest.raises(SystemBodiesError):
        sb.check(rows)


def test_a_moon_sits_at_its_planet_plus_its_own_offset():
    places = sb.absolute_positions_au(_system())
    assert places["s"] == (0.0, 0.0, 0.0)
    assert places["p1"] == (5.0, 0.0, 0.0)
    assert places["m1"] == pytest.approx((5.0, 0.003, 0.0))
    assert places["b"] == (0.0, 0.0, 0.0)


def test_a_subtree_counts_as_the_sum_of_its_parts_and_the_total_counts_each_body_once():
    rows = _system()
    assert sb.subtree_mass_kg(rows, "p1") == pytest.approx(JUPITERISH + 1e23)
    assert sb.total_mass_kg(rows) == pytest.approx(SUN + JUPITERISH + 1e23)
    assert sb.subtree_mass_kg(rows, "s") == pytest.approx(sb.total_mass_kg(rows))


def test_the_barycentre_of_one_planet_is_the_two_body_formula():
    rows = _system()[:2]
    expected = JUPITERISH * 5.0 / (SUN + JUPITERISH)
    assert sb.barycentre_au(rows) == pytest.approx((expected, 0.0, 0.0), rel=1e-12)
    assert sb.central_offset_au(rows) == pytest.approx((-expected, 0.0, 0.0), rel=1e-12)


def test_a_moon_moves_the_barycentre_by_its_share_of_the_total_mass():
    rows = _system()
    expected_x = (JUPITERISH * 5.0 + 1e23 * 5.0) / (SUN + JUPITERISH + 1e23)
    assert sb.barycentre_au(rows)[0] == pytest.approx(expected_x, rel=1e-12)
    assert sb.barycentre_au(rows)[1] == pytest.approx(1e23 * 0.003 / (SUN + JUPITERISH + 1e23), rel=1e-9)


def test_a_round_centred_belt_does_not_move_the_star():
    rows = _system() + _ring_of_members(12, 1e20, 2.7, 3e22)
    assert sb.belt_wobble_au(rows) < 1e-9 * 2.7
    # and the belt's mass is in the total
    assert sb.total_mass_kg(rows) == pytest.approx(SUN + JUPITERISH + 1e23 + 12e20 + 3e22)


def test_one_heavy_member_on_one_side_moves_the_star_by_the_two_body_amount():
    heavy = 1e24
    rows = [_star(), Body("b", "s", "belt", 0.0), Body("bm", "b", "belt_member", heavy, 1000.0, 2000.0, 2.7, 0.0, 0.0, "dwarf planet", "L1", (2.7, 0.0, 0.0))]
    assert sb.belt_wobble_au(rows) == pytest.approx(heavy * 2.7 / (SUN + heavy), rel=1e-12)


def test_the_extended_remainder_of_a_belt_sits_at_its_centre_and_does_not_move_the_star_alone():
    rows = [_star(), Body("b", "s", "belt", 0.0), Body("bx", "b", "extended", 4e23, tier="L3")]
    assert sb.belt_wobble_au(rows) == 0.0


def test_the_belt_residual_is_zero_when_the_rows_add_up_and_signed_when_they_do_not():
    rows = _system() + _ring_of_members(4, 1e20, 2.7, 6e20)
    assert sb.belt_mass_residual(rows, "b", 1e21) == pytest.approx(0.0, abs=1e-12)
    assert sb.belt_mass_residual(rows, "b", 5e20) == pytest.approx(1.0)
    assert sb.belt_mass_residual(rows, "b", 2e21) == pytest.approx(-0.5)
    with pytest.raises(SystemBodiesError):
        sb.belt_mass_residual(rows, "p1", 1e21)
    with pytest.raises(SystemBodiesError):
        sb.belt_mass_residual(rows, "b", 0.0)


def test_no_mass_to_take_a_barycentre_of_is_an_error():
    rows = [_star(0.0)]
    with pytest.raises(SystemBodiesError):
        sb.barycentre_au(rows)


def test_elements_to_a_position():
    assert sb.relative_position_au(2.0, 0.0, 0.0, 0.0, 0.0, 0.0) == pytest.approx((2.0, 0.0, 0.0))
    assert sb.relative_position_au(2.0, 0.0, 0.0, 0.0, 0.0, 90.0) == pytest.approx((0.0, 2.0, 0.0), abs=1e-12)
    assert sb.relative_position_au(2.0, 0.0, 90.0, 0.0, 0.0, 90.0) == pytest.approx((0.0, 0.0, 2.0), abs=1e-12)
    assert sb.relative_position_au(2.0, 0.25, 0.0, 0.0, 0.0, 0.0) == pytest.approx((1.5, 0.0, 0.0))
    assert sb.relative_position_au(2.0, 0.25, 0.0, 0.0, 0.0, 180.0) == pytest.approx((-2.5, 0.0, 0.0), abs=1e-9)
    # the node turns the whole orbit about z
    assert sb.relative_position_au(2.0, 0.0, 0.0, 90.0, 0.0, 0.0) == pytest.approx((0.0, 2.0, 0.0), abs=1e-12)
    # the distance is preserved by any orientation
    place = sb.relative_position_au(3.0, 0.3, 37.0, 112.0, 250.0, 71.0)
    ecc = sb.solve_eccentric_anomaly(math.radians(71.0), 0.3)
    assert math.hypot(*place) == pytest.approx(3.0 * (1 - 0.3 * math.cos(ecc)), rel=1e-12)


@pytest.mark.parametrize("a,e", [(0.0, 0.1), (-1.0, 0.1), (1.0, 1.0), (1.0, -0.1), (float("nan"), 0.1)])
def test_bad_elements_are_refused(a, e):
    with pytest.raises(SystemBodiesError):
        sb.relative_position_au(a, e, 0.0, 0.0, 0.0, 0.0)
