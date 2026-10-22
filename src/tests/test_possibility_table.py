# tests/test_possibility_table.py

"""
GEN.472: the gated possibility table and its loader (`planetgen.generation.possibility_table`): the class gate and its
condition predicates, a weight function, normalisation to 100 percent, the empty-gate and all-zero errors, the refusal of a
row without a source, the open and stand-in report, the Dirichlet fractions, the shared components table,
the registry, and determinism at any worker count. Also `draw.gamma` and `draw.dirichlet`, which the fractions use.
"""

import math
import statistics

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.generation.possibility_table import (ANY_CLASS, OPEN, OPEN_STATE, PARTLY_VERIFIED, STAND_IN, Condition,
                                                    PossibilityError, Registry, Row, Subject, Table, TableError)
from planetgen.util import draw


def row(name, classes=("rocky",), weight=1.0, **kw):
    kw.setdefault("source", "Test et al. 2020")
    return Row(name, classes, weight, **kw)


def open_row(name, classes=("rocky",), **kw):
    return Row(name, classes, OPEN, "open", state=OPEN_STATE, **kw)


def probs(table, body_class="rocky", **quantities):
    return {r.name: p for r, p in table.probabilities(Subject(body_class, quantities))}


# --- draw.gamma and draw.dirichlet ---------------------------------------------------------------------------------


def test_gamma_draws_have_the_mean_and_variance_of_the_shape():
    for shape in (0.4, 1.0, 2.5, 9.0):
        stream = draw.Stream(f"gamma:{shape}")
        values = [stream.gamma(shape) for _ in range(20000)]
        assert statistics.fmean(values) == pytest.approx(shape, rel=0.05)
        assert statistics.pvariance(values) == pytest.approx(shape, rel=0.12)
        assert min(values) >= 0.0


def test_a_dirichlet_draw_adds_up_to_one_and_its_mean_follows_the_alphas():
    stream = draw.Stream("dirichlet")
    alphas = [6.0, 3.0, 1.0]
    draws = [stream.dirichlet(alphas) for _ in range(6000)]
    for d in draws:
        assert math.fsum(d) == pytest.approx(1.0, abs=1e-12)
        assert all(x >= 0 for x in d)
    for i, alpha in enumerate(alphas):
        assert statistics.fmean(d[i] for d in draws) == pytest.approx(alpha / sum(alphas), abs=0.01)


def test_gamma_and_dirichlet_are_the_same_from_the_same_seed_and_refuse_bad_shapes():
    assert draw.Stream(9).gamma(2.0) == draw.Stream(9).gamma(2.0)
    assert draw.Stream(9).dirichlet([1, 2, 3]) == draw.Stream(9).dirichlet([1, 2, 3])
    assert draw.Stream(9).dirichlet([1, 2, 3]) != draw.Stream(10).dirichlet([1, 2, 3])
    with pytest.raises(ValueError):
        draw.Stream(1).gamma(0.0)
    with pytest.raises(ValueError):
        draw.Stream(1).dirichlet([])
    with draw.bound(4):
        assert sum(draw.dirichlet([1.0, 1.0])) == pytest.approx(1.0)


def test_a_dirichlet_with_tiny_alphas_still_returns_a_proportion_vector():
    out = draw.Stream(3).dirichlet([1e-9, 1e-9, 1e-12])
    assert math.fsum(out) == pytest.approx(1.0, abs=1e-12)


# --- The gate -----------------------------------------------------------------------------------------------------


def test_a_row_is_drawn_only_for_the_classes_it_supports():
    table = Table("t", 1, [row("a", ("rocky",)), row("b", ("rocky", "icy")), row("c", ("gas-giant",))])
    assert probs(table, "rocky") == {"a": 0.5, "b": 0.5}
    assert probs(table, "icy") == {"b": 1.0}
    assert probs(table, "gas-giant") == {"c": 1.0}
    assert table.classes() == ["gas-giant", "icy", "rocky"]


def test_any_class_supports_every_class():
    table = Table("t", 1, [row("a", (ANY_CLASS,)), row("b", ("icy",))])
    assert probs(table, "icy") == {"a": 0.5, "b": 0.5}
    assert probs(table, "never-heard-of") == {"a": 1.0}
    assert table.classes() == ["icy"]


def test_condition_predicates_gate_on_stored_quantities():
    table = Table("t", 1, [
        row("ice", ("rocky",), conditions=(Condition("formation_distance_au", minimum=2.7),)),
        row("rock", ("rocky",), conditions=(Condition("formation_distance_au", maximum=3.0),)),
        row("oxidised", ("rocky",), conditions=(Condition("redox", equals="oxidised"),))])
    assert probs(table, formation_distance_au=1.0, redox="reduced") == {"rock": 1.0}
    assert probs(table, formation_distance_au=2.8, redox="reduced") == {"ice": 0.5, "rock": 0.5}
    assert probs(table, formation_distance_au=5.0, redox="oxidised") == {"ice": 0.5, "oxidised": 0.5}
    assert probs(table, formation_distance_au=2.7, redox="x")["ice"] == 0.5     # a bound is inclusive


def test_a_condition_on_a_quantity_the_subject_lacks_is_an_error_and_not_a_pass():
    table = Table("t", 1, [row("ice", conditions=(Condition("snow_line_au", maximum=3.0),))])
    with pytest.raises(KeyError, match="snow_line_au"):
        table.probabilities(Subject("rocky", {}))


# --- The weights --------------------------------------------------------------------------------------------------


def test_the_weights_of_the_rows_that_pass_are_normalised_to_100_percent():
    table = Table("t", 1, [row("a", weight=2.0), row("b", weight=1.0), row("c", weight=0.5), row("d", ("icy",))])
    p = probs(table)
    assert p == pytest.approx({"a": 2 / 3.5, "b": 1 / 3.5, "c": 0.5 / 3.5}, abs=1e-15)
    assert math.fsum(p.values()) == pytest.approx(1.0, abs=1e-12)


def test_a_weight_may_be_a_function_of_the_subject():
    table = Table("t", 1, [row("metal-rich", weight=lambda s: 10 ** s["fe_h"]), row("metal-poor", weight=1.0)])
    assert probs(table, fe_h=0.0) == {"metal-rich": 0.5, "metal-poor": 0.5}
    assert probs(table, fe_h=1.0)["metal-rich"] == pytest.approx(10 / 11)


def test_a_row_with_a_weight_of_zero_is_dropped():
    table = Table("t", 1, [row("a", weight=0.0), row("b", weight=lambda s: s.get("share", 0.0)), row("c")])
    assert probs(table) == {"c": 1.0}
    assert probs(table, share=1.0) == {"b": 0.5, "c": 0.5}


def test_a_weight_function_that_returns_a_bad_number_is_an_error():
    for bad in (-1.0, math.nan, math.inf, "x", True):
        table = Table("t", 1, [row("a", weight=lambda s, bad=bad: bad)])
        with pytest.raises(TableError, match="weight"):
            table.probabilities(Subject("rocky"))


def test_open_rows_share_what_the_sourced_weights_leave_equally():
    table = Table("t", 1, [row("a", weight=0.5), open_row("x"), open_row("y")])
    assert probs(table) == pytest.approx({"a": 0.5, "x": 0.25, "y": 0.25})
    only_open = Table("t", 1, [open_row("x"), open_row("y"), open_row("z")])
    assert probs(only_open) == pytest.approx({"x": 1 / 3, "y": 1 / 3, "z": 1 / 3})


def test_open_rows_of_another_gate_do_not_take_a_share():
    table = Table("t", 1, [row("a", ("rocky",), 0.5), open_row("x", ("icy",)), row("b", ("rocky",), 0.5)])
    assert probs(table, "rocky") == {"a": 0.5, "b": 0.5}
    assert probs(table, "icy") == {"x": 1.0}


def test_open_rows_get_nothing_when_the_sourced_weights_use_up_the_whole():
    table = Table("t", 1, [row("a", weight=1.0), open_row("x")])
    assert probs(table) == {"a": 1.0}


def test_open_rows_with_sourced_weights_over_one_are_an_error():
    table = Table("t", 1, [row("a", weight=0.8), row("b", weight=0.6), open_row("x")])
    with pytest.raises(PossibilityError, match="over 1"):
        table.probabilities(Subject("rocky"))


# --- The errors ---------------------------------------------------------------------------------------------------


def test_a_class_no_row_supports_is_an_error_that_names_the_class_and_the_table():
    table = Table("gases", 1, [row("a", ("rocky",))])
    with pytest.raises(PossibilityError) as caught:
        table.draw("icy-moon")
    assert caught.value.table == "gases" and caught.value.body_class == "icy-moon"
    assert "icy-moon" in str(caught.value) and "gases" in str(caught.value)


def test_a_table_whose_passing_weights_are_all_zero_is_an_error():
    table = Table("gases", 1, [row("a", weight=0.0), row("b", weight=lambda s: 0.0)])
    with pytest.raises(PossibilityError, match="every weight is zero"):
        table.draw("rocky")


def test_a_row_without_a_source_is_refused_and_a_row_has_no_licence_field():
    with pytest.raises(TableError, match="no source"):
        row("a", source="  ")
    with pytest.raises(TypeError):
        row("a", licence="CC BY 4.0")


def test_a_row_with_no_classes_or_no_name_is_refused():
    with pytest.raises(TableError, match="empty gate"):
        row("a", classes=())
    with pytest.raises(TableError, match="no name"):
        row(" ")


def test_the_open_state_and_the_open_weight_go_together():
    with pytest.raises(TableError, match="open"):
        Row("a", ("rocky",), OPEN, "Test 2020")                            # an open weight with a source and state
    with pytest.raises(TableError, match="open"):
        Row("a", ("rocky",), 0.3, "open", state=OPEN_STATE)
    with pytest.raises(TableError, match="open"):
        Row("a", ("rocky",), OPEN, "open")                          # the state defaults to verified


def test_a_weight_must_be_a_finite_number_of_zero_or_more():
    for bad in (-0.1, math.nan, math.inf, True, "heavy"):
        with pytest.raises(TableError, match="weight"):
            row("a", weight=bad)


def test_a_partly_verified_row_names_its_limit_and_a_stand_in_names_its_replacement():
    with pytest.raises(TableError, match="limit"):
        row("a", state=PARTLY_VERIFIED)
    assert row("a", state=PARTLY_VERIFIED, limit="valid 0.5 to 5 au").limit
    with pytest.raises(TableError, match="replacement"):
        row("a", state=STAND_IN)
    assert row("a", state=STAND_IN, replacement="GEN.473").replacement == "GEN.473"
    with pytest.raises(TableError, match="state"):
        row("a", state="probably fine")


def test_a_table_with_no_rows_a_repeated_row_or_a_bad_version_is_refused():
    with pytest.raises(TableError, match="no rows"):
        Table("t", 1, [])
    with pytest.raises(TableError, match="twice"):
        Table("t", 1, [row("a"), row("a", ("icy",))])
    for version in (0, -1, 1.5, True):
        with pytest.raises(TableError, match="version"):
            Table("t", version, [row("a")])


# --- The draw -----------------------------------------------------------------------------------------------------


def test_a_draw_follows_the_normalised_weights():
    table = Table("t", 1, [row("a", weight=6.0), row("b", weight=3.0), row("c", weight=1.0)])
    stream = draw.Stream("pick")
    counts = {"a": 0, "b": 0, "c": 0}
    with draw.bound(stream):
        for _ in range(20000):
            counts[table.draw("rocky").name] += 1
    for name, share in (("a", 0.6), ("b", 0.3), ("c", 0.1)):
        assert counts[name] / 20000 == pytest.approx(share, abs=0.015)


def test_a_draw_reports_the_probability_it_was_made_with():
    table = Table("t", 1, [row("a", weight=3.0), row("b", weight=1.0)])
    with draw.bound(1):
        picked = table.draw("rocky")
    assert picked.probability == pytest.approx({"a": 0.75, "b": 0.25}[picked.name])
    assert dict(picked.probabilities) == pytest.approx({"a": 0.75, "b": 0.25})


def test_a_draw_is_the_same_from_the_same_seed_in_any_order_of_other_draws():
    table = Table("t", 1, [row(c, weight=1.0) for c in "abcdefgh"])

    def run(seed, noise):
        with draw.bound(seed):
            if noise:
                draw.random()
            return [table.draw("rocky").name for _ in range(30)]

    assert run("s", False) == run("s", False)
    assert run("s", False) != run("t", False)
    assert run("s", True) != run("s", False)                            # the stream is the only input


def test_the_table_is_walked_in_its_fixed_row_order_and_a_reorder_changes_a_seed_s_draw():
    rows = [row("a"), row("b"), row("c")]
    one, other = Table("t", 1, rows), Table("t", 1, list(reversed(rows)))
    with draw.bound(5):
        first = [one.draw("rocky").name for _ in range(40)]
    with draw.bound(5):
        second = [other.draw("rocky").name for _ in range(40)]
    assert first != second
    assert [r.name for r, _p in one.probabilities(Subject("rocky"))] == ["a", "b", "c"]


# --- The fractions inside a row -----------------------------------------------------------------------------------


def test_a_sourced_kernel_gives_the_fractions_directly_scaled_to_one():
    r = row("a", fractions={"olivine": 3.0, "pyroxene": 1.0})
    assert r.draw_fractions() == {"olivine": 0.75, "pyroxene": 0.25}


def test_a_range_only_row_draws_its_fractions_from_a_dirichlet_around_the_weights():
    r = row("a", fractions={"olivine": 0.6, "pyroxene": 0.3, "metal": 0.1}, concentration=40.0)
    draws = []
    with draw.bound("dirichlet-row"):
        for _ in range(4000):
            f = r.draw_fractions()
            assert list(f) == ["olivine", "pyroxene", "metal"]
            assert math.fsum(f.values()) == pytest.approx(1.0, abs=1e-12)
            draws.append(f)
    assert statistics.fmean(f["olivine"] for f in draws) == pytest.approx(0.6, abs=0.01)
    assert statistics.fmean(f["metal"] for f in draws) == pytest.approx(0.1, abs=0.01)
    assert len({f["olivine"] for f in draws}) > 3000
    with draw.bound("dirichlet-row"):
        assert r.draw_fractions() == draws[0]


def test_a_larger_concentration_spreads_the_fractions_less():
    def spread(c):
        r = row("a", fractions={"x": 0.5, "y": 0.5}, concentration=c)
        with draw.bound("spread"):
            return statistics.pstdev(r.draw_fractions()["x"] for _ in range(2000))

    assert spread(200.0) < spread(20.0) < spread(2.0)


def test_the_concentration_follows_from_the_spread_of_a_fraction():
    c = pt.concentration_from_spread(0.3, 0.05)
    assert 0.3 * 0.7 / (c + 1) == pytest.approx(0.05 ** 2)
    with pytest.raises(ValueError):
        pt.concentration_from_spread(0.3, 0.5)
    with pytest.raises(ValueError):
        pt.concentration_from_spread(1.0, 0.1)
    with pytest.raises(ValueError):
        pt.concentration_from_spread(0.3, 0.0)


def test_bad_fractions_are_refused():
    for bad in ({"a": -1.0}, {"a": 0.0, "b": 0.0}, {"a": math.nan}):
        with pytest.raises(TableError, match="fractions"):
            row("a", fractions=bad)
    with pytest.raises(TableError, match="concentration needs"):
        row("a", concentration=5.0)
    for c in (0.0, -2.0, math.inf):
        with pytest.raises(TableError, match="positive"):
            row("a", fractions={"x": 1.0}, concentration=c)
    with pytest.raises(TableError, match="no fractions"):
        row("a").draw_fractions()


# --- The shared components table ----------------------------------------------------------------------------------


def test_a_table_with_fractions_must_name_the_shared_components_table():
    with pytest.raises(TableError, match="component list owned by one kind"):
        Table("belts", 1, [row("a", fractions={"olivine": 1.0})])


def test_the_registry_checks_every_component_against_the_one_shared_table():
    registry = Registry()
    registry.register(Table(pt.SHARED_COMPONENTS, 1, [row("olivine", (ANY_CLASS,)), row("water-ice", (ANY_CLASS,))]))
    registry.register(Table("belts", 1, [row("a", fractions={"olivine": 1.0, "water-ice": 1.0})],
                            components=pt.SHARED_COMPONENTS))
    registry.validate()
    registry.register(Table("moons", 1, [row("b", fractions={"olivine": 1.0, "unobtainium": 1.0})],
                            components=pt.SHARED_COMPONENTS))
    with pytest.raises(TableError, match="unobtainium"):
        registry.validate()


def test_a_table_may_not_take_its_components_from_a_list_of_its_own():
    registry = Registry()
    registry.register(Table(pt.SHARED_COMPONENTS, 1, [row("olivine", (ANY_CLASS,))]))
    registry.register(Table("comets", 1, [row("a", fractions={"olivine": 1.0})], components="comet-materials"))
    with pytest.raises(TableError, match="one shared table"):
        registry.validate()
    solo = Registry()
    solo.register(Table("belts", 1, [row("a", fractions={"olivine": 1.0})], components=pt.SHARED_COMPONENTS))
    with pytest.raises(TableError, match="not registered"):
        solo.validate()


# --- The registry, the report and the version -----------------------------------------------------------------------


def test_tables_register_themselves_and_a_name_is_taken_once():
    registry = Registry()
    t = registry.register(Table("a", 1, [row("x")]))
    registry.register(Table("b", 1, [row("y")]))
    assert registry.tables()[0] is t and [x.name for x in registry.tables()] == ["a", "b"]
    assert "a" in registry and len(registry) == 2 and registry.get("b").name == "b"
    with pytest.raises(TableError, match="already registered"):
        registry.register(Table("a", 1, [row("z")]))
    with pytest.raises(KeyError):
        registry.get("c")
    registry.unregister("a")
    assert "a" not in registry


def test_the_global_registry_starts_with_no_test_tables():
    assert "t" not in pt.REGISTRY and isinstance(pt.REGISTRY, Registry)


def test_the_report_counts_the_states_and_lists_every_open_weight_and_stand_in():
    registry = Registry()
    registry.register(Table("atmospheres", 3, [
        row("a"), open_row("x", ("rocky", "icy")),
        row("s", state=STAND_IN, replacement="GEN.473"),
        row("p", state=PARTLY_VERIFIED, limit="below 10 bar")]))
    text = registry.report()
    assert "atmospheres (version 3): 1 verified, 1 partly verified, 1 open, 1 stand-in" in text
    assert "Open weights: 1\n  atmospheres: x [icy, rocky]" in text
    assert "Stand-ins: 1\n  atmospheres: s [rocky] replaced by GEN.473" in text
    assert "limit: below 10 bar" in text
    assert registry.report() == text


def test_the_digest_changes_with_any_row_change_and_the_versions_list_it():
    base = [row("a", weight=1.0), row("b", weight=2.0, fractions={"x": 1.0})]
    t1 = Table("t", 1, base, components="m")
    assert t1.digest() == Table("t", 1, list(base), components="m").digest()
    changed = [
        [row("a", weight=1.5), base[1]],
        [base[0], row("b", weight=2.0, fractions={"x": 1.0, "y": 1.0})],
        [base[0], row("b", ("icy",), 2.0, fractions={"x": 1.0})],
        [base[0]],
        [base[1], base[0]],
    ]
    digests = {Table("t", 1, rows, components="m").digest() for rows in changed}
    assert t1.digest() not in digests and len(digests) == len(changed)
    assert Table("t", 2, base, components="m").digest() != t1.digest()
    registry = Registry()
    registry.register(t1)
    assert registry.versions() == {"t": (1, t1.digest())}


def test_state_counts_name_every_state():
    table = Table("t", 1, [row("a"), open_row("x")])
    assert table.state_counts() == {"verified": 1, "partly verified": 0, "open": 1, "stand-in": 0}
