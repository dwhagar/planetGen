"""
GEN.1682: the neutron stars and black holes below the phenomenon scatter's cut that a backfill ring's sectors hold before they are filled. The draw is the fill's own (`below_cut_draws`) on the sector's own
stream, the rows land in the scatter table at points inside the sector's cell, a sector that holds some is not drawn again, and the fill builds them and draws none of its own.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy.spherical import sector_address_at, sector_position_pc
from planetgen.generation import phenomenon_scatter as scatter
from planetgen.generation import ring_remnants, run_galaxy, run_plan
from planetgen.generation import scatter_companions
from tests import first_fill
from tests.test_phenomenon_scatter import EDGE_PC, E_VALUE, SHAPE, _rows, _scatter_args, _seed_galaxy

CUT = 20.0
"""float: The cut the tests scatter at: the fixed presets run from 8 up, and 20 leaves every neutron star and a good share of the stellar black holes below it."""
RING = [(2, 0, 3), (3, 1, 2), (4, 0, 5)]
"""list: Sectors of the toy galaxy the tests take for backfill ring sectors."""


def _cell_center(address):
    return sector_position_pc(*address, EDGE_PC)


def _draw(address, seed=7, expected=20000.0, cut=CUT):
    return scatter.ring_remnant_rows(address, _cell_center(address), SHAPE, expected, cut, seed, EDGE_PC)


def test_the_ring_rows_are_the_fills_own_draw_in_the_sectors_own_cell():
    address = (2, 0, 3)
    rows = [dict(zip(scatter.PHENOMENON_SCATTER_COLUMNS, row)) for row in _draw(address)]
    draws, _rng = scatter.below_cut_draws(address, _cell_center(address), SHAPE, 20000.0, CUT, 7)
    wanted = {(kind, subtype): count for kind, subtype, bounds, count in draws if bounds is not None}
    assert sum(wanted.values()) > 20, "the toy sector should hold enough remnants for the test to mean something"
    got = {}
    for row in rows:
        got[(row["kind"], row["subtype"])] = got.get((row["kind"], row["subtype"]), 0) + 1
        assert (row["shell_index"], row["band_index"], row["wedge_index"]) == address
        assert tuple(sector_address_at((row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0), EDGE_PC)) == address
    assert set(got) <= set(wanted)
    # a point the cell's redraws could not find is left out, so a count is at most the draw's and nearly always equal to it
    assert all(got.get(key, 0) <= count for key, count in wanted.items())
    assert sum(got.values()) >= 0.95 * sum(wanted.values())
    for row in rows:
        low, high = next(bounds for kind, subtype, bounds, _count in draws if (kind, subtype) == (row["kind"], row["subtype"]))
        assert low <= row["mass_solar"] <= high
        assert row["mass_solar"] < CUT


def test_the_ring_rows_repeat_with_the_seed_and_the_sector_and_none_come_from_a_cut_under_every_mass_law():
    address = (2, 0, 3)
    assert _draw(address) == _draw(address)
    assert _draw(address, seed=8) != _draw(address)
    assert _draw((3, 1, 2)) != _draw(address)
    assert _draw(address, cut=1.0) == []


def test_a_neutron_star_row_below_the_cut_is_a_neutron_stars_mass():
    rows = [dict(zip(scatter.PHENOMENON_SCATTER_COLUMNS, row)) for row in _draw((2, 0, 3))]
    neutron = [row for row in rows if row["kind"] == "neutron-star"]
    assert neutron and all(1.1 <= row["mass_solar"] <= 2.2 for row in neutron)
    assert all(row["subtype"] != "magnetar" for row in rows)


def _backfill(mysql_config, addresses, mass=2.0):
    """Marks `addresses` as sectors a mass backfill took down (`sector_stats.bright_mass_sol`), as a real backfill's chunk does."""
    conn = store.get_connection(mysql_config)
    try:
        entries = [(address, 1.0, E_VALUE) for address in addresses]
        store.lock_sector_stats(conn, entries)
        store.set_sector_bright_masses(conn, {address: (mass, None) for address in addresses})
        conn.commit()
    finally:
        conn.close()


def _remnants_in(mysql_config, address):
    return [row for row in _rows(mysql_config, "shell_index = %d AND band_index = %d AND wedge_index = %d AND mass_solar < %r" % (*address, CUT))
            if row["kind"] in store.BELOW_CUT_REMNANT_KINDS]


def test_the_phenomena_unit_draws_the_ring_sectors_remnants_and_counts_their_classes(mysql_config):
    _seed_galaxy(mysql_config)
    _backfill(mysql_config, RING)
    run_plan.scatter_phenomena(_scatter_args(mysql_config, "--phenomenon-min-mass", str(int(CUT))))
    conn = store.get_connection(mysql_config)
    try:
        held = store.addresses_holding_below_cut_remnants(conn, CUT)
        classes = {(row["kind"], row["subtype"]): row["unbuilt"] for row in store.phenomenon_scatter_classes(conn)}
        total = conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter WHERE mass_solar < ? AND kind IN ('neutron-star', 'black-hole') AND COALESCE(subtype, '') <> 'magnetar'", (CUT,)).fetchone()["n"]
        stored = {kind: count for kind, count in ((row["kind"], row["n"]) for row in conn.execute("SELECT kind, COUNT(*) AS n FROM phenomenon_scatter WHERE mass_solar < ? GROUP BY kind", (CUT,)).fetchall())}
    finally:
        conn.close()
    assert held == set(RING), "each ring sector holds some, and a sector no backfill touched holds none"
    assert total == len(sum((_remnants_in(mysql_config, address) for address in RING), []))
    assert classes[("neutron-star", "")] >= stored["neutron-star"], "the classes count the rows, so the Phenomena table counts them as unbuilt"
    for address in RING:
        assert _remnants_in(mysql_config, address), address


def test_a_charted_sector_and_a_sector_already_holding_remnants_are_not_drawn(mysql_config):
    _seed_galaxy(mysql_config)
    _backfill(mysql_config, RING)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("INSERT INTO chart_ledger (shell_index, band_index, wedge_index) VALUES (?, ?, ?)", RING[0])
        conn.commit()
        assert store.backfilled_uncharted_addresses(conn) == sorted(RING[1:])
        skeleton = store.get_galaxy_shape(conn)
        limits = store.scatter_limits_of(conn, skeleton)
        rows = ring_remnants.rows_for(conn, skeleton, limits, 7, CUT, RING[1:], held=set())
        again = ring_remnants.rows_for(conn, skeleton, limits, 7, CUT, RING[1:], held={RING[1]})
        assert rows and {tuple(row[:3]) for row in rows} == set(RING[1:])
        assert {tuple(row[:3]) for row in again} == {RING[2]}
        assert ring_remnants.rows_for(conn, skeleton, limits, 7, None, RING[1:]) == []
    finally:
        conn.close()


def test_a_backfill_after_the_scatter_draws_the_remnants_in_its_own_transaction_and_the_fill_makes_no_second_draw(mysql_config, monkeypatch):
    first_fill.without_the_deepest_floor(monkeypatch)
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config, "--phenomenon-min-mass", str(int(CUT))))
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        limits, seed = run_plan._scatter_limits_and_seed(conn, skeleton)
        assert not _remnants_in(mysql_config, RING[0])
        run_plan._draw_sector_masses(conn, skeleton, {address: 2.0 for address in RING}, limits, seed)
        held = store.addresses_holding_below_cut_remnants(conn, CUT)
        assert held <= set(RING) and RING[0] in held
        # the same backfill again, and the same chunk: nothing is drawn twice
        before = len(_rows(mysql_config))
        run_plan._draw_sector_masses(conn, skeleton, {address: 5.0 for address in RING}, limits, seed)
        assert len(_rows(mysql_config)) == before
    finally:
        conn.close()
    # the fill builds the rows at their stored points and draws no remnants of its own
    context = run_galaxy._fill_context(run_galaxy._default_generation_args(config=mysql_config), RING[0])
    assert context.below_cut is None
    assert [row["kind"] for row in context.phenomenon_rows if row["kind"] in store.BELOW_CUT_REMNANT_KINDS]
    other = (1, 0, 1)
    assert run_galaxy._fill_context(run_galaxy._default_generation_args(config=mysql_config), other).below_cut is not None


def test_a_backfill_while_no_scatter_is_recorded_draws_none(mysql_config, monkeypatch):
    first_fill.without_the_deepest_floor(monkeypatch)
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        limits, seed = run_plan._scatter_limits_and_seed(conn, skeleton)
        run_plan._draw_sector_masses(conn, skeleton, {address: 2.0 for address in RING}, limits, seed)
        assert conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter").fetchone()["n"] == 0
    finally:
        conn.close()


def _cuts(baseline):
    """Cuts of the toy galaxy whose core and halo cuts are far from the baseline's, so a sector (all of them baseline bands here) must take the baseline's."""
    from planetgen.galaxy.scatter_limits import ScatterLimits, TierLimit
    from planetgen.generation.phenomenon_cuts import Cuts
    tiers = (3.0, baseline, 30.0)
    return Cuts(tiers, ScatterLimits(SHAPE, EDGE_PC, tuple(TierLimit(1.0, cut) for cut in tiers)))


def test_gen377_a_ring_sector_draws_what_lies_under_its_own_bands_cut(mysql_config):
    """The cut is the sector's band's tier cut: the sectors of the toy galaxy are baseline bands, so `Cuts` with other cuts in the core and the halo draws what the baseline's cut alone does."""
    _seed_galaxy(mysql_config)
    _backfill(mysql_config, RING)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        limits = store.scatter_limits_of(conn, skeleton)
        by_number = ring_remnants.rows_for(conn, skeleton, limits, 7, 12.0, RING, held=set())
        by_cuts = ring_remnants.rows_for(conn, skeleton, limits, 7, _cuts(12.0), RING, held=set())
        assert by_cuts == by_number and by_cuts
        assert all(row[12] < 12.0 for row in by_cuts if row[3] in store.BELOW_CUT_REMNANT_KINDS)
        assert by_cuts != ring_remnants.rows_for(conn, skeleton, limits, 7, 20.0, RING, held=set())
    finally:
        conn.close()


def test_gen377_a_cell_holds_below_cut_remnants_when_a_row_is_under_its_own_cut(mysql_config):
    _seed_galaxy(mysql_config)
    _backfill(mysql_config, RING)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        limits = store.scatter_limits_of(conn, skeleton)
        rows = ring_remnants.rows_for(conn, skeleton, limits, 7, 20.0, RING, held=set())
        from planetgen.generation import scatter_companions
        scatter_companions.insert_phenomenon_scatter(conn, rows)
        conn.commit()
        assert store.addresses_holding_below_cut_remnants(conn, _cuts(20.0)) == store.addresses_holding_below_cut_remnants(conn, 20.0) == set(RING)
        assert store.addresses_holding_below_cut_remnants(conn, _cuts(1.0)) == set(), "every row is over a cut of 1 solar mass, whatever the other tiers' cuts are"
        masses = sorted(row[12] for row in rows if row[3] in store.BELOW_CUT_REMNANT_KINDS)
        middle = masses[len(masses) // 2]
        some = store.addresses_holding_below_cut_remnants(conn, _cuts(middle))
        assert some <= set(RING) and some == store.addresses_holding_below_cut_remnants(conn, middle)
    finally:
        conn.close()
