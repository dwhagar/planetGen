"""
GEN.262: companions in the scatter passes (`planetgen.generation.scatter_companions`,
docs/design/systems/multi-star-systems.md 6 and 9). A star row always draws its member tree; a black hole, neutron star or
planetary nebula row has one at the table's odds, screened by a uniform that is a pure function of the row's seed;
the trees are stored under their rows with the barycentre bookkeeping, deleted with the rows, and the phenomenon
scatter's class rate is divided by `1 + p_placed` so the density of the class is unchanged.
"""

import math
import random
from types import SimpleNamespace
from collections import Counter

import pytest

from planetgen import tuning
from planetgen.db import scatter_systems, store
from planetgen.galaxy.scatter_limits import ScatterLimits
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import multistar, phenomenon_scatter, scatter_companions
from planetgen.physics import constants
from planetgen.physics.units import ly_to_pc
from tests.toy_shape import toy_shape

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
EXTENTS = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3)]
STAR = scatter_companions._STAR_COLUMN
PHEN = scatter_companions._PHENOMENON_COLUMN
EPOCH = 1.7e9


def _seeds(count, start=1):
    rng = random.Random(start)
    return [rng.getrandbits(63) for _ in range(count)]


def _phenomenon_row(kind, subtype, seed, mass=15.0, shell=3, band=0, wedge=0):
    row = [None] * len(store.PHENOMENON_SCATTER_COLUMNS)
    for name, value in (("shell_index", shell), ("band_index", band), ("wedge_index", wedge), ("kind", kind),
                        ("subtype", subtype), ("position_x_mpc", 1000), ("position_y_mpc", 2000),
                        ("position_z_mpc", 3000), ("seed", seed), ("mass_solar", mass)):
        row[PHEN[name]] = value
    return tuple(row)


def _star_row(seed, initial_mass=20.0, age_gy=0.004, mass=None, yerkes="V"):
    row = [None] * len(store.BRIGHT_STAR_COLUMNS)
    mass = initial_mass if mass is None else mass
    for name, value in (("shell_index", 3), ("band_index", 0), ("wedge_index", 0), ("position_x_mpc", 1000),
                        ("position_y_mpc", 2000), ("position_z_mpc", 3000), ("population", "young"),
                        ("star_type", "O5V Blue Main Sequence Star"), ("yerkes_class", yerkes),
                        ("mass_kg", mass * constants.SOLAR_MASS_TO_KG), ("radius_km", 1.0e7),
                        ("temperature_k", 4.0e4), ("luminosity_w", 1.0e32), ("age_gy", age_gy),
                        ("lifespan_gy", 0.01), ("initial_mass_sol", initial_mass), ("phase_end_age_gy", 0.011),
                        ("seed", seed)):
        row[STAR[name]] = value
    return tuple(row)


def _hit_seeds(kind, subtype, wanted, state=None, start=1):
    """`wanted` seeds whose screen names a companion (of `state`, when given) for the class."""
    found, rng = [], random.Random(start)
    while len(found) < wanted:
        seed = rng.getrandbits(63)
        hit = scatter_companions.phenomenon_companion(kind, subtype, seed)
        if hit is not None and (state is None or hit[2] == state):
            found.append(seed)
    return found


# --- The screen ---------------------------------------------------------------------------------------------------


def test_the_screen_is_a_pure_uniform_of_the_seed_and_the_question():
    assert scatter_companions.screen_unit(12345) == scatter_companions.screen_unit(12345)
    assert scatter_companions.screen_unit(12345) == pytest.approx(0.1330796686614273, abs=1e-15)
    assert scatter_companions.screen_unit(12345, "recycled") == pytest.approx(0.6020397784669149, abs=1e-15)
    units = [scatter_companions.screen_unit(seed) for seed in _seeds(20000)]
    assert all(0.0 <= u < 1.0 for u in units)
    assert abs(sum(units) / len(units) - 0.5) < 4 * math.sqrt(1.0 / 12.0 / len(units))
    deciles = Counter(int(u * 10) for u in units)
    assert all(abs(count - 2000) < 4 * math.sqrt(2000 * 0.9) for count in deciles.values())
    other = [scatter_companions.screen_unit(seed, "recycled") for seed in _seeds(20000)]
    correlation = sum((a - 0.5) * (b - 0.5) for a, b in zip(units, other)) / len(units) * 12.0
    assert abs(correlation) < 4 / math.sqrt(len(units))


def _within_four_sigma(count, total, share):
    mean = total * share
    return abs(count - mean) <= 4.0 * math.sqrt(total * share * (1.0 - share)) + 1.0


@pytest.mark.parametrize("kind, subtype, placed_class", [
    ("black-hole", "stellar", "black-hole"), ("black-hole", "intermediate", "intermediate-mass-black-hole"),
    ("planetary-nebula", None, "planetary-nebula-star")])
def test_the_share_of_rows_with_each_companion_state_is_the_tables_within_four_sigma(kind, subtype, placed_class):
    total = 300_000
    counts = Counter()
    for seed in _seeds(total, start=11):
        hit = scatter_companions.phenomenon_companion(kind, subtype, seed)
        counts[None if hit is None else hit[2]] += 1
    table = multistar.MEMBER_STATE_BY_PLACED[placed_class]
    for state, (odds, _low, _high) in table.items():
        assert _within_four_sigma(counts[state], total, odds), (state, counts[state], total * odds)
    assert _within_four_sigma(counts[None], total, 1.0 - sum(odds for odds, _l, _h in table.values()))
    assert set(counts) <= set(table) | {None}


def test_a_neutron_star_row_is_a_millisecond_pulsar_at_the_pulsar_chances_and_takes_that_tables_odds():
    total = 400_000
    seeds = _seeds(total, start=21)
    recycled = [s for s in seeds if scatter_companions.screen_unit(s, "recycled")
                < phenomenon_scatter.NEUTRON_STAR_RECYCLED_SHARE]
    share = phenomenon_scatter.NEUTRON_STAR_RECYCLED_SHARE
    assert _within_four_sigma(len(recycled), total, share)
    states = Counter()
    for seed in recycled:
        hit = scatter_companions.phenomenon_companion("neutron-star", None, seed)
        assert hit is None or hit[1] is True
        states[None if hit is None else hit[2]] += 1
    # a recycled pulsar has a companion about 78% of the time (a white dwarf mostly)
    assert states["white-dwarf"] > 0.45 * len(recycled)
    young = [s for s in seeds[:100_000] if s not in set(recycled)]
    star = sum(1 for s in young if (scatter_companions.phenomenon_companion("neutron-star", None, s) or (0, 0, 0))[2] == "star")
    assert _within_four_sigma(star, len(young), 0.02)


@pytest.mark.parametrize("kind, subtype", [
    ("supernova-remnant", None), ("hypervelocity-star", None), ("quasar", None), ("black-hole", "supermassive")])
def test_a_class_without_a_tree_is_never_a_hit(kind, subtype):
    assert all(scatter_companions.phenomenon_companion(kind, subtype, seed) is None for seed in _seeds(5000))
    assert scatter_companions.phenomenon_tree(_phenomenon_row(kind, subtype, 99)) is None


# --- The trees ----------------------------------------------------------------------------------------------------


def test_a_hit_draws_the_companion_the_screen_named():
    expected_kind = {"black-hole": lambda m: m.kind == "black-hole", "neutron-star": lambda m: m.kind == "neutron-star",
                     "star": lambda m: m.kind not in multistar.COMPACT_KINDS}
    tried = Counter()
    kept = Counter()
    for state in ("star", "black-hole", "neutron-star"):
        for seed in _hit_seeds("black-hole", "stellar", 120, state=state, start=31):
            tried[state] += 1
            tree = scatter_companions.phenomenon_tree(_phenomenon_row("black-hole", "stellar", seed))
            if tree is not None:
                kept[state] += 1
                assert tree.members[0].kind == "black-hole" and tree.member_count >= 2
                assert any(expected_kind[state](member) for member in tree.members[1:]), state
    # a pair whose orbit finds no stable, long-lived place is left as a single object
    assert all(kept[state] > 0.6 * tried[state] for state in tried), (tried, kept)


def test_a_rows_tree_is_a_function_of_its_seed_and_nothing_else():
    seeds = _hit_seeds("black-hole", "stellar", 12, start=41)
    rows = [_phenomenon_row("black-hole", "stellar", seed) for seed in seeds]
    first = [scatter_companions.phenomenon_tree(row) for row in rows]
    second = [scatter_companions.phenomenon_tree(row) for row in reversed(rows)][::-1]
    assert first == second
    moved = [scatter_companions.phenomenon_tree(_phenomenon_row("black-hole", "stellar", seed, shell=40, band=2, wedge=9))
             for seed in seeds]
    assert first == moved                   # the address is not an input
    stars = [_star_row(seed) for seed in _seeds(6, start=43)]
    assert [scatter_companions.star_tree(row) for row in stars] == [scatter_companions.star_tree(row) for row in stars]
    assert scatter_companions.star_tree(_star_row(1)) != scatter_companions.star_tree(_star_row(2))


def test_the_tree_is_drawn_from_the_members_stream_of_the_seed():
    seed = _hit_seeds("black-hole", "stellar", 1, state="star", start=47)[0]
    row = _phenomenon_row("black-hole", "stellar", seed, mass=16.0)
    rng = scatter_companions.draw.Stream(f"{seed}:members")
    age = rng.uniform(*tuning.COMPACT_REMNANT_AGE_RANGE_GY) * 1.0e9
    by_hand = multistar.draw_member_tree(multistar.Placed("black-hole", 16.0), age, rng, companion="star")
    assert scatter_companions.phenomenon_tree(row) == (by_hand if by_hand.nodes else None)


def test_a_planetary_nebulas_central_star_is_a_white_dwarf_with_a_star_or_white_dwarf_companion():
    seen = Counter()
    for seed in _hit_seeds("planetary-nebula", None, 150, start=51):
        tree = scatter_companions.phenomenon_tree(_phenomenon_row("planetary-nebula", None, seed, mass=None))
        if tree is not None:
            assert tree.members[0].kind == "white-dwarf"
            assert tuning.WD_MASS_RANGE_SOL[0] <= tree.members[0].mass_solar <= tuning.WD_MASS_RANGE_SOL[1]
            seen[tree.members[1].kind] += 1
    assert sum(seen.values()) > 100 and set(seen) <= {"main-sequence", "giant", "white-dwarf", "brown-dwarf",
                                                      "pre-main-sequence"}


def test_a_star_row_draws_the_star_count_of_its_spectral_class():
    counts = [scatter_companions.star_tree(_star_row(seed)) for seed in _seeds(600, start=61)]
    stars = [1 if tree is None else tree.star_count for tree in counts]
    weights = multistar.star_count_weights(20.0)
    total = sum(weights)
    expected = sum((k + 1) * w / total for k, w in enumerate(weights))
    variance = sum((k + 1 - expected) ** 2 * w / total for k, w in enumerate(weights))
    mean = sum(stars) / len(stars)
    # the compact companion adds a member that is not a star, so the stars alone follow the table
    assert abs(mean - expected) < 3.0 * math.sqrt(variance / len(stars)) + 0.05, (mean, expected)
    assert all(tree is None or tree.members[0].kind == "main-sequence" for tree in counts)


def test_a_white_dwarf_row_is_the_placed_object_of_the_white_dwarf_table():
    trees = [scatter_companions.star_tree(_star_row(seed, initial_mass=3.0, age_gy=5.0, mass=0.65, yerkes="VII"))
             for seed in _seeds(400, start=67)]
    with_company = [tree for tree in trees if tree is not None]
    assert all(tree.members[0].kind == "white-dwarf" for tree in with_company)
    assert _within_four_sigma(len(with_company), len(trees), 0.30) or len(with_company) > 0.1 * len(trees)


# --- Rates --------------------------------------------------------------------------------------------------------


def test_a_rows_rate_is_the_class_rate_over_one_plus_the_same_class_companions_per_row():
    assert phenomenon_scatter.row_density_divisor("black-hole", "stellar") == pytest.approx(1.002)
    young, recycled = 1.8e-4, 2.0e-3
    share = phenomenon_scatter.NEUTRON_STAR_RECYCLED_SHARE
    assert phenomenon_scatter.row_density_divisor("neutron-star", None) == pytest.approx(
        1.0 + (1.0 - share) * young + share * recycled)
    for kind, subtype in (("black-hole", "intermediate"), ("planetary-nebula", None), ("supernova-remnant", None)):
        assert phenomenon_scatter.row_density_divisor(kind, subtype) == 1.0
    rates = {(kind, subtype): rate for kind, subtype, rate in phenomenon_scatter._scattered_rates(1.0)}
    for kind, subtype in phenomenon_scatter.SCATTER_CLASSES:
        row_rate = rates[(kind, subtype)]
        object_rate = (phenomenon_scatter.class_rate_per_star(kind, subtype)
                       * phenomenon_scatter.share_above(kind, subtype, 1.0))
        assert row_rate * phenomenon_scatter.row_density_divisor(kind, subtype) == pytest.approx(object_rate)


def test_rows_and_their_same_class_companions_keep_the_classs_object_density():
    # Over many seeded black-hole rows, the black-hole companions the screen names add p_placed objects per row.
    total = 600_000
    companions = sum(1 for seed in _seeds(total, start=71)
                     if (scatter_companions.phenomenon_companion("black-hole", "stellar", seed) or (0, 0, ""))[2]
                     == "black-hole")
    p_placed = multistar.same_class_companions("black-hole")
    assert _within_four_sigma(companions, total, p_placed)
    rows = 1.0 / (1.0 + p_placed)
    assert rows + rows * p_placed == pytest.approx(1.0)


def test_the_scatter_draws_the_same_rows_it_drew_before_companions():
    rows = list(phenomenon_scatter.scatter_shell(SHAPE, 8, 3, EDGE_PC, 600.0, 5, min_mass_solar=1.0))
    assert rows and len(rows[0]) == len(store.PHENOMENON_SCATTER_COLUMNS)
    assert rows == list(phenomenon_scatter.scatter_shell(SHAPE, 8, 3, EDGE_PC, 600.0, 5, min_mass_solar=1.0))


# --- Storage ------------------------------------------------------------------------------------------------------


def _members(conn, table, row_id):
    return conn.execute("SELECT * FROM scatter_members WHERE source = ? AND source_id = ? ORDER BY member_no",
                        (scatter_systems.SOURCES[table], row_id)).fetchall()


def _nodes(conn, table, row_id):
    return conn.execute("SELECT * FROM scatter_nodes WHERE source = ? AND source_id = ? ORDER BY node_no",
                        (scatter_systems.SOURCES[table], row_id)).fetchall()


def _count(conn, table):
    return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]


def test_hits_are_stored_under_the_ids_their_rows_got_and_the_rest_stay_single(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            before = [_phenomenon_row("black-hole", "stellar", s, mass=14.0 + i % 5, band=i % 3)
                      for i, s in enumerate(_seeds(30, start=81))]
            store.insert_phenomenon_scatter(conn, before)            # rows from another writer, already in the table
        hits = _hit_seeds("black-hole", "stellar", 25, start=83)
        mixed = []
        for i, seed in enumerate(hits):
            mixed.append(_phenomenon_row("black-hole", "stellar", seed, mass=15.0 + i % 4, band=i % 3))
            mixed.extend(_phenomenon_row("black-hole", "stellar", s, mass=16.0)
                         for s in _seeds(2, start=1000 + i) if scatter_companions.phenomenon_companion(
                             "black-hole", "stellar", s) is None)
        with conn:
            written = scatter_companions.insert_phenomenon_scatter(conn, mixed, EPOCH, batch_size=17)
        assert written.rows == len(mixed) and 0 < written.systems <= len(hits)
        rows = {r["seed"]: r for r in conn.execute("SELECT * FROM phenomenon_scatter").fetchall()}
        with_members = 0
        for seed, row in rows.items():
            members, nodes = _members(conn, "phenomenon_scatter", row["id"]), _nodes(conn, "phenomenon_scatter", row["id"])
            if row["member_count"] is None:
                assert not members and not nodes and row["system_mass_solar"] is None
                continue
            with_members += 1
            assert seed in hits
            tree = scatter_companions.phenomenon_tree(next(r for r in mixed if r[PHEN["seed"]] == seed))
            assert row["member_count"] == len(members) == tree.member_count
            assert len(nodes) == len(tree.nodes) == len(members) - 1
            assert row["system_mass_solar"] == pytest.approx(sum(m["mass_solar"] for m in members))
            assert members[0]["kind"] == "black-hole" and members[0]["mass_solar"] == pytest.approx(row["mass_solar"])
            assert {n["epoch_unix"] for n in nodes} == {EPOCH}
        assert with_members == written.systems
    finally:
        conn.close()


def test_a_star_pass_stores_every_trees_members_under_the_right_star(mysql_config):
    rows = list(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, 2000.0, 2500.0, 7))
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            written = scatter_companions.insert_bright_stars(conn, rows, EPOCH, batch_size=300)
        assert written.rows == len(rows) == _count(conn, "bright_stars")
        stored = {r["seed"]: r for r in conn.execute("SELECT * FROM bright_stars").fetchall()}
        expected = {row[STAR["seed"]]: scatter_companions.star_tree(row) for row in rows}
        assert written.systems == sum(tree is not None for tree in expected.values()) > 100
        for seed, tree in expected.items():
            row = stored[seed]
            if tree is None:
                assert row["member_count"] is None and not _members(conn, "bright_stars", row["id"])
                continue
            members = _members(conn, "bright_stars", row["id"])
            assert row["member_count"] == tree.member_count == len(members)
            assert row["system_mass_solar"] == pytest.approx(tree.system_mass_solar)
            assert [m["seed"] for m in members] == [m.seed for m in tree.members]
            assert members[0]["mass_solar"] == pytest.approx(row["mass_kg"] / constants.SOLAR_MASS_TO_KG)
            assert not {m["kind"] for m in members} & {"planet", "moon", "asteroid-belt"}
    finally:
        conn.close()


def test_two_rows_with_one_seed_each_get_the_same_tree(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            written = scatter_companions.insert_bright_stars(conn, [_star_row(77), _star_row(77)], EPOCH)
        ids = [r["id"] for r in conn.execute("SELECT id FROM bright_stars ORDER BY id").fetchall()]
        assert written.rows == 2 and written.systems == (2 if scatter_companions.star_tree(_star_row(77)) else 0)
        if written.systems:
            first, second = (_members(conn, "bright_stars", i) for i in ids)
            assert [m["mass_solar"] for m in first] == [m["mass_solar"] for m in second]
    finally:
        conn.close()


def test_a_batch_without_hits_writes_the_rows_and_nothing_else(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        silent = [s for s in _seeds(200, start=91) if scatter_companions.phenomenon_companion("black-hole", "stellar", s) is None]
        with conn:
            written = scatter_companions.insert_phenomenon_scatter(
                conn, [_phenomenon_row("black-hole", "stellar", s) for s in silent])
        assert written == (len(silent), 0)
        assert _count(conn, "scatter_members") == _count(conn, "scatter_nodes") == 0
    finally:
        conn.close()


def test_the_default_epoch_is_a_fixed_origin_until_an_orbit_update_has_set_the_clock(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        assert store.get_orbit_epoch_unix(conn) is None
        assert scatter_companions.default_epoch_unix(conn) == tuning.SCATTER_PHASE_ORIGIN_UNIX
        with conn:
            store.orbit_clock(conn, 1.0)
        epoch = store.get_orbit_epoch_unix(conn)
        if epoch is not None:
            assert scatter_companions.default_epoch_unix(conn) == float(epoch)
    finally:
        conn.close()


# --- Members go with their rows -----------------------------------------------------------------------------------


def _seed_with_members(conn, mass_pass=3, light_pass=3):
    rows = ([_star_row(s, initial_mass=20.0) for s in _seeds(mass_pass * 6, start=101)]
            + [_star_row(s, initial_mass=6.0, age_gy=0.05) for s in _seeds(light_pass * 6, start=103)])
    with conn:
        written = scatter_companions.insert_bright_stars(conn, rows, EPOCH)
    assert written.systems > 4
    return rows


def test_clearing_the_star_scatter_clears_its_members_but_not_the_phenomenon_scatters(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        _seed_with_members(conn)
        with conn:
            scatter_companions.insert_phenomenon_scatter(
                conn, [_phenomenon_row("black-hole", "stellar", s) for s in _hit_seeds("black-hole", "stellar", 6)], EPOCH)
        both = conn.execute("SELECT source, COUNT(*) AS n FROM scatter_members GROUP BY source").fetchall()
        assert {r["source"] for r in both} == {1, 2}
        store.clear_bright_stars(conn)
        left = conn.execute("SELECT source, COUNT(*) AS n FROM scatter_members GROUP BY source").fetchall()
        assert {r["source"] for r in left} == {2}
        assert not conn.execute("SELECT 1 FROM scatter_nodes WHERE source = 1").fetchall()
        store.clear_phenomenon_scatter(conn)
        assert _count(conn, "scatter_members") == _count(conn, "scatter_nodes") == 0
    finally:
        conn.close()


def _limits():
    """The scatter limits of a galaxy scattered at the defaults (no stored skeleton here, so the core and halo are
    the defaults': the star rows seeded are far from the 14 to 15.4 solar mass window they differ in)."""
    skeleton = SimpleNamespace(shape=SHAPE, edge_pc=EDGE_PC, core_min_luminosity_sol=None, core_mass_limit_sol=None,
                               halo_min_luminosity_sol=tuning.HALO_MIN_LUMINOSITY_SOL)
    return ScatterLimits.of(skeleton, 9000.0, 14.0)


@pytest.mark.parametrize("part, gone, kept", [("mass", 20.0, 6.0), ("luminosity", 6.0, 20.0)])
def test_clearing_one_star_pass_clears_only_that_passs_members(mysql_config, part, gone, kept):
    conn = store.get_connection(mysql_config)
    try:
        _seed_with_members(conn)
        store.clear_bright_star_pass(conn, part, _limits())
        initial = {r["initial_mass_sol"] for r in conn.execute("SELECT initial_mass_sol FROM bright_stars").fetchall()}
        assert initial == {kept}
        owners = {r["source_id"] for r in conn.execute("SELECT DISTINCT source_id FROM scatter_members").fetchall()}
        assert owners and owners <= {r["id"] for r in conn.execute("SELECT id FROM bright_stars").fetchall()}
        assert owners == {r["source_id"] for r in conn.execute("SELECT DISTINCT source_id FROM scatter_nodes").fetchall()}
    finally:
        conn.close()


def test_deleting_the_unbuilt_stars_of_a_cell_deletes_their_members_too(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        rows = [_star_row(s) for s in _seeds(40, start=111)]
        elsewhere = [r[:STAR["shell_index"]] + (9,) + r[STAR["shell_index"] + 1:] for r in
                     [_star_row(s) for s in _seeds(40, start=113)]]
        with conn:
            scatter_companions.insert_bright_stars(conn, rows + elsewhere, EPOCH)
        before = _count(conn, "scatter_members")
        store.delete_unbuilt_bright_stars(conn, (3, 0, 0))
        conn.commit()
        assert 0 < _count(conn, "scatter_members") < before
        survivors = {r["id"] for r in conn.execute("SELECT id FROM bright_stars").fetchall()}
        assert {r["source_id"] for r in conn.execute("SELECT source_id FROM scatter_members").fetchall()} <= survivors
        assert {r["source_id"] for r in conn.execute("SELECT source_id FROM scatter_nodes").fetchall()} <= survivors
    finally:
        conn.close()


def test_the_nightly_update_moves_the_barycentre_and_writes_no_node(mysql_config):
    from tests.test_galactic_motion import SCATTER_RADIUS_PC, _seed_scatter_grid
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        row = list(_star_row(5))
        row[STAR["position_x_mpc"]], row[STAR["position_y_mpc"]], row[STAR["position_z_mpc"]] = (
            int(SCATTER_RADIUS_PC * 1000), 0, 0)
        row[STAR["shell_index"]] = row[STAR["band_index"]] = row[STAR["wedge_index"]] = 0
        seeds = [s for s in _seeds(40, start=121) if scatter_companions.star_tree(_star_row(s)) is not None][:3]
        rows = []
        for s in seeds:
            copy = list(row)
            copy[STAR["seed"]] = s
            rows.append(tuple(copy))
        with conn:
            written = scatter_companions.insert_bright_stars(conn, rows, EPOCH)
        assert written.systems == 3
        nodes_before = conn.execute("SELECT * FROM scatter_nodes ORDER BY source_id, node_no").fetchall()
        clock = store.orbit_clock(conn, 1.0e6)
        positions = {r["id"]: r["position_x_mpc"] for r in conn.execute("SELECT id, position_x_mpc FROM bright_stars").fetchall()}
        with conn:
            store.advance_scatter_points(conn, clock)
        moved = {r["id"]: r["position_x_mpc"] for r in conn.execute("SELECT id, position_x_mpc FROM bright_stars").fetchall()}
        assert moved != positions
        assert conn.execute("SELECT * FROM scatter_nodes ORDER BY source_id, node_no").fetchall() == nodes_before
    finally:
        conn.close()


# --- A plan run ---------------------------------------------------------------------------------------------------


def _members_by_row_seed(conn, table):
    """`{row seed: [(kind, mass, member seed), ...]}` of every stored system, for comparing two runs."""
    systems = {}
    for row in conn.execute(
            f"SELECT t.seed AS seed, m.member_no AS member_no, m.kind AS kind, m.mass_solar AS mass, m.seed AS member_seed"
            f" FROM {table} t JOIN scatter_members m ON m.source = ? AND m.source_id = t.id"
            " ORDER BY t.seed, m.member_no", (scatter_systems.SOURCES[table],)).fetchall():
        systems.setdefault(row["seed"], []).append((row["kind"], row["mass"], row["member_seed"]))
    return systems


def test_a_phenomenon_plan_stores_the_companions_and_a_rerun_stores_the_same_ones(mysql_config):
    from tests.test_phenomenon_scatter import _scatter_args, _seed_galaxy
    from planetgen.generation import run_plan
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    conn = store.get_connection(mysql_config)
    try:
        first = _members_by_row_seed(conn, "phenomenon_scatter")
        rows = conn.execute("SELECT kind, subtype, member_count, system_mass_solar FROM phenomenon_scatter").fetchall()
        with_members = [r for r in rows if r["member_count"]]
        assert len(first) == len(with_members) > 20
        assert {r["kind"] for r in with_members} <= {"black-hole", "neutron-star", "planetary-nebula"}
        assert not [r for r in rows if r["kind"] in ("hypervelocity-star", "quasar", "supernova-remnant")
                    and r["member_count"]]
        assert _count(conn, "scatter_members") == sum(len(v) for v in first.values())
        assert _count(conn, "scatter_nodes") == _count(conn, "scatter_members") - len(first)
    finally:
        conn.close()
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    conn = store.get_connection(mysql_config)
    try:
        assert _members_by_row_seed(conn, "phenomenon_scatter") == first
        # the rerun replaced the rows, and with them their members: nothing is left over from the first run
        assert _count(conn, "scatter_members") == sum(len(v) for v in first.values())
    finally:
        conn.close()


def test_a_star_plan_stores_the_companions_and_the_phenomena_plan_leaves_them_alone(mysql_config):
    from tests.test_bright_star_scatter import THRESHOLD, _scatter_args, _seed_galaxy
    from tests.test_phenomenon_scatter import _scatter_args as phenomenon_args
    from planetgen.generation import run_plan
    _seed_galaxy(mysql_config)
    summary = run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--workers", "1"))
    conn = store.get_connection(mysql_config)
    try:
        systems = _members_by_row_seed(conn, "bright_stars")
        assert summary["total"] > 20 and len(systems) > 0.4 * summary["total"]
        assert all(len(members) >= 2 for members in systems.values())
        before = _count(conn, "scatter_members")
    finally:
        conn.close()
    run_plan.scatter_phenomena(phenomenon_args(mysql_config))     # source 2 rows come and go; source 1 stays
    conn = store.get_connection(mysql_config)
    try:
        assert _members_by_row_seed(conn, "bright_stars") == systems
        assert _count(conn, "scatter_members") > before
        assert THRESHOLD > 0
    finally:
        conn.close()
