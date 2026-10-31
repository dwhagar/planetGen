# tests/test_db_column_round_trip.py

"""
TEST.11: every column of every content table is written by the code that
saves its rows and read back by the code that loads them. A small galaxy
is built through the real `planetgen` commands (`rich_galaxy_support`),
then every column in information_schema must hold a value in some row
(so a new column the INSERT forgot, left NULL, fails), and must be read
by the loaders and the API's detail endpoints (so a column nothing ever
loads fails). The exceptions below each say why.
"""

import pytest

from planetgen import tuning
from planetgen.physics import rings
from planetgen.db import store
from tests.column_read_support import columns_read, tracking_reads
from tests.rich_galaxy_support import rich_galaxy  # noqa: F401  (fixture)
from tests.publicids import pid

pytestmark = pytest.mark.db

BOOKKEEPING_TABLES = {
    "schema_migrations", "alembic_version", "generation_run_arguments", "phenomenon_scatter", "phenomenon_scatter_classes", "sector_system_counts", "system_interests", "nebula_field", "galaxy_components", "facility_types", "facility_type_fields", "facility_field_values", "galaxy_block_stats", "galaxy_block_lights", "region_stats", "object_timing_stats", "object_timing_samples", "id_blocks", "id_counters", "system_name_registry", "sector_name_registry", "population_state",
    "orbit_simulation_state", "galaxy_column", "galaxy_layer", "sector_stats", "bright_stars",
    "nearest_systems", "sector_paths", "sector_path_knots", "scatter_members", "scatter_nodes",
    "orbit_nodes", "dark_members", "object_index", "galaxy_containers", "vessels", "constellations", "request_comments", "request_notifications", "region_density_overrides", "object_locks", "change_plans", "change_plan_ops", "edit_journal", "sector_tree_nodes", "species_limit_values", "galaxy_shape_pins", "nebula_block_stats", "nebula_block_classes", "nebula_block_state",
    "orbital_events", "orbital_event_bodies",
}
"""Tables no object loader reads row by row: migration and id bookkeeping,
the name registries, generation progress, and the galaxy skeleton and
bright-star scatter (read in bulk by the Galaxy Map's tile queries,
covered by test_galaxymap3d.py/test_bright_star_scatter.py), the galaxy's components (read by `store.get_galaxy_components`, covered by test_migration_galaxy_components.py) and the saved
sector paths (read by `db/sector_paths.py`, covered by test_sector_paths_db.py) and the orbit tree's nodes (DB.53: the
migration builds them, the generator and the loaders follow in the DB.53 chain; covered by test_migration_orbit_nodes.py) and the dark members of a tree (GEN.1013: few systems have one; covered by test_gen_stored_trees.py) and the nebula totals per block of the Galaxy Map's ladder (MAP.172: the plan keeps them with the field; covered by test_nebula_blocks_db.py)."""

CHANCE = "chance"
"""Reason for a column whose value in this galaxy depends on the draw."""

CONTAINMENT = ("inside_nebula_id", "inside_remnant_id")
SECTOR_PLACEMENT = ("sector_id", "center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc", "quadrant")

NULL_IN_THIS_GALAXY = {
    # DB.81: the stellar mass a plan pinned, which a galaxy whose seed drew everything leaves NULL (test_shape_pins.py
    # writes and reads it back).
    ("galaxy_shape", "pinned_stellar_mass_solar"): "NULL when the seed drew the mass (DB.81)",
    # DB.44: nothing writes the checksum until DB.9's storage protection does; a sector is saved with `parity_dirty` 1.
    ("sectors", "storage_checksum"): "written by DB.9's storage protection (DB.44 adds the column only)",
    # DB.60: nothing writes the density columns until a sector fill reads the region overrides (MAP.119).
    **{("sectors", column): "written by the sector fill under a density override (DB.60 adds the columns only)"
       for column in ("density_multiplier_used", "override_id")},
    # DB.66: the bar and the thick disk were never stored, so a spiral leaves them NULL ("as the structure draw gives them").
    **{("galaxy_shape", column): "NULL means as the structure draw gives it from the seed (DB.66)" for column in (
        "bar_present", "bar_axis_ratio_b", "bar_axis_ratio_c", "bar_angle_rad", "thick_disk_amplitude")},
    # DB.61: nothing writes the cluster link until the cluster membership of a star is generated (GEN.209).
    ("stars", "cluster_id"): "written when a star's cluster is known (DB.61 adds the column only)",
    # GEN.267: only `planetgen system --stars N` stores a forced star count; sector and galaxy runs refuse it.
    ("system_configs", "stars"): "written by a single system run with --stars (sector and galaxy runs refuse it)",
    # Only something sitting inside a nebula or supernova remnant.
    **{(table, column): CHANCE for table in store.CONTAINABLE_TABLES for column in CONTAINMENT},
    # A quasar is placed only at a galaxy's own center, so whether the run
    # wrote a placed one is chance; the CLI's is standalone (TEST.84).
    **{("quasars", column): CHANCE for column in SECTOR_PLACEMENT},
    ("quasars", "jet_length_ly"): CHANCE,
    # The CLI can't place a supernova remnant, and whether it has a core is chance.
    **{("supernova_remnants", column): "standalone remnant" for column in SECTOR_PLACEMENT},
    **{("supernova_remnants", column): CHANCE for column in
       ("compact_remnant_kind", "compact_remnant_black_hole_id", "compact_remnant_neutron_star_id")},
    # `plan --no-bright-stars`: the scatter is slow and tested on its own.
    ("galaxy_shape", "bright_star_min_luminosity_sol"): "no scatter", ("galaxy_shape", "bright_star_seed"): "no scatter",
    ("galaxy_shape", "bright_star_mass_limit_sol"): "no scatter",
    ("galaxy_shape", "phenomenon_scatter_seed"): "no scatter",
    ("galaxy_shape", "phenomenon_min_mass_solar"): "no scatter",
    ("galaxy_shape", "core_min_luminosity_sol"): "follows the galaxy-wide limit until set (DB.54)",
    ("galaxy_shape", "core_mass_limit_sol"): "follows the galaxy-wide limit until set (DB.54)",
    # High-pressure ice lies only under the deepest oceans.
    # GEN.1057: a rogue planet has no moon that fell and no comet that passed, so it has no ring (the columns are stored for a ring set by hand).
    **{("rogue_planets", column): "a rogue planet has no ring (GEN.1057)" for column in rings.BODY_FIELDS},
    **{(table, "hp_ice_km"): CHANCE for table in ("planets", "moons", "rogue_planets")},
    # Only a heap of grains lighter than SMALL_BODY_MASS_KG has a porosity, and a planet that small is a rare draw (a moon's is set by its origin).
    ("planets", "porosity"): CHANCE,
    # Whether any moon of this small galaxy holds an ocean is chance.
    **{("moons", column): CHANCE for column in
       ("ocean_class", "ocean_depth_km", "ocean_ph", "phosphorus", "water_activity", "tidal_flux_w_m2")},
    # A moon with a gas envelope (a class I, J or T moon) is rare: the envelope split is stored for every body that has one.
    **{("moons", column): CHANCE for column in ("envelope_helium_fraction", "envelope_z_fraction")},
    # Only a world with no liquid at its surface can hold an aquifer, and whether any of this small galaxy's does is chance.
    **{(table, column): CHANCE for table in ("planets", "moons") for column in ("aquifer_depth_km", "aquifer_thickness_km")},
    # GEN.1035: only a world with a liquid at its surface stores one, and the depth of a lava surface is not known.
    **{(table, column): CHANCE for table in ("planets", "moons") for column in ("surface_liquid", "surface_liquid_fraction", "surface_liquid_depth_km")},
    ("planets", "moon_origin"): "a planet has no origin as a moon (DB.70)",
    **{("moons", column): "a moon has no circumplanetary disk of its own (DB.70)" for column in (
        "moon_disk_lifetime_myr", "moon_disk_alpha_over_f", "moon_disk_f_inf", "moon_disk_ice_line_km")},
    # Facilities on stars, moons and asteroid fields, and the field's galaxy position.
    **{("facilities", column): "facility hosts used" for column in
       ("star_id", "moon_id", "asteroid_field_id", "sector_id", "center_x_pc", "center_y_pc", "center_z_pc",
        "galactic_radius_pc")},
    # A life-bearing moon (rare) and its paragraphs and spectrum.
    **{("moons", column): CHANCE for column in
       ("atm_density", "atm_molar_density", "scale_height_km", "life_class", "life_class_variant", "life_class_best", "life_class_best_score",
        "life_animal_grade", "life_complexity",
        "life_stage", "evolutionary_pace", "life_milestone_age_gyr", "life_system_age_gyr", *tuning.LIFE_SCHEDULE_FIELDS)},
    **{("moon_reflection_spectrum", column): CHANCE for column in
       ("id", "moon_id", "spectrum_type", "position", "value")},
    # A giant moon (rare) and its interior (DB.89).
    **{("moons", column): CHANCE for column in (
        "interior_kind", "core_mass_earth", "dilute_core_mass_earth", "metallization_pressure_gpa", "helium_depletion",
        "has_helium_shell")},
    # A moon has a field only when it is over `magnetism.GANYMEDE_MIN_MASS_EARTH` and wins `GANYMEDE_CHANCE`, and only then a magnetopause (GEN.843).
    ("moons", "magnetopause_rp"): CHANCE,
    # A moon has no ring (GEN.732).
    **{("moons", column): CHANCE for column in rings.BODY_FIELDS},
    **{("moon_interior_layers", column): CHANCE for column in (
        "id", "moon_id", "layer_index", "kind", "state", "top_pressure_gpa", "bottom_pressure_gpa", "mass_earth",
        "top_radius_fraction", "bottom_radius_fraction", "hydrogen", "helium", "iron", "rock", "water_liquid", "ice_x",
        "nh3_ice", "ch4_ice")},
    # The same for a planet: whether any planet of this small galaxy bears life is chance.
    **{("planets", column): CHANCE for column in
       ("evolutionary_pace", "flavor_key", "life_class", "life_class_variant", "life_class_best", "life_class_best_score",
        "life_animal_grade", "life_complexity", "life_milestone_age_gyr", "life_stage",
        "life_system_age_gyr")},
    # The ages from complex cells on exist only for a world whose air holds or builds free oxygen (GEN.750), which is
    # chance for the few living worlds of this small galaxy (TEST.402).
    **{("planets", column): CHANCE for column in
       ("life_age_complex_cells_gyr", "life_age_multicellularity_gyr", "life_age_technological_civilization_gyr")},
    # A young star's disk (GEN.787): whether any star of this small galaxy is young enough to have one is chance.
    **{("stars", column): CHANCE for column in ("disk_kind", "disk_mass_solar", "disk_inner_au", "disk_outer_au")},
    **{("planet_reflection_spectrum", column): CHANCE for column in
       ("id", "planet_id", "spectrum_type", "position", "value")},
    # The CLI places an interstellar comet near the galaxy's sectors only sometimes.
    **{("interstellar_comets", column): CHANCE for column in SECTOR_PLACEMENT},
    # Only a lit class (C to G) has host stars (GEN.150), and the seeded galaxy's few nebulae may all be dark.
    **{("nebulae", column): CHANCE for column in ("age_myr", "q_total", "hosts")},
    # The cell key is derived from the center, so it is NULL wherever the center is.
    ("facilities", "cell_key"): "facility hosts used", ("quasars", "cell_key"): CHANCE,
    ("supernova_remnants", "cell_key"): "standalone remnant",
    ("interstellar_comets", "cell_key"): CHANCE,
    ("star_systems", "system_flavor_text"): CHANCE,
    ("star_systems", "runaway_class"): CHANCE, ("star_systems", "runaway_speed_kms"): CHANCE,
    ("system_configs", "age"): "only with --age",
    # GEN.178, GEN.981: a moon and a rocky planet have an induced field; a giant has none.
    ("planets", "induced_moment_a_m2"): CHANCE,
    # v89 (DB.28): the facility type, affiliation and ratings the POP.10 pass and admin pages set.
    **{("facilities", column): "not set yet" for column in (
        "type_id", "affiliation_polity_id", "affiliation_name", "rating_crime", "rating_housing",
        "rating_resources", "rating_maintenance", "rating_health", "ratings_source", "ratings_note",
        "ratings_set_at")},
    # DB.55 and DB.53 added these; the code items that fill them follow (the chain of GEN.294).
    **{("asteroid_belts", column): "a ring about its star: no motion or clock of its own (GEN.337)"
       for column in ("epoch_unix", "next_update_due")},
    **{("quasars", column): CHANCE for column in ("epoch_unix", "next_update_due")},
    **{("asteroid_fields", column): "formed by an event, which has no table yet (GEN.110, DB.56)"
       for column in ("formed_at", "origin_event_id")},
    **{("sectors", column): "worked out by GEN.109" for column in ("max_mass_id", "max_mass_kg")},
    **{("stars", column): "a member of a built multi-star system, filled by the DB.53 chain" for column in (
        "epoch_unix", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms")},
    ("planets", "host_node_id"): "a circumbinary planet, filled by the DB.53 chain",
    # DB.73: a free planet or comet is always `free`, so its parent_id stays NULL; only a moon of a free planet names one.
    **{(table, "parent_id"): "a free object has no parent row" for table in ("interstellar_comets", "rogue_planets")},
    ("moons", "rogue_planet_id"): "only a moon of a free planet (test_parent_link)",
    # GEN.764: the seeded galaxy runs no scatter, and its population pass finds no civilization.
    **{("scatter_epochs", column): "no scatter in this galaxy (GEN.764)" for column in
       ("kind", "generator_epoch", "run_id", "recorded_at")},
    ("species", "generator_epoch"): "no species in this galaxy (GEN.764)",
    ("planets", "parent_shadow_share"): "only a moon is in its planet's shadow (GEN.767)",
}
"""`(table, column)` allowed to be NULL in every row of the seeded galaxy,
with why: each is set only by something this small galaxy doesn't
contain (the cases are exercised in their own tests)."""

NEVER_READ = {
    # Row order in a child table: the loaders read them ORDER BY position.
    **{(table, "position"): "sort key" for table in (
        "asteroid_belt_composition", "asteroid_field_composition", "comet_composition",
        "interstellar_comet_composition", "planet_reflection_spectrum", "moon_reflection_spectrum")},
    # Child tables' own ids and parent keys: matched in the WHERE clause.
    **{(table, "id"): "child row id" for table in (
        "asteroid_belt_composition", "asteroid_field_composition", "comet_composition",
        "interstellar_comet_composition", "planet_reflection_spectrum", "moon_reflection_spectrum",
        "planet_interior_layers", "moon_interior_layers", "planet_gases", "moon_gases", "planet_life_class_scores", "moon_life_class_scores", "system_config_slots", "system_configs")},
    **{(table, parent): "parent key" for table, parent in (
        ("comet_composition", "comet_id"), ("interstellar_comet_composition", "comet_id"),
        ("asteroid_field_composition", "field_id"),
        ("planet_reflection_spectrum", "planet_id"),
        ("moon_reflection_spectrum", "moon_id"), ("system_config_slots", "config_id"),
        ("planet_interior_layers", "planet_id"), ("moon_interior_layers", "moon_id"),
        ("system_owners", "star_system_id"), ("nebula_shape_balls", "nebula_id"))},
    ("star_systems", "schema_version"): "the version that wrote the row, for diagnosis",
    **{("sectors", column): "read by DB.9's storage protection (DB.44 adds the columns only)"
       for column in ("storage_checksum", "parity_dirty")},
    ("stars", "cluster_id"): "read when a star's cluster is shown (DB.61 adds the column only)",
    **{("sectors", column): "read by the rebuild of a density override (DB.60 adds the columns only)"
       for column in ("density_multiplier_used", "override_id")},
    ("sectors", "leaf_digest"): "read by the fingerprint tree (GEN.305)",
    ("sectors", "run_id"): "read in a query by the version check, not by a loader (DB.16)",
    # GEN.764: the stage epochs are read by `galaxy/epoch_record.py` and `bring_to_epoch`, not by a page's loader.
    **{("scatter_epochs", column): "read by galaxy/epoch_record (GEN.764)" for column in
       ("kind", "generator_epoch", "run_id", "recorded_at")},
    ("species", "generator_epoch"): "read by galaxy/epoch_record (GEN.764)",
    ("galaxy_shape", "id"): "singleton key", ("galaxy_shape", "bright_star_seed"): "only to repeat a scatter",
    ("galaxy_shape", "phenomenon_scatter_seed"): "only to repeat a scatter",
    ("galaxy_shape", "phenomenon_min_mass_solar"): "read by a sector fill, not a loader",
    ("facilities", "galactic_radius_pc"): "an index column; pages place a facility by its center",
    ("facilities", "cell_key"): "an index column MySQL fills; the region query reads it (NAV.43)",
    **{("star_systems", column): "the point a system sits at in the galaxy frame; the region query and the orbit update read it"
       for column in ("cell_key", "galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc")},
    # GEN.106: the orbit update's own clock, read only by planetgen.cli.orbits.
    **{(table, column): "the orbit update's clock" for table, column in (
        ("star_systems", "epoch_unix"), ("star_systems", "next_update_due"),
        ("star_systems", "binary_epoch_unix"), ("star_systems", "binary_next_update_due"),
        ("facilities", "epoch_unix"), ("facilities", "next_update_due"))},
    ("moon_reflection_spectrum", "spectrum_type"): CHANCE, ("moon_reflection_spectrum", "value"): CHANCE,
    **{(table, "layer_index"): "sort key" for table in ("planet_interior_layers", "moon_interior_layers")},
    # No moon of this small galaxy is a giant (DB.89): its layers are read only when one is.
    **{("moon_interior_layers", column): CHANCE for column in (
        "kind", "state", "top_pressure_gpa", "bottom_pressure_gpa", "mass_earth", "top_radius_fraction", "bottom_radius_fraction",
        "hydrogen", "helium", "iron", "rock", "water_liquid", "ice_x", "nh3_ice", "ch4_ice")},
    ("planet_reflection_spectrum", "spectrum_type"): CHANCE, ("planet_reflection_spectrum", "value"): CHANCE,
    **{(table, stamp): "row timestamps (the Galaxy Map's change feed reads modified_at in bulk)"
       for table in ("star_systems", "facilities") for stamp in ("created_at", "modified_at")},
    ("polities", "created_at"): "row timestamp", ("species", "created_at"): "row timestamp",
    # DB.6: what made the galaxy and its runs, for a rebuild (OPS.12) and the update history (OPS.13). The
    # galaxy's own record is read back by the version check (ADM.50), so it is not listed here.
    # PERF.11: galaxy-wide stats, read by `store.galaxy_density_ratio`, not by an object loader.
    ("galaxy_shape", "density_ratio_avg"): "density stats (PERF.11)",
    ("galaxy_shape", "density_ratio_samples"): "density stats (PERF.11)",
    **{("generation_runs", column): "run history (DB.6)" for column in
       ("id", "command", "run_seed", "galaxy_seed", "started_at", "finished_at", "outcome")},
    # DB.55 and DB.53 columns nothing loads yet (the gravity sum, GEN.109 and the DB.53 chain read them).
    ("facilities", "mass_kg"): "read by the gravity sum (GEN.294 chain)",
    **{("sectors", column): "read by GEN.109" for column in ("max_mass_id", "max_mass_kg")},
    **{("stars", column): "read by the DB.53 chain" for column in ("epoch_unix",)},
    **{("stars", column): CHANCE for column in ("parent_node_id", "slot")},    # read back for a system of 3 stars or more
    # DB.73: the link is stored for the readers of 8.4 (the unified object); a system loader rebuilds the tree from parent_node_id.
    **{("stars", column): "stored for the parent-link readers (DB.73)" for column in ("parent_type", "parent_id")},
    # GEN.329: read back only for a system saved as a runaway, which this galaxy may not hold.
    ("stars", "flags"): CHANCE,
}
"""`(table, column)` the loaders and detail endpoints don't read, with why."""


def _client(config):
    from planetgen.web.app import create_app
    from planetgen.api.config import Config

    class TestConfig(Config):
        MYSQL_CONFIG = config
        WRITE_MYSQL_CONFIG = config
        CONTROL_MYSQL_CONFIG = config
        RATELIMIT_ENABLED = False

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def _content_columns(conn):
    rows = conn.execute(
        "SELECT c.TABLE_NAME AS t, c.COLUMN_NAME AS c FROM information_schema.COLUMNS c"
        " JOIN information_schema.TABLES t ON t.TABLE_SCHEMA = c.TABLE_SCHEMA AND t.TABLE_NAME = c.TABLE_NAME"
        " WHERE c.TABLE_SCHEMA = DATABASE() AND t.TABLE_TYPE = 'BASE TABLE'"
        " ORDER BY c.TABLE_NAME, c.ORDINAL_POSITION"
    ).fetchall()
    return [(row["t"], row["c"]) for row in rows if row["t"] not in BOOKKEEPING_TABLES]


def test_every_column_is_written(rich_galaxy):
    conn = store.get_connection(rich_galaxy)
    try:
        always_null = [
            (table, column) for table, column in _content_columns(conn)
            if conn.execute(f"SELECT 1 FROM `{table}` WHERE `{column}` IS NOT NULL LIMIT 1").fetchone() is None
        ]
    finally:
        conn.close()
    unexpected = sorted(set(always_null) - set(NULL_IN_THIS_GALAXY))
    assert not unexpected, f"NULL in every row (never written?): {unexpected}"
    stale = sorted(key for key, why in NULL_IN_THIS_GALAXY.items() if why != CHANCE and key not in always_null)
    assert not stale, f"written after all; drop from NULL_IN_THIS_GALAXY: {stale}"


def test_every_column_is_read_back(rich_galaxy, monkeypatch):
    client = _client(rich_galaxy)
    conn = store.get_connection(rich_galaxy)
    try:
        def ids(table):
            return [row["id"] for row in conn.execute(f"SELECT id FROM {table} ORDER BY id").fetchall()]

        from planetgen.db import query

        with tracking_reads(monkeypatch) as reads:
            pages = [f"/api/sectors/{pid('sector', i, conn)}" for i in ids("sectors")]
            pages += [f"/api/systems/{pid('system', i, conn)}{part}" for i in ids("star_systems")
                      for part in ("", "/owner", "/facilities")]
            pages += [f"/api/phenomena/{label}/{pid(label, i, conn)}" for table, label, *_ in query._PHENOMENON_TABLES
                      for i in ids(table)]
            pages += [f"/api/nebulae/{pid('nebula', i, conn)}/shape" for i in ids("nebulae")]
            pages += [f"/api/facilities/{pid('facility', i, conn)}" for i in ids("facilities")]
            pages += [f"/api/species/{i}" for i in ids("species")]
            pages += [f"/api/polities/{i}" for i in ids("polities")]
            pages += ["/api/galaxy/shape", "/api/population", "/api/territories"]
            for page in pages:
                assert client.get(page).status_code == 200, page
            for sector_id in ids("sectors"):
                store.load_sector(conn, sector_id)
            for system_id in ids("star_systems"):
                store.load_star_system(conn, system_id)
            columns = _content_columns(conn)
            read = columns_read(reads, {table for table, _column in columns})
    finally:
        conn.close()
    unread = sorted((table, column) for table, column in columns if column not in read[table])
    unexpected = sorted(set(unread) - set(NEVER_READ))
    assert not unexpected, f"never read back: {unexpected}"
    stale = sorted(key for key, why in NEVER_READ.items() if why != CHANCE and key not in unread)
    assert not stale, f"read after all; drop from NEVER_READ: {stale}"
