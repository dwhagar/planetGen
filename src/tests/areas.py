# tests/areas.py

"""
Which part of planetGen a test file covers, and which files are heavy (TEST.139).

`conftest.py` reads this to add two kinds of marker, so a run can pick a slice
without moving any file (`pytest -m area_web`, `scripts/test.py area web db`):

- `area_<name>`: the first rule in `AREA_RULES` whose pattern matches the file's
  name without `test_` and `.py`. A file no rule matches is `area_other`; there
  is no error, so a new test file never fails a lane's run for lack of a line
  here, but `test_test_areas.py` lists the files that landed in `other` so they
  can be given a rule.
- `heavy`: the files in `HEAVY_FILES`, which took over 100 s each in the
  2026-10-14 baseline (docs/design/test-suite-rebuild.md). They stay out of the
  `quick` slice.

The areas are coarse on purpose. What a change really touches is decided by the
recorded dependency map (`depmap.py`, `pytest --changed`), not by this table.
"""

import re

AREA_RULES = (
    ("wiki", r"wikiclient"),
    ("names", r"names?($|_)|naming|lazy_names|dedupe_names|gated_phoneme|body_names|parallel_names"),
    ("db", r"db_|database|sql_|schema|migrate|update_database|spatial_position_db|sector_paths_db|count_cache|"
           r"bughunt_db|settings|alembic|id_blocks|tilecache|caches|slow_read|query_db"),
    ("api", r"api|ratelimit|login|password|two_factor|admin_auth|corridor|route_edge|activity"),
    ("web", r"web|fuzz_web|page|systempage|systemmap|starmap|galaxymap|galaxydrill|galaxyprisms|galaxystages|"
            r"galaxy_urls|galaxy_viewport|navmap|nav_|tile|datatable|generatebuttons|icons|mdrender|js_unit|"
            r"class_pages|console|pagination|number_format|speed_period|distance_format|fmt_|local_time|"
            r"connectionnote|system_render|system_scene|warmup|eta_|progress_|navigation|render_edges|text_rendering"),
    ("physics", r"astropy|atmosphere|climate|comet|disk|evolution|galactic|habitability|heliopause|heliosphere|"
                r"hydrosphere|kepler|light|luminosity|orbit|planet|quasar|radiation|remnant|spin|state_vectors|"
                r"surface|math|body_positions|gen_kepler|moons|facilit|wide_binary|binary_pair|rogue|spatial_position|"
                r"physical|physics_edges|planet_inputs"),
    ("generation", r"gen_|generation|galaxy_|sector|phenomen|nebula|bright|star|scatter|object|parallel|"
                   r"work_queue|worker|redis_queue|prevalence|population|mass_backfill|fingerprint|reproducible|"
                   r"neighborhood|directives|single_system|skeleton|poisson|redo|span|stages|steps|step_registry|"
                   r"near|search|save_progress|failed_run|validation|examples|full_matrix|space_sector|"
                   r"system_velocity|system_|systems|planets|benchmark|sector_|serialization|end_to_end"),
    ("ops", r"admin_script|edge_|deploy|install|bump|version|todo|update_|modulegraph|debug_log|job_tree|"
            r"private_dirs|admin|state|activity_log|cli|utils_and_config"),
)

HEAVY_FILES = frozenset({
    "test_neighborhood_edge.py", "test_galaxy_gen.py", "test_parallel_generation.py", "test_prevalence.py",
    "test_api.py", "test_phenomenon_scatter.py", "test_work_queue_failures.py", "test_galaxy_urls_combined.py",
    "test_db_old_schemas.py", "test_web_generate.py", "test_api_thin_routes.py", "test_bright_star_scatter.py",
    "test_work_queue.py", "test_admin_edits.py", "test_gen_resume.py", "test_api_auth_sweep.py",
    "test_search_edge_cases.py", "test_db_persistence.py", "test_db_schema_version_detect.py",
})

_COMPILED = tuple((area, re.compile(pattern)) for area, pattern in AREA_RULES)


def area_of(file_name):
    """The area of a test file by its name (`test_web_pages.py` -> `web`), or `other`."""
    stem = file_name[:-3] if file_name.endswith(".py") else file_name
    if stem.startswith("test_"):
        stem = stem[len("test_"):]
    for prefix in ("fuzz_", "bughunt_"):
        if stem.startswith(prefix):
            stem = stem[len(prefix):]
    for area, pattern in _COMPILED:
        if pattern.match(stem) or pattern.match(stem.split("_", 1)[0]):
            return area
    return "other"


def is_heavy(file_name):
    return file_name in HEAVY_FILES
