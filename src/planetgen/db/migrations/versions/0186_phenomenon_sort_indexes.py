"""The Phenomena table's sorts by placed, descriptor and radius read in index order (PERF.211).

Schema v186. No column changes and no stored value changes: each of the eight phenomenon tables gets an index ending in `name` for the placed sort
(`(center_x_pc IS NOT NULL)`) and one for the descriptor sort, and the three tables with a light-year radius one for the radius sort. The
placed and descriptor keys of some tables are expressions (functional indexes); the compact remnants lead with `star_id`, as the list leaves out the
ones that belong to a star system. The indexes on columns build online (`LOCK=NONE`); the server will not build one on an expression that way, so those
let readers carry on and make writers wait until they are built (`LOCK=SHARED`). A stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0186"
down_revision = "0185"
branch_labels = None
depends_on = None

PLACED = "(center_x_pc IS NOT NULL)"
COLLATE = " COLLATE utf8mb4_unicode_ci"
DESCRIPTORS = {
    "nebulae": "nebula_type",
    "asteroid_fields": "density",
    "black_holes": "(CASE WHEN has_accretion_disk THEN 'accreting' ELSE 'quiescent' END)" + COLLATE,
    "neutron_stars": "pulsar_type",
    "supernova_remnants": "morphology",
    "rogue_planets": ("(CASE WHEN mass_bin = 'brown-dwarf' THEN 'brown dwarf' WHEN planet_type = 'g' THEN 'gas giant' "
                      "ELSE 'terrestrial' END)" + COLLATE),
    "interstellar_comets": "(CASE WHEN is_active THEN 'active' ELSE 'dormant' END)" + COLLATE,
    "quasars": "(CASE WHEN is_radio_loud THEN 'radio-loud' ELSE 'radio-quiet' END)" + COLLATE,
}
LEADING = {"black_holes": ("star_id",), "neutron_stars": ("star_id",)}
RADIUS_TABLES = ("nebulae", "asteroid_fields", "supernova_remnants")


def upgrade():
    connection = op.get_bind()
    for table, descriptor in DESCRIPTORS.items():
        lead = list(LEADING.get(table, ()))
        expressions = {f"idx_{table}_placed_name": lead + [PLACED, "name"]}
        if descriptor.startswith("("):
            expressions[f"idx_{table}_descriptor_name"] = lead + [descriptor, "name"]
        helpers.add_expression_indexes(connection, table, expressions)
        if not descriptor.startswith("("):
            helpers.add_index_online(connection, table, f"idx_{table}_descriptor_name", tuple(lead + [descriptor, "name"]))
        if table in RADIUS_TABLES:
            helpers.add_index_online(connection, table, f"idx_{table}_radius_name", ("radius_ly", "name"))
