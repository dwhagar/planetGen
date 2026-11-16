# tests/pair_overlap_recipes.py

"""
The recipes of the pair proof (OPS.242, `tests/pair_overlap_support.py`): for each pair of commands listed in `queue/independence.json`, how to put a toy galaxy in the state both stages work on
and which command lines make them.
"""

from planetgen.db import store
from tests.bughunt_support import mysql_argv, run_cli
from tests.pair_overlap_support import PairRecipe
from tests.test_bright_star_scatter import SMALL_E_VALUE, THRESHOLD, _seed_galaxy


def _prepare(context):
    """The toy galaxy of the bright star scatter tests (`test_bright_star_scatter._seed_galaxy`), five times as dense, so that every scatter draws rows."""
    _seed_galaxy(context.mysql_config, e_value=SMALL_E_VALUE * 5)


def _scatters(context):
    """The command lines that run the stages of the pairs for real, so the capture sees the params each is cut with."""
    for flags in (["--mass", "--phenomenon-min-mass", "8"], ["--phenomena"], ["--nebula-field"]):
        run_cli("scatter", [*flags, "--bright-star-min-luminosity", str(THRESHOLD), *mysql_argv(context.mysql_config)])


def _through_owner(context, child, parent, source):
    """
    `(name, sql)` of the rows of `child` (`scatter_members`, `scatter_nodes`) of the scatter rows `parent` (`bright_stars`, `phenomenon_scatter`) holds, as the owner's `uid` and every column of
    the row but the owner's row id (`source_id`): the id is the order the rows landed in, which no output depends on, and the `uid` is what names the row.
    """
    conn = store.get_connection(context.mysql_config)
    try:
        columns = [row["column_name"] for row in conn.execute(
            "SELECT column_name AS column_name FROM information_schema.columns WHERE table_schema = DATABASE() AND table_name = ? AND column_name <> 'source_id' ORDER BY ordinal_position",
            (child,)).fetchall()]
    finally:
        conn.close()
    return (f"{child} of {parent}", f"SELECT p.uid AS owner_uid, {', '.join(f'c.`{name}`' for name in columns)} FROM {child} c JOIN {parent} p ON c.source = {source} AND c.source_id = p.id")


def _scatter_views(context):
    from planetgen.db import scatter_systems
    return tuple(_through_owner(context, child, parent, scatter_systems.SOURCES[parent]) for child in ("scatter_members", "scatter_nodes") for parent in ("bright_stars", "phenomenon_scatter"))


SKELETON = {"prepare": _prepare, "capture": _scatters, "views": _scatter_views,
            "ignore_columns": (("bright_stars", "id"), ("phenomenon_scatter", "id"), ("scatter_members", "source_id"), ("scatter_nodes", "source_id"))}
IDS = ("`id` of a scatter row is its row number, the order the rows landed in; `source_id` of a member or a node is that number, so the members and nodes are compared through the owner's `uid`, "
       "which the row's seed and rank give.")

PAIRS = {
    ("scatter-phenomena", "draw-nebula-field"): PairRecipe(
        "scatter-phenomena", "draw-nebula-field", **SKELETON,
        note="The phenomena scatter writes its own table, the members and nodes of what it drew and the light parts of the blocks; the nebula field is its own table and its stage's own rows. " + IDS),
    ("scatter-mass", "scatter-phenomena"): PairRecipe(
        "scatter-mass", "scatter-phenomena", **SKELETON,
        note="The mass pass writes `bright_stars`, the phenomena scatter `phenomenon_scatter`; both write members, nodes, block light parts and the ids of their own rows (`id_counters`: a counter each, "
             "the stars' and the phenomena's), and each reads the other's table only to give the blocks their light. " + IDS),
}
"""dict: `{(first, second): recipe}`, the commands in plan order, one for each pair in `queue/independence.json`."""


def _mass_then_luminosity(context):
    for flags in (["--mass", "--phenomenon-min-mass", "8"], ["--luminosity"]):
        run_cli("scatter", [*flags, "--bright-star-min-luminosity", str(THRESHOLD), *mysql_argv(context.mysql_config)])


CONTROL = PairRecipe(
    "scatter-mass", "scatter-luminosity", prepare=_prepare, capture=_mass_then_luminosity, views=_scatter_views, ignore_columns=SKELETON["ignore_columns"],
    note="Not a pair: the luminosity pass leaves out the sectors the mass pass marked, read when it is cut, so overlapped (cut before the mass pass has stored a star) it draws what plan order does not. "
         "The control of the proof: it must fail.")
"""PairRecipe: two stages that are not independent, to show that the proof fails when it should."""
