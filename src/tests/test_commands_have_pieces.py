"""OPS.196, first part: no command goes silent by being written as one loop of its own. Every registered command gives the engine `pieces` and `do` (`queue/engine.py`), which count them into the
unit's beat and check for a cancel or a pause between two, or is listed in `NOT_YET` with the item of the family that converts it. The list only gets shorter: a command that has been converted and
is still listed fails here, and so does a name that is not a registered command. The second part (every command run on the small test galaxy with its longest stretch of no movement held to a
budget, and the `whole` flag deleted) comes with the last family (reports/job-engine-items-stage-1-2-2026-11-13.md, item P6)."""

from planetgen.queue import commands

NOT_YET = {
    # OPS.198: scatters and fields
    **dict.fromkeys(("scatter-mass", "scatter-phenomena", "scatter-luminosity", "scatter-band", "draw-nebula-field", "fill-map-luminosity"), "OPS.198"),
    # OPS.200: clears and drops, a step to each statement through ctx.sql
    **dict.fromkeys(("clean-nebula-field", "clear-scatter", "empty-phenomenon-scatter", "drop-phenomenon-indexes", "clear-bright-stars", "empty-nebula-field"), "OPS.200"),
    # OPS.201: stamped walkers by id range
    **dict.fromkeys(("rescore-habitability", "refill-effect-ledger"), "OPS.201"),
    # OPS.202: galaxy and sectors
    **dict.fromkeys(("scatter-backfill", "refresh-run-blocks", "generate-start-sector", "generate-neighborhood", "generate-sectors", "link-neighbors", "settle-sector-paths",
                     "generate-unplaced-sector", "generate-system"), "OPS.202"),
    # OPS.203: population
    **dict.fromkeys(("run-population", "run-territories"), "OPS.203"),
    # OPS.204: the checks
    **dict.fromkeys(("check-math", "check-datasets", "check-revision", "check-tables", "check-tables-extended", "check-orphans", "check-ids", "check-values", "check-counts",
                     "check-version-keys", "check-storage", "check-systems", "check-object-index"), "OPS.204"),
    # OPS.205: warm-up
    **dict.fromkeys(("warm-cluster-field", "warm-object-index", "warm-galaxy-map", "warm-galaxy-nodes", "warm-galaxy-deep-nodes", "warm-lists", "warm-sectors", "warm-systems",
                     "warm-phenomena", "warm-search", "end-map-rebuild-maintenance"), "OPS.205"),
    # OPS.206: reset
    "reset-galaxy": "OPS.206",
}
"""dict: The commands that still have a `unit` function of their own, each with the item that converts it."""


def test_every_command_gives_the_engine_its_pieces_or_waits_for_the_item_that_converts_it():
    unlisted = sorted(kind for kind, spec in commands.REGISTRY.items() if spec.pieces is None and kind not in NOT_YET)
    assert not unlisted, f"commands with a loop of their own and no item to convert them (give them `pieces` and `do`, or list them in NOT_YET with the item): {unlisted}"


def test_a_command_that_was_converted_leaves_the_list():
    converted = sorted(kind for kind in NOT_YET if kind in commands.REGISTRY and commands.REGISTRY[kind].pieces is not None)
    assert not converted, f"converted to pieces and still listed as not yet (delete them from NOT_YET): {converted}"


def test_the_list_names_only_registered_commands():
    gone = sorted(kind for kind in NOT_YET if kind not in commands.REGISTRY)
    assert not gone, f"listed but not registered: {gone}"


def test_the_commands_converted_so_far_are_the_skeleton_and_the_lights():
    assert {kind for kind, spec in commands.REGISTRY.items() if spec.pieces is not None} >= {"plan-skeleton", "finish-phenomena", "refresh-block-lights"}
