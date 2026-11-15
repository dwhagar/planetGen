# planetgen/queue/needs.py

"""
What each command needs and makes (OPS.238, first part: the declarations). A command that `makes` an artifact leaves it stored; a command that `needs` one cannot run before something has made it. An artifact is a
name for a stored result of the galaxy, not a table: "skeleton" is the density skeleton the plan stores, "sectors" the sectors the generators save, "bright-stars-cleared" the emptied bright-star scatter that a mass pass
draws into. The names are the vocabulary of this module (`ARTIFACTS`); a test holds every command to it, every need to a maker, and the fixed order of a galaxy run (`generation.stages.galaxy_stages`) to the
declarations, so the order a planner chooses (the second part of OPS.238) is checked against the order the runs have today. Nothing reads this to start a run yet: `commands.REQUIRES` still holds the stages that wait
for another, and goes when the planner reads these.

A need is met in a run when no command of the run makes it (it is stored from an earlier run: a scatter of one pass needs the skeleton, which `planetgen plan` made) or a command of the run that makes it comes earlier.
"""

ARTIFACTS = {
    "skeleton": "the galaxy's density skeleton (`plan-skeleton`)",
    "start-sector": "the sector the galaxy starts in, which the neighborhood is found around",
    "sectors": "the sectors the generators saved",
    "neighbor-links": "the links between neighboring sectors",
    "sector-paths": "the sector paths that are saved last",
    "bright-stars-cleared": "the bright-star scatter emptied before a pass draws it again",
    "scatter-limits": "the limits the star scatter runs at, stored before the backfill",
    "bright-stars-backfill": "the massive stars drawn from the neighborhood",
    "bright-stars-mass": "the mass pass of the bright stars",
    "bright-stars-luminosity": "the luminosity pass of the bright stars",
    "bright-stars-band": "the dimmer band of bright stars",
    "phenomena-cleared": "the phenomenon scatter emptied before it is drawn again",
    "phenomenon-indexes-off": "the phenomenon table's indexes taken off before the shells go in",
    "phenomena-shells": "the phenomena of the shells, before the special ones and the indexes",
    "phenomena": "the phenomenon scatter finished: the special phenomena, the indexes and the record",
    "nebula-field-cleared": "the nebula field table emptied before it is drawn again",
    "nebula-field": "the nebula field",
    "population": "the species, civilizations and polities",
    "territories": "which polity owns which system",
    "run-blocks": "the Galaxy Map's block totals up to what the run stored",
    "block-lights": "the Galaxy Map's block light totals",
}
"""dict: Every artifact a command can need or make, with what it is."""

MAKES = {
    "plan-skeleton": ("skeleton",),
    "generate-start-sector": ("start-sector", "sectors"),
    "generate-sectors": ("start-sector", "sectors"),
    "generate-neighborhood": ("sectors",),
    "generate-unplaced-sector": ("sectors",),
    "link-neighbors": ("neighbor-links",),
    "settle-sector-paths": ("sector-paths",),
    "clear-bright-stars": ("bright-stars-cleared",),
    "store-scatter-limits": ("scatter-limits",),
    "scatter-backfill": ("bright-stars-backfill",),
    "scatter-mass": ("bright-stars-mass",),
    "scatter-luminosity": ("bright-stars-luminosity",),
    "scatter-band": ("bright-stars-band",),
    "empty-phenomenon-scatter": ("phenomena-cleared",),
    "drop-phenomenon-indexes": ("phenomenon-indexes-off",),
    "scatter-phenomena": ("phenomena-shells",),
    "finish-phenomena": ("phenomena",),
    "empty-nebula-field": ("nebula-field-cleared",),
    "draw-nebula-field": ("nebula-field",),
    "run-population": ("population",),
    "run-territories": ("territories",),
    "refresh-run-blocks": ("run-blocks",),
    "refresh-block-lights": ("block-lights",),
}
"""dict: The artifacts a command makes, by command kind; a command that is not here makes none (the clears of a scatter, the warm-ups, the checks, the fills)."""

NEEDS = {
    "generate-start-sector": ("skeleton",),
    "generate-sectors": ("skeleton",),
    "generate-neighborhood": ("skeleton", "start-sector"),
    "generate-unplaced-sector": ("skeleton",),
    "link-neighbors": ("sectors",),
    "settle-sector-paths": ("sectors",),
    "store-scatter-limits": ("skeleton", "bright-stars-cleared"),
    "scatter-backfill": ("scatter-limits", "sectors"),
    "scatter-mass": ("skeleton", "bright-stars-cleared"),
    "scatter-luminosity": ("bright-stars-mass",),
    "scatter-band": ("bright-stars-luminosity",),
    "drop-phenomenon-indexes": ("phenomena-cleared",),
    "scatter-phenomena": ("skeleton", "phenomenon-indexes-off"),
    "finish-phenomena": ("phenomena-shells",),
    "draw-nebula-field": ("skeleton", "nebula-field-cleared"),
    "run-population": ("sectors",),
    "run-territories": ("population",),
    "refresh-run-blocks": ("sectors",),
    "refresh-block-lights": ("bright-stars-mass", "bright-stars-luminosity", "phenomena", "run-blocks"),
}
"""dict: The artifacts a command needs, by command kind; a command that is not here needs none that a run could lack."""


def makes_of(kind):
    """The artifacts command `kind` makes, `()` for none."""
    return MAKES.get(kind, ())


def needs_of(kind):
    """The artifacts command `kind` needs, `()` for none."""
    return NEEDS.get(kind, ())


def makers_of(artifact):
    """The kinds of command that make `artifact`, in the order of `MAKES`."""
    return [kind for kind, made in MAKES.items() if artifact in made]


def edges():
    """`{(needer, maker)}`: every pair of commands where the first needs what the second makes."""
    return {(kind, maker) for kind, needed in NEEDS.items() for artifact in needed for maker in makers_of(artifact)}
