# planetgen/queue/needs.py

"""
What each command needs and makes (OPS.238, first part: the declarations). A command that `makes` an artifact leaves it stored; a command that `needs` one cannot run before something has made it. An artifact is a
name for a stored result of the galaxy, not a table: "skeleton" is the density skeleton the plan stores, "sectors" the sectors the generators save, "bright-stars-cleared" the emptied bright-star scatter that a mass pass
draws into. The names are the vocabulary of this module (`ARTIFACTS`); a test holds every command to it, every need to a maker, and the fixed order of a galaxy run (`generation.stages.galaxy_stages`) to the
declarations, so the order a planner chooses (the second part of OPS.238) is checked against the order the runs have today. The order of a run is planned from these (OPS.238, second part): `requires_of` is the kinds a command waits for, `required_stage_keys` the earlier steps of a run it waits for, and `validate_order`
refuses a run whose steps are listed so that a need is only met after the step that needs it.

The stored plan (`planrows.record_plan`) says which needs of a run are unmet (`unmet`: no earlier step makes them and the rows of the galaxy's earlier runs do not say they exist, `planrows.made_artifacts`, `exists_after`). The checks and the warm families declare nothing: a check reads the database as it is and a warm family builds its own input.

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
    "phenomena-special": "the special phenomena (the nucleus and the hypervelocity stars) and the backfill rings' below-cut neutron stars and black holes written with their classes",
    "phenomena-indexed": "the phenomenon table's indexes built and the hypervelocity stars' epoch stamped",
    "phenomena": "the phenomenon scatter finished: its seed, mass cut and epoch recorded",
    "nebula-field-cleared": "the nebula field table emptied before it is drawn again",
    "nebula-field": "the nebula field",
    "species": "the species the planets' scan found",
    "civilizations": "the civilizations dated and the polities founded or dissolved",
    "technology": "the map's civilization rows carrying their technology level",
    "population": "the species, civilizations and polities, counted: the population pass over",
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
    "write-special-phenomena": ("phenomena-special",),
    "build-phenomenon-indexes": ("phenomena-indexed",),
    "record-phenomenon-scatter": ("phenomena",),
    "empty-nebula-field": ("nebula-field-cleared",),
    "draw-nebula-field": ("nebula-field",),
    "scan-species": ("species",),
    "date-civilizations": ("civilizations",),
    "refresh-technology": ("technology",),
    "count-population": ("population",),
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
    "write-special-phenomena": ("phenomena-shells", "bright-stars-backfill"),
    "build-phenomenon-indexes": ("phenomena-special",),
    "record-phenomenon-scatter": ("phenomena-indexed",),
    "draw-nebula-field": ("skeleton", "nebula-field-cleared"),
    "scan-species": ("sectors",),
    "date-civilizations": ("species",),
    "refresh-technology": ("civilizations",),
    "count-population": ("technology",),
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


def requires_of(kind):
    """The kinds of command that make what command `kind` needs, in the order of `MAKES` (not `kind` itself): the kinds a stage of `kind` waits for when the run has them."""
    found = []
    for artifact in needs_of(kind):
        for maker in makers_of(artifact):
            if maker != kind and maker not in found:
                found.append(maker)
    return tuple(found)


def required_stage_keys(steps, index):
    """The keys of the earlier steps of `steps` that step `index` waits for: those of a kind in `requires_of` its kind."""
    wanted = set(requires_of(steps[index].kind))
    return [step.key for step in steps[:index] if step.kind in wanted]


def validate_order(steps):
    """
    Refuses a list of steps in which a need is met only after the step that has it: for each artifact a step needs that some step of the list makes, a maker comes before the step. A need no
    step of the list makes is a stored result of an earlier run. The wrong order is found where the list is made and not by a run that reads rows that are not there yet.

    Raises:
        ValueError: Naming the step, the artifact and the later steps that make it.
    """
    for index, step in enumerate(steps):
        for artifact in needs_of(step.kind):
            makers = [place for place, other in enumerate(steps) if place != index and artifact in makes_of(other.kind)]
            if makers and min(makers) > index:
                raise ValueError(f"step {step.key!r} ({step.kind}) needs {artifact}, which only {', '.join(repr(steps[place].key) for place in makers)} (later in the list) makes")
    return steps


RESETS = frozenset({"reset-galaxy", "plan-skeleton"})
"""frozenset: The commands that leave nothing a galaxy had made before them (`reset-galaxy` empties it, `plan-skeleton` plans it again with a new seed): `exists_after` starts again from their stage."""


def exists_after(stages):
    """
    The artifacts a galaxy has when `stages` (`(kind, made)` per finished stage, oldest first: the stages whose rows say they ended, with the artifacts their row says they made) are over. A stage of a command in `RESETS`
    first forgets everything the stages before it made.
    """
    found = set()
    for kind, made in stages:
        if kind in RESETS:
            found.clear()
        found.update(made)
    return found


def unmet(steps, exists):
    """
    `[(step, artifact)]`: each need of a step that no earlier step of the list makes and `exists` (the artifacts the ledger's rows say the galaxy has, `planrows.made_artifacts`) does not hold. A need that a later step of the list
    makes is `validate_order`'s to refuse; here it is unmet as well, since nothing has made it when the step starts.
    """
    lacking = []
    for index, step in enumerate(steps):
        have = set(exists)
        for other in steps[:index]:
            if other.kind in RESETS:
                have.clear()
            have.update(makes_of(other.kind))
        lacking.extend((step, artifact) for artifact in needs_of(step.kind) if artifact not in have)
    return lacking
