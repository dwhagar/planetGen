# planetgen/queue/tokens.py

"""
Named resources (OPS.232, handle D7 of `reports/job-engine-items-stage-3-6-2026-11-13.md`; design sections 3.4 and 3.5): what a job must hold besides memory and cores.

A claim carries `tokens`, the names of the resources it must hold while it runs. The gate (`memorygate.decide`) lets a claim go only when no claim held already holds one of its tokens, so a token is
held by one job at a time and the job that waits for it is told which job holds it (the Queue page shows that as the reason). Three kinds of name:

* `ALONE` (`alone`): the host for itself. It is the claim `memorygate.Claim.solo` was (OPS.145, the restarted unit that runs alone): it asks for the whole memory budget of the bulk and maintenance
  lanes, waits for every claim of those lanes to end and is let by no other claim of them while it is held (`holds_host`). The site's lanes keep their places.
* `MYSQL_WRITE` (`mysql_write`): the one writer of the heavy kind. A command that writes so much that two at once only slow each other down (`CommandSpec.tokens`) holds it for the length of a unit.
* `warm_set(name)` (`warm:<name>`): the one build of a warm set at a time. A template process (OPS.226) builds a warm set (the bright-star sampling tables, 1.5 GB and 52 s); two builds of the same set
  at once would hold two copies, so the build holds the set's token.

A token is a name, so a new resource is a new name and no new rule. The rules the first two claims had (`Claim.exclusive` for the cut keeps its own: the whole lane to itself, which a token cannot say, as
it is the absence of every other claim of the lane) are not changed.
"""

ALONE = "alone"
MYSQL_WRITE = "mysql_write"
WARM_PREFIX = "warm:"


def warm_set(name):
    """The token of the warm set `name`."""
    return f"{WARM_PREFIX}{name}"


def holds_host(tokens):
    """Whether `tokens` hold `ALONE`: the claim runs alone (OPS.145)."""
    return ALONE in tokens


def clash(tokens, others):
    """
    The reason `tokens` (of a claim that asks) cannot be held while the claims `others` (`memorygate.Holder`s) hold theirs, or `None`. `ALONE` is not a token that clashes here (its rule is the host's:
    `memorygate.decide`); every other name is held by one claim at a time.
    """
    for token in tokens:
        if token == ALONE:
            continue
        for holder in others:
            if token in holder.claim.tokens:
                return f"waits for {holder.claim.kind} to release {token}"
    return None
