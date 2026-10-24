# planetgen/generation/fill_stages.py

"""
The stages of a sector fill, each on its own random stream (GEN.251)
=====================================================================

`run_sector.generate_sector` used to draw everything from the one stream its sector was bound to, in a fixed order, so
adding a kind in the middle moved every draw after it. A fill is now a list of stages
(docs/design/systems/phenomena-catalogue-and-fill.md section 3), and each stage draws from a stream of its own, seeded from
the sector's stream and the stage's name:

====  ==============  ================================================================================
      stage           what
====  ==============  ================================================================================
A     pre-placed      the bright and massive stars placed by the scatter become full systems
B     scattered       the phenomenon scatter's rows become the objects they stand for
C     system-pass     the objects of B and G that can hold planets roll for them (`system_pass`)
D     systems         the cell's other systems
E     per-star-kinds  the kinds rolled per star (rogue planets, brown dwarfs, comets), one stream per kind
F     nebulae         the galaxy's cloud field, then the nebulae grown around the sector's hot stars
G     below-cut       the neutron stars and black holes below the scatter's mass cut
H     flags           runaway and hypervelocity flags, then the galactic nucleus (stream `flags:nucleus`)
J     ejecta          the neutron stars and black holes the galaxy ejected (GEN.238), one stream per kind; it
                      runs after G, and as the last stage added it moves no other stage's draws
====  ==============  ================================================================================

Stage I is not a stream: `count` keeps the objects each stage made by kind in `sector.stage_counts`, which the
representation audit of the catalogue's section 7 reads and `describe` prints.

A stage's stream is `draw.Stream(f"{base}:{stage}")` (`f"{base}:{stage}:{sub}"` for a stream within a stage), `base`
being one 64-bit number drawn from the sector's own stream when the fill begins, so two fills of one sector (the
galaxy binds the sector's stream from its address) draw the same, and a directive's redraws differ. Adding a kind to
the end of a stage, or a stage, therefore moves no draw of any other stage; the objects built from their own seeds
(A and B, and the system pass) were already independent of order.
"""

import contextlib

from planetgen.util import draw

PRE_PLACED = "pre-placed"
SCATTERED = "scattered"
SYSTEM_PASS = "system-pass"
SYSTEMS = "systems"
PER_STAR_KINDS = "per-star-kinds"
NEBULAE = "nebulae"
BELOW_CUT = "below-cut"
FLAGS = "flags"
EJECTA = "ejecta"

STAGES = (PRE_PLACED, SCATTERED, SYSTEM_PASS, SYSTEMS, PER_STAR_KINDS, NEBULAE, BELOW_CUT, EJECTA, FLAGS)
"""tuple: The stage names in the order a fill runs them (A to H, J between G and H)."""

BASE_BITS = 64
"""int: Bits of the number a fill's stage streams are seeded from."""


def begin(sector):
    """Starts the fill of `sector`: draws its stage base from the stream bound now and keeps the counts clean."""
    sector.stage_base = draw.getrandbits(BASE_BITS)
    sector.stage_counts.clear()
    return sector.stage_base


def stream_seed(sector, name, sub=None):
    """The seed of the stream of stage `name` (and `sub` within it), or `None` for a sector `begin` was not called
    on."""
    if sector.stage_base is None:
        return None
    return f"{sector.stage_base}:{name}" if sub is None else f"{sector.stage_base}:{name}:{sub}"


def stage(sector, name, sub=None):
    """
    A context manager that binds the draws to the stream of stage `name` (and `sub`) for the length of the block. A
    sector no fill has begun on (a caller that builds the pieces by hand) keeps drawing from the stream it has.
    """
    seed = stream_seed(sector, name, sub)
    return contextlib.nullcontext() if seed is None else draw.bound(seed)


def count(sector, name, kind, n=1):
    """Records that stage `name` made `n` more objects of `kind` (stage I)."""
    if n:
        sector.stage_counts[(name, kind)] += n


def counts(sector):
    """`{stage: {kind: n}}` of what the fill made, the stages in fill order."""
    by_stage = {name: {} for name in STAGES}
    for (name, kind), n in sorted(sector.stage_counts.items()):
        by_stage.setdefault(name, {})[kind] = n
    return {name: kinds for name, kinds in by_stage.items() if kinds}


def describe(sector):
    """One line per stage that made something: `"scattered: 2 black-hole, 1 neutron-star"`."""
    return [f"{name}: " + ", ".join(f"{n} {kind}" for kind, n in kinds.items())
            for name, kinds in counts(sector).items()]
