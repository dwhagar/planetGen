# planetgen/physics/capture_windows.py

"""
The stretches of a system's life in which a planet could catch small bodies (GEN.1082, docs/design/systems/captured-asteroid-moons.md sections 3 and 9)
================================================================================================================================================

A planet does not pick up the bodies of a belt that lies beside it (`irregular_population`): it caught the irregular moons it has in short windows in which something took the bodies' energy
away. Two of the windows are properties of the system, not of one planet, and the planet architecture (`planet_architecture`, GEN.728) already holds what they need:

* **The instability.** Two giants that end close enough to scatter each other (`planet_architecture.Lost.thrown`: a pair of neighbours within `GIANT_CLEARING_SPACING` mutual Hill radii that holds
  a giant, the lighter thrown out) pass each other many times before one leaves. Each such event with a giant on both sides is one instability, with `ENCOUNTERS_PER_INSTABILITY`
  close giant-giant encounters, the number a reference instability had in the three histories of [R2158] (`irregular_population.REFERENCE_ENCOUNTERS`). A system with no such pair had none,
  and its planets get no captures from this window. The encounters are the system's, one number: the planets that stayed share them (the Hill volume of each decides how much it caught,
  `irregular_population.make_pool`). The same event tells the scatter step of the outer disk (`scatter.scatter_shares`) how much the giants stirred each other: `instability_level`.
* **The gas era.** The disk's gas lasted the circumplanetary disk's lifetime that each giant already keeps (`moon_composition.DiskSystem.lifetime_myr`, stored in the planet's `moon_disk_lifetime_myr`): the
  gas window needs no new number, and `make_pool` reads it from the planet.

The binary-exchange window and the temporary one are the planet's own (`irregular_population.draw_exchange`, `TEMPORARY_MIN_PLANET_RADII`). Pure functions of the architecture; nothing is drawn.

Use:

    windows = from_architecture(architecture)      # Windows(instabilities, encounters)
    instability_level(windows)                     # 0 to 1, for scatter.scatter_shares
"""

from collections import namedtuple

ENCOUNTERS_PER_INSTABILITY = 130.0
"""float: The close giant-giant encounters of one instability: the mean of the three histories of [R2158] (59, 280, 49 for Jupiter), the reference instability of `irregular_population.REFERENCE_ENCOUNTERS`, which the irregular-moon wiring reads and this must equal. [Stand-in CW1, tier T2, tweak L1.]"""

FULL_INSTABILITIES = 1
"""int: The number of instabilities that count as a full one for the scatter step (`instability_level` is 1): the reference instability the scatter step's hot population is scaled to is a single event. [Stand-in CW2, tier Marginal educated guess from needed output, tweak L3.]"""

Windows = namedtuple("Windows", ("instabilities", "encounters"))
"""The system's capture windows: how many giant-giant instabilities the architecture holds and the close giant-giant encounters they come to."""

NO_WINDOWS = Windows(0, 0.0)


def from_architecture(architecture):
    """The `Windows` of an architecture: each pair in `lost.thrown` with a giant on both sides is one instability."""
    events = sum(1 for lighter, heavier in architecture.lost.thrown if lighter.runaway and heavier.runaway)
    return Windows(events, events * ENCOUNTERS_PER_INSTABILITY)


def instability_level(windows):
    """0 to 1, how much the giants of the system scattered each other: the instabilities over `FULL_INSTABILITIES`, at most 1."""
    return min(windows.instabilities / FULL_INSTABILITIES, 1.0)
