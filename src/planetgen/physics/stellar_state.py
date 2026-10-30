# planetgen/physics/stellar_state.py

"""
The evolutionary state of a star and the age since it stopped burning (GEN.786 and GEN.785, definitions).

The star sheet carries two columns nothing fills: `state` ("State", wanted by the sheet and by the post-main-sequence work of GEN.935) and
`decay_age_yr` ("Age since it stopped burning"). This module defines them from what the star already stores, as pure functions, so the
wiring (a call in `Star.__init__` once Star 1's GEN.935 phases are in the star epoch, then the sheet's `KNOWN_EMPTY` entries removed) is one
line each.

**State** (`evolutionary_state`), a short word that fits the 16-character column:

| Yerkes class and phase | State |
|---|---|
| `protostar` | `protostar` |
| `PMS` (pre-main-sequence) | `pre-main-seq` |
| `V`, `VI` (a subdwarf is a main-sequence star of low metals) | `main-sequence` |
| `IV` | `subgiant` |
| `III` with no phase, `II` | `red-giant` |
| `III` with phase `horizontal-branch` | `helium-burning` |
| `III` with phase `asymptotic-giant` | `asymptotic-giant` |
| `IB`, `IAB`, `IA` | `supergiant` |
| `IA+`, `0` | `hypergiant` |
| `W` (Wolf-Rayet) | `wolf-rayet` |
| `VII`, `D` | `white-dwarf` |
| a brown dwarf (no class `V`, a kind) | `brown-dwarf` |

**Age since it stopped burning** (`decay_age_yr`): defined where a star has stopped: a white dwarf's cooling age, the star's age less the
age at which it became a remnant (`lifespan_gy`; `physics/remnant_age.md` gives a neutron star's or black hole's own age the same way). A star
still burning anything (hydrogen in the core or a shell, helium) has none, and a brown dwarf never burned hydrogen; the sheet shows the row
only where it is not none. The sheet rule (`describe/sheet.py` APPLIES) says so in a sentence when the column is wired.
"""

from planetgen.physics import pre_main_sequence, wolf_rayet

STATES = ("protostar", "pre-main-seq", "main-sequence", "subgiant", "red-giant", "helium-burning", "asymptotic-giant", "supergiant",
          "hypergiant", "wolf-rayet", "white-dwarf", "brown-dwarf")
"""tuple: Every state, each within the 16 characters of the `stars.state` column."""

_BY_CLASS = {
    pre_main_sequence.PROTOSTAR: "protostar", pre_main_sequence.PRE_MAIN_SEQUENCE: "pre-main-seq",
    "V": "main-sequence", "VI": "main-sequence", "IV": "subgiant", "II": "red-giant",
    "IB": "supergiant", "IAB": "supergiant", "IA": "supergiant", "IA+": "hypergiant", "0": "hypergiant",
    wolf_rayet.CLASS: "wolf-rayet", "VII": "white-dwarf", "D": "white-dwarf",
}
_BY_PHASE = {"horizontal-branch": "helium-burning", "asymptotic-giant": "asymptotic-giant"}


def evolutionary_state(yerkes_class, phase=None, brown_dwarf=False):
    """The state word of a star of this Yerkes class and phase (`stellar_evolution.evolve_star`'s `phase`, or None). `brown_dwarf` is a
    body below the hydrogen-burning limit, which has no luminosity class of the stars'. Raises ValueError for a class the table lacks."""
    if brown_dwarf:
        return "brown-dwarf"
    if yerkes_class == "III":
        return _BY_PHASE.get(phase, "red-giant")
    try:
        return _BY_CLASS[yerkes_class]
    except KeyError:
        raise ValueError(f"no evolutionary state for the Yerkes class {yerkes_class!r}") from None


def decay_age_yr(state, age_gy, lifespan_gy):
    """The age, years, since the star stopped burning: for a white dwarf its age less `lifespan_gy` (the age at which it became one),
    never negative; None for every state that is still burning or never burned."""
    if state != "white-dwarf" or lifespan_gy is None:
        return None
    return max(age_gy - lifespan_gy, 0.0) * 1.0e9
