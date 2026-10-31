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
| `III` or `II` with no phase | `red-giant` |
| `III` or `II` with phase `horizontal-branch` or `core-helium-burning` | `helium-burning` |
| `III` or `II` with phase `asymptotic-giant` | `asymptotic-giant` |
| `III` or `II` with phase `classical-cepheid` | `cepheid` |
| `IB`, `IAB`, `IA` | `supergiant` |
| `IA+`, `0` | `hypergiant` |
| `W` (Wolf-Rayet) | `wolf-rayet` |
| `VII`, `D` | `white-dwarf` |
| a brown dwarf (no class `V`, a kind) | `brown-dwarf` |

**Age since it stopped burning** (`decay_age_yr`): defined where a star has stopped: a white dwarf's cooling age in years, the star's age less
the age at which its progenitor ended its life (`stellar_evolution.cooling_age_gy`, the age `evolve_star` cools the white dwarf from). A star
still burning anything (hydrogen in the core or a shell, helium) has none, and a brown dwarf never burned hydrogen; the sheet rule
(`describe/sheet.py` APPLIES) says so in a sentence. `set_state` fills both on a `Star`, from the stored class, type text, age, initial mass and
[Fe/H], so a star read back from storage and one just generated agree.
"""

from planetgen.physics import pre_main_sequence, stellar_evolution, wolf_rayet

STATES = ("protostar", "pre-main-seq", "main-sequence", "subgiant", "red-giant", "helium-burning", "asymptotic-giant", "cepheid",
          "supergiant", "hypergiant", "wolf-rayet", "white-dwarf", "brown-dwarf")
"""tuple: Every state, each within the 16 characters of the `stars.state` column."""

_BY_CLASS = {
    pre_main_sequence.PROTOSTAR: "protostar", pre_main_sequence.PRE_MAIN_SEQUENCE: "pre-main-seq",
    "V": "main-sequence", "VI": "main-sequence", "IV": "subgiant", "II": "red-giant",
    "IB": "supergiant", "IAB": "supergiant", "IA": "supergiant", "IA+": "hypergiant", "0": "hypergiant",
    wolf_rayet.CLASS: "wolf-rayet", "VII": "white-dwarf", "D": "white-dwarf",
}
_BY_PHASE = {"horizontal-branch": "helium-burning", "core-helium-burning": "helium-burning", "asymptotic-giant": "asymptotic-giant",
             "classical-cepheid": "cepheid"}


def phase_of(star_type):
    """The `stellar_evolution` phase a star's type text names (`"G8III Yellow Horizontal Branch Giant Star"` -> `"horizontal-branch"`), or
    None: the text is what a stored star keeps of the phase."""
    for phase, name in stellar_evolution.PHASE_NAMES.items():
        if f"{name} Star" in (star_type or ""):
            return phase
    return None


def evolutionary_state(yerkes_class, phase=None, brown_dwarf=False):
    """The state word of a star of this Yerkes class and phase (`stellar_evolution.evolve_star`'s `phase`, or None). `brown_dwarf` is a
    body below the hydrogen-burning limit, which has no luminosity class of the stars'. Raises ValueError for a class the table lacks."""
    if brown_dwarf:
        return "brown-dwarf"
    if yerkes_class in ("II", "III") and phase in _BY_PHASE:
        return _BY_PHASE[phase]
    try:
        return _BY_CLASS[yerkes_class] if yerkes_class != "III" else "red-giant"
    except KeyError:
        raise ValueError(f"no evolutionary state for the Yerkes class {yerkes_class!r}") from None


def decay_age_yr(state, cooling_age_gy):
    """The age, years, since the star stopped burning: a white dwarf's `cooling_age_gy`; None for every state that is still burning or
    never burned."""
    if state != "white-dwarf" or cooling_age_gy is None:
        return None
    return cooling_age_gy * 1.0e9


def set_state(star):
    """Sets `state` and `decay_age_yr` on a `Star` from its `yerkes_class`, `type`, `age`, `initial_mass_sol` and `metallicity_fe_h`
    (call once the [Fe/H] is set)."""
    star.state = evolutionary_state(star.yerkes_class, phase_of(star.type))
    cooling = None
    if star.state == "white-dwarf" and star.initial_mass_sol:
        cooling = stellar_evolution.cooling_age_gy(star.initial_mass_sol, star.age, star.metallicity_fe_h)
    star.decay_age_yr = decay_age_yr(star.state, cooling)
