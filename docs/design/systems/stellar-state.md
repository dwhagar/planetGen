# Evolutionary state and age since it stopped burning (GEN.786 and GEN.785, definitions)

Status: built and wired (GEN.785 and GEN.786 wiring, rides a train: it fills two columns of every star). Code: `physics/stellar_state.py`. Tests: `tests/test_stellar_state.py`. No stand-in values: the module holds words and one subtraction.

## State (GEN.786)

`evolutionary_state(yerkes_class, phase, brown_dwarf)` gives one of thirteen words, each within the 16 characters of `stars.state`: protostar, pre-main-seq, main-sequence, subgiant, red-giant, helium-burning, asymptotic-giant, cepheid, supergiant, hypergiant, wolf-rayet, white-dwarf, brown-dwarf. A subdwarf (VI) is a main-sequence star of low metals. A class III or II star with phase `horizontal-branch` or `core-helium-burning` is helium-burning, with `asymptotic-giant` asymptotic-giant and with `classical-cepheid` a cepheid (the phase names are `stellar_evolution`'s and are too long for the column). A stored star keeps its phase only in its type text (`PHASE_NAMES`), so `phase_of` reads it back from there. A class the table lacks raises. The test sweeps `evolve_star` over masses and ages so a new class or phase the evolution gains fails the test, not the sheet.

## Age since it stopped burning (GEN.785)

`decay_age_yr(state, cooling_age_gy)`: for a white dwarf, the cooling age in years, which is the star's age less the end of its progenitor's life (`stellar_evolution.cooling_age_gy`, the age `evolve_star` cools it from, at least `WD_MIN_COOLING_AGE_GY`); none for every state still burning and for a brown dwarf. The sheet rule says so: the row applies to a white dwarf only. A neutron star's or black hole's own age follows the remnant-age note (`remnant-age.md`).

## Wired

`Star._finish_init` calls `stellar_state.set_state(star)` (no draw) after the [Fe/H]; `store.STAR_FIELDS` writes `state` and `decay_age_yr`; the sheet's open gaps for `state` and `decay_age_yr` are gone and `describe/sheet.py` APPLIES gives `decay_age_yr` to white dwarfs only. `stars.state` is the evolutionary word; `scatter_members.state` (a placed member's Yerkes class, `multistar.Member`) is a different vocabulary in a different table. `lethal_event_rate_per_gyr` (the rest of GEN.785) is not touched here. Output changes (two new columns on every star), so it rides a train: `--affects fill`, migration no (both columns exist).
