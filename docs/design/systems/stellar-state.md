# Evolutionary state and age since it stopped burning (GEN.786 and GEN.785, definitions)

Status: module built, nothing in generation calls it (the wiring rides train 3 after Star 1's GEN.935). Code: `physics/stellar_state.py`. Tests: `tests/test_stellar_state.py`. No stand-in values: the module holds words and one subtraction.

## State (GEN.786)

`evolutionary_state(yerkes_class, phase, brown_dwarf)` gives one of twelve words, each within the 16 characters of `stars.state`: protostar, pre-main-seq, main-sequence, subgiant, red-giant, helium-burning, asymptotic-giant, supergiant, hypergiant, wolf-rayet, white-dwarf, brown-dwarf. A subdwarf (VI) is a main-sequence star of low metals. A class III star with phase `horizontal-branch` is helium-burning and with `asymptotic-giant` is asymptotic-giant (the phase names are `stellar_evolution`'s and are too long for the column). A class the table lacks raises. The test sweeps `evolve_star` over masses and ages so a new class or phase the evolution gains fails the test, not the sheet.

## Age since it stopped burning (GEN.785)

`decay_age_yr(state, age_gy, lifespan_gy)`: for a white dwarf, the star's age less the age at which it became one (`lifespan_gy`), in years and never negative; none for every state still burning and for a brown dwarf. A neutron star's or black hole's own age follows the remnant-age note (`remnant-age.md`).

## Wiring

`Star.__init__` sets `state` and `decay_age_yr` from the evolved star; the sheet's `KNOWN_EMPTY` entries for `state` and `decay_age_yr` are removed, and `lethal_event_rate_per_gyr` stays with its own item. The wiring waits for GEN.935, which adds the post-main-sequence phases this table must cover (a new phase gets one line in `_BY_PHASE`).
