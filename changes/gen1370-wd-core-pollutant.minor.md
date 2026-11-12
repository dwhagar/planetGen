**Headline:** The oxygen and neon in a white dwarf's core and the elements a polluted white dwarf takes from an accreted asteroid can now be worked out in one place; nothing reads them yet.

### Changed

- GEN.1370 (slice 0, N16 of the 8.1 real-material set): new `physics/white_dwarf_core.py`, pure. `central_oxygen` gives the central oxygen mass fraction of a carbon-oxygen white dwarf from its mass (the ten models of [R5238], linear between and held outside, kept inside the range two pulsating white dwarfs fit [R5239]); `neon22_fraction` gives the 22Ne from the carbon, nitrogen and oxygen of the progenitor's birth composition (0.0152 at the proto-solar mix); `core_fractions` gives the three central fractions, adding to 1. `pollutant` takes the stored materials of a belt body and returns its element mass fractions, the ratios to magnesium of the elements a white dwarf's spectrum shows (O, Mg, Si, Fe, Ca, Al, Ti, Cr, Mn, Ni, Na, K, S, C, P) and the class of the fragment: core-rich, crust-rich, volatile-rich or primitive [R1832]. Stand-ins WDC1 to WDC8 (tiers T2 to T6, tweak L0 and L1) are in section 6.285 of the gap-filling plan; the alloys, organics and six solid-solution minerals of the materials table that state no formula take the shares and end members of `tuning.py`.
- Nothing in the generator reads it yet (wiring register): storing the core fractions in the white dwarf's layers and the pollutant in its pollution detail is the lock wiring of GEN.1370.

generation-output: unchanged
