**Headline:** The program can now work out how far a planet's tilt wanders under the pull of its star, moons, sibling planets and a companion star; nothing reads it yet, so no planet changes.

### Changed

- GEN.1454 (part 1): new `physics/obliquity_variation.py`. From a planet's spin, its moons, the star and a companion it works out the precession of the spin axis; set against the nodal modes of the secular exchange (GEN.1450) it gives the range of the planet's obliquity, whether the tilt sits in a resonance and whether two resonances overlap into a chaotic zone. The Earth with its Moon comes out at 22.6 to 24.3 degrees against the 22.1 to 24.5 measured, and without the Moon it swings about five times as far. It is a lower bound for the chaos of the full theory. Registry section 6.271 (OB1 to OB5). Part 2 will use it in the settle pass, the seasonal swing and the oddity table.

generation-output: unchanged
generation-affects: none
