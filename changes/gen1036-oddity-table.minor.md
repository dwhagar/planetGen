**Headline:** The weird-worlds oddity table and its loader exist, empty for now, with a weirdness dial and a stream of its own for every oddity roll.

### Added
- **The oddity table (GEN.1036).** `physics/data/oddities.json` holds the vocabulary for rare, physically plausible outcomes (the lever each one changes and the inputs it may touch, and the quantities a gate or a limit may read) and an empty list of rows; `physics/oddities.py` loads it and refuses a row without a source, a gate, limits, a rate, a tier or a tweak level. Nothing calls it yet: the weird-worlds modules add rows and GEN.1047 wires them into the generator, so no world changes.
- **A weirdness dial (GEN.1036).** `tuning.WEIRDNESS_DIAL` (1.0; 0 turns every oddity off) and `tuning.WEIRDNESS_RATE_SCALE` (a scale by oddity id) multiply a row's rate, and every roll comes from a stream keyed on the oddity and the body, so adding, removing or turning off one oddity moves no other draw.
