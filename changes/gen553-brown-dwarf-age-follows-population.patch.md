**Headline:** A free-floating brown dwarf is as old as the stars of the place it is in, not any age from 0.1 to 10 billion years.

### Changed
- A brown dwarf made by a galaxy sector's fill takes a population from the sector's mix, like a star, and its age comes from that population's age window spread by its star formation rate (GEN.546's `age_law`). A brown dwarf in the thin disk's young population is under 100 million years old; one in the halo is 10 to 13 billion. Its radius, temperature and spectral type follow from the age as before.
- A brown dwarf made without a place (the single-object generator) draws its age over the whole disk's window, 0 to 10 billion years, in place of 0.1 to 10.
- `FillContext.apply_population` is the part of `FillContext.apply` that gives a body its population, age law, [Fe/H] peaks and place terms. First step of GEN.553; the rest of the body draw follows.
- `docs/design/brown-dwarf-state.md` describes the rule.

generation-output: changed

TODO items resolved: none (GEN.553 stays open) / Added: none
