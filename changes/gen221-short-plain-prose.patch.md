**Headline:** The text the program writes about stars, worlds, gas giants, comets and life is shorter and plainer, and no longer claims life or "Earth-like" for a world the habitability index rates differently.

### Changed
- The star paragraph is two short sentences: "The star is about 13.8 Gyr old and should last about 1,321 Gyr." and one about life from `STAR_EVOLUTION`, whose notes are now full sentences with no claim about the universe ("the longest-lived star in the universe" is gone).
- The class texts no longer claim life, vegetation, a percentage of water or "Earth-like": F is "a volcanic world with shallow seas", G "a rocky, barren world", H "a desert world with little liquid water", L "a marginally habitable world", M "a temperate terrestrial world", O "an ocean world with most of its surface covered in liquid water". Because the text is the same whatever the index allows, `NO_LIFE_CLASS_DESCRIPTIONS`, `MICROBIAL_ONLY_CLASS_DESCRIPTIONS` and the moon re-wording in `life._set_description` are gone; the life sentences say what the index allows.
- The life sentences read "Life based on melanin could arise on this planet, and it would evolve at a normal pace." and "In visible light that life would look charcoal-black." (it was "look none (charcoal-black)"); the gas giant depth sentences drop the polytrope jargon; the evolutionary timeline is written in the describer from the stored state with its ages through the Gyr ladder (it said "13.80 Billion Years"); comet and binary sentences are shorter.
- New tests: the star sentence, the timeline and the class texts claim nothing the index can contradict.

TODO items resolved: GEN.221
