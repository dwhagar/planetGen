**Headline:** The flavor lines of planets, moons and systems are now a gated tree in a data file, walked so only the lines that fit the body are offered.

### Changed
- **The flavor text is a tree (GEN.1430).** `physics/data/flavor_tree.json` holds the 82 lines of the four old lists in two scopes (system, planet), with gates on stored conditions (a ground to stand on, multicellular or technological life) and a tag for the ancient works family; `physics/flavor_tree.py` walks it with the gate evaluator that `planet_class_fit.json` and `oddities.json` already use, takes out the lines the system already showed (this replaces the recency list), and draws one by weight on the body's own stream, so the line no longer depends on the order the bodies were made in. The four lists in `tuning.py` and `life.flavor_key` are gone; stored rows keep their keys, which are the line ids. A world with multicellular life is now offered the life lines and the ground lines together. This changes the generator output, so it takes a generator epoch (the Train Manager assigns it).

generation-output: changed
generation-affects: fill
