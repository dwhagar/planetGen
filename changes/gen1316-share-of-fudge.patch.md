**Headline:** the class pages say what the fudge is and how soft each class edge is, and a new ranking lists every stand-in by its share of the width of each class quantity.

### Changed

- `/classes/planet` says that the width of a class edge is the model's doubt about the value (what the stand-ins that make it could move it by), not the variety of the worlds; each planet class page gains a "Soft edges" line ("mean surface temperature at least 275 K: soft by 0.925 K"), read from the fudge table through the new `class_fit.softness_words`. `planet_classes.quantity_entry` is the label, unit and percent flag that `condition_words` already used, now shared.
- `scripts/build_fudge_ranking.py` builds `docs/design/research/stand-in-fudge-ranking.md` from the fudge table, its matrix and the registry (reusing `build_fudge_table.column`): the quantities with what makes their width, and every stand-in with at least 1 percent of any quantity's variance, its tier, tweak and unreliability score, its shares and how far the width falls if it were exact. A test fails when the page differs.
- `condition-based-planet-classes.md` gains the section "What the fudge is (GEN.1316)".
- No stored value changes; the class table and the fudge table are untouched.

generation-output: unchanged
generation-affects: none
