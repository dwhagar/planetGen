**Headline:** class texts now agree with what the classifier and the generated bodies do: dead world, hot dense, anoxic, carbon dioxide, marginal thin-air, small hot and hot thin-air worlds, plus two stale docs and a stale class-page note.

### Changed

- `physics/data/planet_classes.json`: class C no longer says "bare rock" only (about half the stored dead worlds carry ice or an ice-covered sea); class N says its air is usually mostly carbon dioxide, often with water vapor, not "carbon dioxide and sulfur" (sulfur is 3 percent of the stored air); classes F and G no longer say "mild" or "warm but not hot" (the conditions run from cool to 400 K); class L no longer says its air "is not toxic" (its carbon dioxide runs to tens of kilopascals); class B no longer says "not an ocean world" (nothing in its conditions bars an ocean); class Z no longer says "strong gravity" (its conditions start at half a gravity). The class page, planet page and system map read these words, so nothing stored changes.
- Boss (9:01 PM EST): the code bugs stay as they are for now and the shown text follows what the code does today. T (terrestrial) says its air is mostly hydrogen and can be thin; U, Y and A say the liquid or sulfur dioxide could stand or is present, not that seas, acid or volcanoes are stored; X says it has no air; O, H, P and R say they have air you could breathe (their conditions ask for oxygen of 16 to 35.5 kPa); S says some are under ice.
- The class page note on rogue planets no longer says classes "list interstellar space" (no class does); a rogue planet is classed by the same conditions.
- `docs/reference/system-file-format.md`: the class table, the check order, the "Life?" column, the class-moves-the-orbit claim and the list of unused letters described the first-match classifier of October; they now describe the best-fit classifier and the 25 letters in use. `docs/design/habitability/condition-based-planet-classes.md`: the rule-file paragraph, class L's name, class B's rule and the reserve line are corrected and R, U, X and Y are added to the table.
- Report: `reports/class-descriptions-vs-code-2026-11-12.md` lists the mismatches that need a code decision.

generation-output: unchanged
generation-affects: none
