**Headline:** The Conditions line of a planet class page shows its pressures, temperatures, gravities and masses in the visitor's unit setting instead of raw pascals and kelvins.

### Changed

- The class page's Conditions line (and the "Also" line of a class with two rows) wrote each number with the unit the rules are kept in, so a visitor who had chosen another unit system still read pascals, kelvins and Earth masses. Each number that has a unit is now a `qty` span with its value in SI (`web/lib/classref.py`, `fmt.measured_html`), the markup `static/units.js` shows in the chosen units; with the default setting the span reads as the site's ladder writes it (16 kPa, 233 K, 19.6 m/s²). A threshold is shown to three figures there (233.15 K as 233 K). The words round the numbers are escaped as before. `planet_classes.condition_words` and `class_fit.rule_words` take an optional `measure` hook for this; without it they write the same text as before.

generation-output: unchanged
