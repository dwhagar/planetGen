**Headline:** The search page writes every gas formula with its subscripts, in the gas rows, the sort menu and the filter chips.

### Changed
- **Formulas on /search (TEST.555 part 1).** The physical-values form of the search page (UX.227) showed `Partial pressure of O2` and `CO2` as plain text, so `test_no_page_shows_a_formula_that_is_not_written_out` failed on the page. A gas row's label and a filter chip now go through the `chem` filter (`O<sub>2</sub>`), and the text of the sort menu's and the "another gas" menu's `<option>`s, which cannot hold tags, carries Unicode subscripts (`O₂`, `C₆H₆`) through the new `chem_text` and `formula_text` filters; the values the form sends stay `O2` and `C6H6`. A test checks the rows, the menus and the chips. No stored value changes.
generation-output: unchanged
