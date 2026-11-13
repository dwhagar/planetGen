**Headline:** A star is called a subdwarf from the same temperature it is stored with, so a stored subdwarf always passes the subdwarf test.

### Changed

- `stellar_evolution.yerkes_class_of` (GEN.943) judged the subdwarf cut (0.3 dex below the Pecaut and Mamajek dwarf sequence) on the unrounded temperature, while `star_params` stores the temperature rounded to the nearest 100 K. A star near the cut could be labelled class VI and then fail the same test on its stored numbers (TEST.727, `test_a_subdwarf_is_an_old_metal_poor_star`). The cut now uses the rounded temperature, so the label and the stored star agree. A few stars close to the cut change between class V and class VI.

generation-output: changed
