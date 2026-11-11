**Headline:** The test that checks the habitability gear sentence against the equipment line now compares them the way the page writes them, formulas set properly.

### Changed

- TEST.155: `test_every_scored_body_has_both_equipment_lines_and_its_per_item_specs_in_its_details` compared the Document's gear sentence, where a formula is written as a formula (UX.149: `CO₂`), with the equipment line's plain letters (`CO2`). It passed while none of the first 8 seeds needed a CO2 scrubber and failed once Train 21's draws gave a world that does (Illughoris II). The line and the item names now go through `model.chem_text` first, as the page does. No code under test changed.

generation-affects: none
