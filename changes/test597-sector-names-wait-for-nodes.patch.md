**Headline:** The sector page label test waits for the sector's map blocks to arrive before it checks which names are shown, so it no longer fails 1 run in 3.

### Fixed
- `test_sector_map_names_its_systems_through_the_label_engine_with_the_system_names_switch` (TEST.597) read the labels as soon as a system name showed, while the page was still loading its map blocks one level at a time (243, 81, 27, 9). A block is drawn whole until its children are in, and a whole block names its region or family for the two to four seconds that takes, so "Fixture Reach" showed in every run for a moment and the test failed whenever it read the labels in that moment. It waits until those names have gone. The map is as designed; no map bug was found.
