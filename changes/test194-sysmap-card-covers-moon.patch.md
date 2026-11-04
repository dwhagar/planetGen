**Headline:** The System Map drill test no longer fails when the planet's info card sits over the moon it clicks next.

### Changed

- `test_system_map_selection_drill_and_measure` expanded the planet's info card to read it, then clicked the first moon. The expanded card is as tall as the window, and the system is a random one, so when the moon lay under the card the click landed on the card and the test waited out its 120 second timeout (TEST.194, seen as a load flake; it failed 1 of 3 runs alone). The test now presses Expand again to shrink the card before it clicks the moon. Four runs alone pass.
