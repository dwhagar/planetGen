**Headline:** The test that opens a sector by zooming in no longer fails when the camera comes up to a neighbouring sector.

### Fixed
- **A flaky browser test (TEST.419).** `test_opening_a_sector_is_one_continuous_view_with_a_filled_panel` zoomed in on one
  sector and expected exactly it to open, but the view opens the nearest charted sector the camera comes up to (MAP.152),
  and the sample database charts sectors side by side at the core, so a run where the view eased a little off the aimed
  centre opened the neighbour. The test now holds the address and the panel to the sector that opened.
