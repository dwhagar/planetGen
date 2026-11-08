**Headline:** The Galaxy Map drill-down browser test no longer reads the map before a breadcrumb click has taken effect.

### Fixed
- `test_galaxy_map_drill_down_by_windows` clicked a breadcrumb and waited only for the map to hold still. The camera leaves for the clicked step once its data has loaded, so on a loaded machine the wait returned on the unchanged window and the test read four steps where two were asked for. The wait now also needs the state to move off what it was before the click, as every other step of the test does.

generation-output: unchanged
