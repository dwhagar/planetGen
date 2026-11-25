**Headline:** The Galaxy Map opens flown in to the largest charted region, not the whole galaxy

### Changed

- MAP.189 (GitHub issue 1232): with sectors charted, the first frame of the Galaxy Map flies in to the largest charted region (the first of the charted-regions list: its window, the camera where the region's sphere fills the view, as a click on the region lands, MAP.319) and shows only that space at the start; the wheel or Whole galaxy still reaches the whole galaxy. An address that asks for something else (a camera, a sector, a system) is left as it was, and `?whole=1` opens the whole galaxy as before; the whole galaxy's own address carries `whole=1` when the map has a region to fly to, so a fresh load of an address shows the view it names (MAP.106). Browser test `test_the_galaxy_map_opens_flown_in_to_the_largest_charted_region_and_zooms_out_to_the_whole_galaxy`; the browser tests that look at the whole galaxy ask for `?whole=1`.
