**Headline:** Galaxy Map: the wheel reaches any star you can pick, and opens its system

### Fixed

- MAP.434: the wheel stalled short of the star it was aimed at. At the floor of the zoom range (the view's distance cannot shrink any more) a step moved nothing, because the camera was carried toward the pointed-at thing only by the share the zoom shrank. The share the floor takes from a step is now made up by moving the camera along the ray to that thing (it stays under the pointer), so the wheel arrives at the star. A faint, drawn, pickable star stopped the wheel 9.4 pc short of it before.
- MAP.435: the wheel could not open the system of a drawn star that lies in a block next to the open sector (the open sector held none of its stars, so there was nothing to enter). The wheel now enters the system of the drawn star under the pointer when the open sector has none near it (`drawnSystemNear`, and `host.systemHref` for its page).

### Changed

- `test_the_wheel_carries_the_map_from_the_whole_galaxy_into_a_system` aims at a star the map draws and can pick (`galaxyPickableStars`), followed by its system's id as the camera comes in, and skips with its reason when the map never draws that star; the system-handover tests wait for the camera to rest between ticks, as the wheel's steps now carry on at the floor and ticks sent faster than the camera outrun the system scene's load.
- `test_web_browser_fixture_maps.py` reads the address without `whole=1` (the flag the tests add to open the whole galaxy, MAP.189, which the address keeps): four tests of the address (`window_stage_zooms_in_and_back_out`, `opens_a_sector_in_place`, `opens_an_uncharted_sector...`, `back_forward_and_url_state`) were red on main since MAP.189 for that one flag.
