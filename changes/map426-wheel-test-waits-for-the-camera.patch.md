**Headline:** The browser test that carries the Galaxy Map from the whole galaxy into a system and back out waits for the camera to come to rest after each wheel tick, so it no longer depends on the speed of the machine.

### Changed

- `test_the_wheel_carries_the_map_from_the_whole_galaxy_into_a_system` (MAP.426) sent a wheel tick every 120 ms whatever the camera was doing. A tick is one zoom step and the page ignores a tick while a stage change is flying (`zoomAcross`), so on a slow machine the ticks were dropped or piled up in a zoom goal and the run stopped at a different place each time (stuck 6 pc from the star, past the star without opening it, or still in the sector after 400 ticks out). The test now waits for the camera distance to repeat (`_at_rest`) after each tick, going in and coming out. The page is unchanged.

generation-output: unchanged
