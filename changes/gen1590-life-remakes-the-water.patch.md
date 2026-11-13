**Headline:** A world whose air life changes now has water that fits the air it ended with.

### Fixed
- **GEN.1590: life changed the air and left the water as it was.** After life moves a body's surface temperature or pressure, the hydrosphere, the ocean cover and chemistry, and the ground are made again through `planets.refresh_at_final_age` (the step the thaw path already used; no second water model), on the body's own stream, before the cloud, dose and class are made again. A body whose temperature and pressure life left alone keeps what it had. The thaw path now makes the class again after it remakes the ground, so the class no longer reads the ice cover it had. In 160 census systems no ocean is stored above its boiling point afterwards (2 before), 27 hydrospheres change (hycean to surface ocean 10, surface ocean to ice-covered ocean 9, hycean to ice-covered ocean 8) and 14 of 1,541 bodies change class.
generation-output: changed
generation-affects: fill
