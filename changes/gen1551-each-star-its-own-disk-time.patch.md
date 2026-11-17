**Headline:** A young star now holds its disk for its own time, drawn from the surveys' exponential, so disks are found between 3 and 10 Myr and not every disk ends at 3 Myr.

### Changed

- GEN.1551 slice 1: a star's class 0 and I disk ends at the ages the stages have (0.15 and 0.5 Myr) and its class II disk at its own dispersal time (`pre_main_sequence.dispersal_yr`, from `disk_dispersal.draw_dispersal_time_yr`: an exponential of 2.5 Myr at 0.5 solar masses, 3 to 1 Myr by mass), in place of the fixed 3 Myr step (`tuning.DISK_CLASS_AGES_YR` no longer holds the II age). The time is drawn from the stream of the star's initial mass and no age, so it does not move when the age does, and it is held at `tuning.PLANET_MIN_STAR_AGE_GY` (10 Myr) so that no star holds a disk and planets together until GEN.1552 lets the star's own time set the planet cut. `Star` passes its initial mass to `draw_disk`. A 0.5 solar mass star at 5 Myr now has a class II disk; a Sun-like star's disk, 0.29 Myr in this draw, is gone at 0.5 Myr. The close-companion shortening (`CLOSE_COMPANION_TAU_MYR`) is in the law and not yet read: a star's companions are attached after its disk is drawn, so the system step that reads it is GEN.1475's to take.

generation-output: changed
generation-affects: scatter, fill
