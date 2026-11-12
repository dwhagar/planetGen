**Headline:** The solar particle event counts typed in the code now match the 265 events of the NOAA list the setup downloads.

### Fixed
- **GEN.805: `solar_particles.NOAA_CUMULATIVE_COUNTS` counted from the downloaded list.** The typed counts (255, 177, 109, 72, 40, 18, 9) came from an older 255-event table; the list that `datasets.noaa_proton_events` reads has 265 events and counts 265, 181, 112, 74, 41, 18, 9. A new test holds the typed numbers to the file. The dataset is marked by design in the wiring register as the calibration reference of the event rate. Only the test-side rate helpers read these counts, so no generated output moves.

generation-affects: none
