**Headline:** An unused function that drew the composition of ring material from the planet's disk is removed; rings take their ice share from `rings.build_ring` as before.

### Changed

- `small_body_composition.draw_ring` had no caller since the rings were rebuilt from history (GEN.1057): a ring's ice share is drawn by its kind in `rings.build_ring` (broad rings 90 to 99 percent ice like Saturn's, the stripped mantle of a differentiated moon), not from a distance in the circumplanetary disk. The function, its test and its wiring-register line are gone; `disk_ice_fraction`, which the moons read, stays. No generated value changes.

generation-output: unchanged
generation-affects: none
