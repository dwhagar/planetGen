**Headline:** The hydrosphere test no longer fails on a halite sea: the water activity a chloride brine can hold now includes the saturated rock-salt solution.

### Fixed

- **`test_generated_bodies_store_a_consistent_hydrosphere` failed on main since GEN.1730.** A halite sea's water activity is the saturated rock-salt solution's 0.753 (set, not drawn), above the 0.30 to 0.65 chloride-brine range the deliquescent brines are drawn in. `hydrosphere.ocean_activity_range(body)` gives the range a body can hold (its class's drawn range, and only a chloride brine that is a saturated halite sea, activity exactly `NACL_SATURATION_ACTIVITY`, up to that figure) and the test checks against it. An ordinary chloride brine stays held to 0.30 to 0.65. No draw and no stored value changes.

generation-output: unchanged
generation-affects: none
