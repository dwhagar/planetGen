**Headline:** The hydrosphere test no longer fails on a halite sea: the water activity a chloride brine can hold now includes the saturated rock-salt solution.

### Fixed

- **`test_generated_bodies_store_a_consistent_hydrosphere` failed on main since GEN.1730.** A halite sea's water activity is the saturated rock-salt solution's 0.753 (set, not drawn), above the 0.30 to 0.65 chloride-brine range the deliquescent brines are drawn in. The test now checks a liquid body against the range it can hold (its class's drawn range; only a chloride brine whose activity is exactly `NACL_SATURATION_ACTIVITY`, a saturated halite sea, may reach that figure; an ordinary chloride brine stays 0.30 to 0.65). No product code changes, no draw and no stored value changes.

generation-output: unchanged
generation-affects: none
