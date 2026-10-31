**Headline:** The hottest a young brown dwarf can be now follows the hydrogen-burning limit at its age (up to 2,950 K) instead of a flat 2,800 K.

### Changed

- `brown_dwarf` and `rogue_surface` cap a body's temperature at the BHAC15 temperature of the hydrogen-burning limit at its age (`substellar_ceiling.cap_k`): 2,950 K at the peak, falling to about 1,630 K at 5 Gy for a body under 0.07 solar masses, in place of the 2,800 K every age shared. The row g1 of the atmosphere recipes reaches 2,950 K. The constant `tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K` is removed. generation-output: changed
