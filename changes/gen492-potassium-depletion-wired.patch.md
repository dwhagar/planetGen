**Headline:** The potassium a rocky planet, moon or dwarf planet keeps is now drawn from the measured K/U of six Solar System bodies, by the body's mass, instead of one flat range.

### Changed

- GEN.492 (GEN.1559 wiring): `radiogenic.draw_content` takes the potassium depletion from `potassium_depletion.draw_depletion` with the body's mass and the CI K/U from `potassium_depletion.ci_k_over_u` (71,429, [R0270]). A planet of 0.03 Earth masses or more is drawn about 0.19 of the CI value (0.15 dex), a smaller body about 0.015 (0.94 dex), never above 1. Before, every body drew a depletion log-uniform from 0.03 to 1 on a K/U of 67,407. K-40 heat and the ground dose of a body change with it, and the body's stored potassium moves, so a galaxy is made again at the epoch of this change. The draw count is the same (one draw), so the later draws of a body keep their place.
- `radiogenic.draw_content` takes the body's mass in Earth masses as its first argument.

### Removed

- `tuning.RADIOGENIC_K_DEPLETION_RANGE` and `radiogenic.CI_CHONDRITE_K_OVER_U` (no shim). The wiring guard no longer lists `potassium_depletion`, its data file and GEN.725 rows KD3 and KD4 as unread.

generation-output: changed
generation-affects: fill
