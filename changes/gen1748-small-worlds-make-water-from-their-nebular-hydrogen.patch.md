**Headline:** Rocky worlds from 0.2 to 1 Earth mass now capture a hydrogen envelope in the disk and make water from it, so their seas hold more.

### Changed

- GEN.1748 (every galaxy): a rocky body of 0.2 to 20 Earth masses draws a nebular envelope (`envelope_capture.MIN_MASS_EARTH` 0.2, from Young et al. [R5286], which finds the result independent of mass above about 0.2); under 1 Earth mass the median start is flat at 0.2 percent of the mass (it fell as the cube of the mass before, with none under 1). A body under 1 Earth mass is bare at the end (`KEEPS_FROM_MASS_EARTH`), its envelope lost with the disk, and the hydrogen it held made water in its magma ocean: `bulk_draw._nebular_water` adds `envelope_capture.nebular_water_fraction` (about 0.9 Earth oceans per Earth mass at the median start, 0 to 5 over the draw) to the delivered water of rocky bodies of 0.2 to 1 Earth mass that are not moons, late cores or belt members, and not where the Class M abundance tier's water dial decided the class. The start comes from the planet's own stream, so no draw of the run's moves.

generation-output: changed
generation-affects: fill
