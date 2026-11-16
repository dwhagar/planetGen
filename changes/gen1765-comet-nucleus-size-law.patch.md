**Headline:** A star's own comets now follow the measured size law of comet nuclei: many small ones, few large.

### Changed
- GEN.1765: a star-bound comet's nucleus diameter was drawn flat between 0.5 and 20 km, which made the median nucleus 10 km. It is now drawn from the cumulative power law of the Jupiter family (the number above a diameter falls as D to the minus 1.92, [R4139, R5296]) between 2.5 km (the edge of the fit) and 20 km, so the median is about 3.6 km and about one comet in six is above 6 km. The mass, density and strength of a comet follow its diameter, so they change with it. New constant `COMET_NUCLEUS_CUMULATIVE_INDEX`; stand-in rows CMT6 (the law, read) and CMT9 (the hand-set upper end, until the belt's largest body sets it). The supply of comets from the system's outer belt is the next slice.

generation-output: changed
generation-affects: fill
