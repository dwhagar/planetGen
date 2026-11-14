**Headline:** a world is rated Blue for its air only if it has the 16 kPa of oxygen of a temperate world (class M), not just enough to stay conscious.

### Changed

- The oxygen rung of the Air rating (GEN.1703): Blue needed an inspired minimum of 8 kPa, so a world with 8.5 kPa of oxygen read Blue for its air and could be Blue overall while not being class M. The Blue rung now also needs the oxygen's own partial pressure to be at least the floor of class M's oxygen condition (16 kPa, read from its row by `class_fit.floor_of`, so the two cannot drift apart). Under it the oxygen scores Green at best; the Green, Yellow and Red edges and the too-much-oxygen side are as they were. A species that lowers its inspired minimum lowers the floor in proportion. The reference page of the ratings lists the new row.

generation-output: changed
generation-affects: fill
