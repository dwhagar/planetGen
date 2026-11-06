**Headline:** A world's Temperature rating now comes from how much of its surface and year a person can be out in without a suit, instead of one average temperature.

### Changed

- The Temperature domain of the habitability score (GEN.1265) is read from the temperature of every place through the year (the series of GEN.1263, the read-offs of GEN.1264): the share of area and time in the Green band, -20 to 45 C, or better, with the wet bulb rule applied to each sample. The share is turned into a score by straight lines between the edges (0.10 of the world livable is the Yellow edge, 0.50 the Green edge, 0.80 the Blue edge) so the geometric mean of the overall score still works. Blue also needs 0.30 of area and time in the Blue band. A sample that needs a suit counts for 0.15 of a livable one, so a world that needs a suit everywhere reads Yellow and one that kills everywhere reads Red. A world that turns once an orbit and is livable only in its terminator belt reads Yellow. The old mean-based score is not a cap: a world with a cold mean and a warm region can now rate better than it did, and an airless world with a mild mean and a scorching day and freezing night rates worse. A sweep of 1976 rocky bodies set the Blue cut (SS15 to SS19 of GEN.725).
- The seven latitude bands (`bandN_mean_k`, `bandN_swing_k`) are read off the same series, and a locked world's ring swing is its real range.
- A body with no spin or orbit stored, and a rogue planet, is scored from one temperature everywhere.
- The temperature series no longer overflows on a hot airless body that turns about as fast as it orbits (it gave NaN).
- Limits printed in the design note: a circular orbit only (`flux_variation` is not in the score), the sea spread evenly, ground ice following the cell mean, a near-resonant slow spinner using the unlocked transport share, and a moon's bands keeping the old climate when the series and the bands disagree on whether it is locked.
- Migration 0205 adds 25 `temp_*` columns to `planets` and `moons`: the livable share and its Blue, Green, cold and hot pieces, the livable area and time, the coldest and hottest place and day, and the best site (its place, mean, coldest, hottest, livable time and longest spell out of the livable band). A galaxy made before gets them empty until it is regenerated; the stage re-run fills them.

generation-output: changed
