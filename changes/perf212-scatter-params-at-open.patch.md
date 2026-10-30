**Headline:** The bright-star and phenomenon scatters now work out their seed, mass limit and clearing when their stage opens, so a run can plan the galaxy and scatter it in one go.

### Changed

- PERF.212: the followers of the bright-star passes and of the phenomenon scatter no longer read the database when they are built. Each stage carries a deferred `scatter-settings` parameter that the wrapper resolves just before the stage opens (the seed, the stellar and compact mass limits, what the luminosity pass clears, the record the last pass stores), and keeps in the run record so the stage's end and a resumed run take the same values even after the pass has emptied what they were read from. The skeleton and "no massive-star pass is stored" refusals are made when the stage opens and `planetgen scatter` raises them as before; the nebula field stage now always follows the phenomena (it has no unit for a galaxy with no outline). Nothing stored changes, no epoch.
