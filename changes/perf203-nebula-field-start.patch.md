**Headline:** Drawing the nebula field can now follow another stage in the same run, the part of one run per Generate job that needs no other change.

### Changed

- PERF.203 step 2: the nebula field stage says "the galaxy's skeleton has never been built" and skips a galaxy with no outline when the stage opens (its `prepare` and `expand`), not when the command starts, so an earlier stage of the same run may be what builds the plan. `planetgen scatter --nebula-field` refuses and draws exactly what it did. The other fills (map luminosity, habitability, population, territories) had no start-time check to move. The scatters, link, settle and the sectors planner read the database to build their stages and are PERF.212 to PERF.214.
