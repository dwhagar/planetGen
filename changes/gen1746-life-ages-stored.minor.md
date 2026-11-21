**Headline:** The life ages of a world, the stretches of its life between mass extinctions, are now kept as rows in new tables `planet_life_ages` and `moon_life_ages` (nothing reads them back yet)

### Added

- GEN.1746 slice 3: the life pass (`life._history_of`) cuts the history of every world with life into its life ages (`life_ages.ages_of`, an age ends at an extinction that kills 30 percent or more of the species, the Boss's pick) once, and the body keeps them as `life_ages` (a list of `life_ages.entries_of` dicts: number, start, end, cause, kill share; `None` for a world with no life). The air history record takes the same ages instead of cutting them again. A world with life but an air under 500 Pa, which gets no air history record, still gets its ages.
- Two new tables, `planet_life_ages` and `moon_life_ages` (migration, schema, models, column units, fingerprint and command touches), one row for each age, written with the body, loaded back with it (`store.body_life_ages`) and kept by an edit. Indexed on the age number and its end, so "the worlds in their third life age" is a short read. A rogue planet has no star age and so no ages. Nothing reads the rows back yet; the API comes after.

generation-output: changed
generation-affects: fill
