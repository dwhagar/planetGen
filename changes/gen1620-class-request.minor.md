**Headline:** `planetgen system` can be asked for a planet class letter or a life class, and says what came of the request.

### Changed
- **`--force-class LETTER` and `--force-life-class ID` (GEN.1620).** The first sets the existing `SystemConfig.FORCE_CLASS` (a planet of that class, drawn again within the run's tries); the second sets a new `SystemConfig.FORCE_LIFE_CLASS` (a number such as `6`, or a number and variant such as `4c`), judged on the finished planets and moons the same way. The result then says which body holds the class and how well its conditions fit it (the class match, with "the nearest class the conditions allow" when it is not a full fit), or that none came out, and the run keeps its last try with the warning it already gives. Both are off by default and a default run never sets them, so a seeded system is unchanged (a test compares the two). A bad letter or life class, or `--force-class` with a system file's slots, is refused. The form control for the one-off page follows (UX.302).

generation-output: unchanged
