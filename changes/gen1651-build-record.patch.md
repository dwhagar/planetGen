**Headline:** the twenty private attributes a build stashed on each body are now one `BuildRecord` on the body (`body.build`), so a reader finds them by name and a restored body reads the same defaults a new one does.

### Changed

- GEN.1651: the numbers a pass leaves for a later pass of the same body (the birth spin and eccentricity, the draws and the light an air began with, the capture of an envelope, the carbon sink, the moons a planet was born with and the rest) were set as private attributes (`_birth_spin`, `_air_replay`, `_late_core` and seventeen more) and read back with `getattr(body, "_x", None)`. They are the fields of `planetgen.generation.build_record.BuildRecord`, one per body, each with a default; `Planet`, `life._trial_copy`, `copy` and pickling carry the record as they carried the attributes (a pickle leaves the carbon sink out, a copy keeps it). Nothing is stored differently and nothing draws a random number.
- The generator battery's digest moves because it walks each body's attributes and the record is one more; the body values on fixed seeds are identical (seeds 0 to 13, compared field by field before and after).

generation-output: unchanged
