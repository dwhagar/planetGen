**Headline:** `planetgen sectors --around X Y Z --radius-pc R` generates a sphere of space around any point, nearest sectors first, and a rerun after a stop builds only what is missing.

### Added
- **The region run command (GEN.1128, second slice).** A new mode of `planetgen sectors` and `planetgen galaxy`: `--around X Y Z --radius-pc R` with `--grow-shells K` (whole sector-edge shells added to the sphere) and `--unit-sectors N` (the most sectors in a unit of work, default 27). The sector stage cuts the plan of `galaxy/region_run.py` into one unit of work per region unit, nearest first, and cuts it again each time the stage opens less the sectors already saved, so a stopped run resumes where it was. The same seed gives the same sectors whatever the point, the radius, the unit size or the order regions are run in. `docs/reference/cli.md` describes it. The chart ledger write (GEN.1127) and the top-up stage (GEN.1120) are not part of this slice.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (GEN.1128 stays open for the chart ledger write, which needs GEN.1127).
