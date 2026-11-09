**Headline:** A merge train's epoch now alters only the stages its items say they alter, instead of all four.

### Changed
- **The train runner no longer defaults the epoch to every stage.** `scripts/train.py golden` with no `--affects` takes the union of the `generation-affects:` lines of the release notes the train adds (new command `train.py affects` prints it, and `train.py checks` reports it). A note with no line counts as every stage and is named; `generation-affects: none` counts as none; `--affects all` is the old default. Before this, the last 17 epochs (125 to 141) all said plan, scatter, fill and population, so each bring-to-epoch redrew every scatter and the block lights, and a plan stage meant a full regeneration.
- The note check (`bump_version.py`) refuses a `generation-affects:` line that names no stage, an unknown stage, or `none` mixed with stages. The rule is in `changes/README.md`.
