**Headline:** The uncharted sector rows of the block light totals stage are now made by every worker at once, a run of shells to a unit, instead of one worker making them all one after the other.

### Changed
- **The stage is cut into units (OPS.153).** `expand_refresh_block_lights` checks the parts (the cheap check of OPS.152) and, when they are current, cuts the `block-lights` stage into the roll-up of the totals (one unit) and the uncharted index (one unit for each run of about 300,000 cells with a part, `uncharted_index.shell_runs`, never splitting a shell, at most `MAX_UNITS - 1` runs). The two halves read the parts and write different tables, so they run together. When the parts have to be made again (`rebuild`, or a failed check) the stage stays the one unit of old, since the index must wait for the parts. The index was one process at about 10,000 rows a second, which at the size of Boss's galaxy was hours.
- **The bar and the log follow the units.** `BlockLightsFollower` weighs the units as the expansion costs them (the roll-up, then the runs), logs the cut and one summary line (`Made N block light totals and M uncharted sector rows.`), and says an index unit's lines only when it warns. `make_block_lights` takes `part` (`"all"`, `"roll"`, `"index"`) and `shells`; `fill-block-totals` and a run by hand still make both in one process.
- The rows made are the rows the single unit made (the same cells, uids, names and centers); a unit that is run again makes none twice (`INSERT IGNORE` on the sector's uid).

generation-output: unchanged
generation-affects: none
