**Headline:** The block light totals stage no longer counts every star in the galaxy at the end of a run, so it goes straight to rolling up the totals.

### Changed
- **No recount at the end of a run (OPS.152).** `block_lights.parts_current` and `totals_current` counted the unbuilt stars over `bright_stars` (a walk of every row of a full galaxy) and scanned `phenomenon_scatter` (no index on `built_at`) when the class counters were missing, which kept the stage silent before its first phase line. They now read small tables only: the parts against `phenomenon_scatter_classes`, empty parts against the first waiting star or phenomenon (`EXISTS`, which stops at one), and the totals against the parts' sums. `exact=True` (and `roll(exact=True)`) keeps the old count; `planetgen fill-block-totals` still rebuilds the parts from the objects, so it remains the repair and the proof. The roll also says its first phase ("Checking the parts the scatters kept against their counters...") before the check.

generation-output: unchanged
generation-affects: none
