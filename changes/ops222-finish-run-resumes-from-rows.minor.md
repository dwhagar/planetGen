**Headline:** Finish the last run carries a stopped run on from what its records say is done: finished units are not done again and a unit that stopped starts after the pieces it finished.

### Changed
- `finish-run` plans the stopped run again and, for the stage that stopped and each stage after it, the new run subtracts what the stopped run's rows say is done (`queue/resumed.py`): a unit the rows say is done is marked done and its stored result kept for the stage's conclusion, and a unit that was not done gets the cursor and the stored piece results of the unit it continues, so it starts after the pieces it finished. A stage whose units were all done is only concluded.
- A unit of a command that cannot resume is queued with `again` when its stage carries on from a stopped run (the sector stages skip the sectors already saved); its row is cut as before.
- `finish-phenomena` resumes after its last finished piece: the drawn special rows are the first piece's stored result, the write and the record commit through the new `ctx.commit_piece()` (data, cursor and result in one fenced transaction), and the indexes are built only where missing.
- The Generate page and the Queue page say that a stopped stage carries on from what the run's records say is done, not from its start.

### Removed
- `finish_run.REWIND_WITH` and `finish_run.AGAIN`: no stage is rerun from its start, and no stage is rewound with the one before it, to make a stopped stage safe to repeat.

generation-output: unchanged
generation-affects: none
