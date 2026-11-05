**Headline:** The sector resume and parallel-run tests now finish in two to four minutes each instead of fifteen to twenty-three, and no longer fail on the sector rows a rogue planet creates; no generated output changes.

### Changed
- **Whole-sector tests made small (TEST.417).** `test_gen_resume.py` (eight cases, 450 to 920 seconds each) and `test_work_queue.py` (1,354 seconds) filled 12 to 54 sectors three times over, and a sector costs about 3 seconds once the first fill's draw down to the backfill's deepest floor is left out (it is half to nine tenths of a sector's cost and the same in every run, so it says nothing about resuming). The resume tests now turn that draw off (`tests/first_fill.py`, as the bright-star scatter tests do) and use shells 0 to 2 (12 to 14 sectors a mode, not 12 to 54); the parallel-run test fills 4 sectors of shell 1 with 2 workers, not 50 of shell 2. Every check is kept: the same sectors after a resume and an uninterrupted run, each filled once, no orphan system or link, nearest-system lists equal to a recompute, every run record done. Eight cases of `test_gen_resume.py` and all of `test_work_queue.py` also failed on main: a rogue planet that lands in a cell nobody charted gets a sector row there (GEN.1111 part 7), so the sector counts were wrong. The snapshot now counts only charted sectors (those in the chart ledger) and checks that a resumed run and an uninterrupted run end with the same uncharted rows.

generation-output: unchanged

TODO items resolved: none (TEST.417 continues).
