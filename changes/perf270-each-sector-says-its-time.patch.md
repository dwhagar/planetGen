**Headline:** A galaxy stage of sectors now says what each sector cost on a line of its own: its address, systems, planets, time and the memory the worker held.

### Changed
- **One timing line for each sector of a stage of 2,000 sectors or fewer (PERF.270).** The starting sector and the neighborhood stage of a New galaxy run had no per-sector line in the run log (each sector's own line was in the debug log), so when Stage 5 took 21,919 s for 109 sectors there was nothing to say which sectors took the time. The line reads `Sector [designation] shell S band B wedge W: 37 systems, 452 planets in 14.2 minutes, the worker at 2,181 MB.`, and the stage's opening line says so. A stage of more than 2,000 sectors keeps the line in the debug log, with the ticker carrying its progress. No generator output changes.

generation-output: unchanged
