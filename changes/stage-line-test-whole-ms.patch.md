**Headline:** The stage-line test expects whole milliseconds, as stage lines have read since TEST.814.

### Fixed
- `test_stages.py`: `test_a_stage_line_gives_the_time_the_stage_took_on_the_units_ladder` still expected "2.86 ms"; since TEST.814 (PR 7355) a duration under a second reads in whole milliseconds (`logrecord.duration_text(0.00286)` is "3 ms"), so the test was red on main. It now expects "3 ms".

generation-output: unchanged
