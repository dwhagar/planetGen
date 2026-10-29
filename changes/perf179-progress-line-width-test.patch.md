**Headline:** The progress bars in a job log are guarded by a test that they are never cut off at the log's width.

### Fixed
- **PERF.179: progress bars cut off in the job log ("elapsed 0:0... remaining 0:00:...").** The cut was rich assuming 80 columns for output that is not a terminal; the bars are drawn 160 wide there (already in place). New tests render the bars to a log the way a job does and fail if the elapsed or remaining time is cut, if the line passes the log width for the longest realistic step name, count and run time, or if a terminal or a `COLUMNS` setting stops being left to rich. No stored value changes.
