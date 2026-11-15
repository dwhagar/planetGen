**Headline:** The stage memory log test reads the figures the way the log now writes them.

### Fixed
- `tests/test_stage_worker_cap.py`: since UX.318 the stage memory line shows gigabytes ("2.37 GB of the 3.05 GB"); the test still looked for "2,370 MB of the 3,050 MB". It reads the new wording. Test only.
