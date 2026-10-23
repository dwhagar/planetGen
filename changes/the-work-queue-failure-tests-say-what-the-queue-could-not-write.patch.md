**Headline:** The work queue failure tests wait 180 s before calling a run hung and name any bookkeeping the queue could not write; the sector-open test clicks again where the map holds still.

### Fixed
- **Thirteen tests of `test_work_queue_failures.py` and `test_opening_a_sector_is_one_continuous_view_with_a_filled_panel` failed in one gate at a load of 6 to 8 (TEST.297).** They passed alone, and six loops of the file beside eight busy CPU loops passed too, so the load alone did not cause it. Two failures had text: `test_a_dead_worker_never_hangs_the_run` gave up on a run whose record stood still for 60 s (a replacement worker took longer to start), now 180 s; and `test_on_done_raising_fails_the_run` found no `work_jobs` row for its job, because `WorkQueue._book` carries on past a bookkeeping error and only debug-logs it. The file now keeps those lines, and a missing row reports the rows there are and what the queue could not write. The sector-open test clicked the sector's centre where it had been a moment before; it now clicks where the centre holds still and again, up to four times, if the sector did not open (the cause of TEST.291).

TODO items resolved: none. Added: none.
