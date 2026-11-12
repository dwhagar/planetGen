**Headline:** Warnings that were logged as ordinary lines now go out as warnings, so a warning reader and error.log see them.

### Changed
- **Eight "Warning:" lines moved from normal level to warning (OPS.167).** They were in the sector and block-totals runs, the plan and the run helpers. A lost MySQL connection in the sky unit and a deadlock retry's final failed try are warnings too, and each deadlock retry is a notice. generation-output: unchanged
