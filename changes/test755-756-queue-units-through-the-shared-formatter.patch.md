**Headline:** The queue's quiet-stage line and the splitter's reason write their times and memory through the shared formatter.

### Fixed
- The "no unit has finished for ..." line and the splitter's reason line spell durations and memory with `format_duration_seconds` and `format_megabytes` (UX.318), so the hand-formatted-units test passes again. The run-record test and the stall-words test expect the same words (1.49 GB, 35 minutes, 1.17 hours).
