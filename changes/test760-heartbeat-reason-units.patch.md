**Headline:** three tests that expected a duration written one way now read it from the formatter the pages use.

### Fixed

- TEST.760: the heartbeat rules test looked for "600" in a quiet reason, the Queue page row test for "2.0 s each" and the update guard test for "30 seconds", while `format_duration_seconds` writes "10 minutes", "2 s" and "30 s" (the limits and clocks agree; only the wording moved). Each test now builds the expected text from the same formatter, so the next change of wording moves them together.
