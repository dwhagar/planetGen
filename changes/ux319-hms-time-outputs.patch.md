**Headline:** Every time the pages, progress strips and logs show now reads as hours, minutes and seconds (1h 23m 45s), never as a decimal such as 1.25 minutes.

### Changed
- `format_duration_seconds` (and its browser twin in `unitladder.js`) shows a duration as whole parts: "45s" under a minute, "1h 23m 45s" under a day, "2d 3h 4m" from a day up to a year; sub-second times keep "850 ms" and a year or more keeps the years ladder.
- The Generate page, queue and admin pages, ETAs, the CLI progress ticker, log lines and the slow-query log all go through that one helper; the loading-page time sliders use the browser twin.
- A new test fails if a duration between a second and a year ever carries a decimal or a decimal unit.

generation-output: unchanged
