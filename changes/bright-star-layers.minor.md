### Added

- Bright stars can now be scattered in stages (PERF.5). `generate.py plan
  --bright-stars-down-to 100` keeps every bright star already placed and
  adds only those from 100 up to the galaxy's current star-fill level (500
  by default), then lowers the level. The Generate page shows the level
  and has a new "Add a dimmer layer of bright stars" panel. Sectors already
  filled are left out, since their own systems already include stars that
  bright. No database change.
