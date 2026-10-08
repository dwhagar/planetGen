### Fixed
- The System Map marks itself ready (`data-ready`) once its click handlers are on, and its browser test waits for that instead of a fixed time, so it no longer fails now and then under load (TEST.104).
