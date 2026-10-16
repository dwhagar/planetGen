**Headline:** A cell makes at most 5,000 systems, so the galactic core and bulge fill in about a minute and a half instead of up to 17 minutes.

### Changed
- **The default cell star ceiling is 5,000 systems, from 50,000 (GEN.281).** A system takes about 20 ms to generate, so a cell of 50,000 took 14 to 17 minutes. A galaxy planned earlier at the old default is moved to 5,000 by the migration; `planetgen plan --cell-star-ceiling` still sets it, and it stays one galaxy setting for the admin configuration to edit.
- **A cell the cap still fills past the ceiling raises its cap.** The nuclear cluster's cell expects 15 million systems, and the 100 solar luminosity cap would still make 40,700 of them. A cell now takes the lowest cap, on a ladder of steps of a quarter, whose stars fit the ceiling; the fill, the sector page and the API read the same cap from the address. A cell the cap already fits keeps it.

### Fixed
- **Filling a dense cell no longer takes quadratic time.** Each pre-placed star was checked for overlap against every earlier system of its cell one at a time in Python, which is hours for a cell of 40,000; it is one NumPy pass now (the sector's existing neighbour arrays).

TODO items resolved: none. Added: none.
