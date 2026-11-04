**Headline:** The galaxy's thin disc now forms its stars from the Milky Way's inside-out infall, so the reference galaxy forms 2.8 solar masses a year, not 4.0.

### Changed

- The thin disc's gas inflow (`galaxy/disk_history.py`) is an exponential infall with a time scale that grows outward (the linear fit of [R5095], limited to 2 and 8 Gy by [R5094]), in place of the growth of the stars a halo holds. The lighter median halo of the halo mass function had pushed that growth late, so the reference galaxy formed 4.0 solar masses a year against the 1.65 measured (GEN.1242) and the disc's present star formation was 2.7 times the main sequence's (GEN.1230). The gas at 8 kpc is 10.3 solar masses per pc2 and the depletion time 3.4 Gy.
- The central molecular zone inflow time is 0.6 Gy (it was 1.0), so the Milky Way shape keeps its measured 0.8 solar masses a year; the magnetar young-star share is 0.18 percent (it was 0.16).
- Removed with it: the halo growth track, the stellar-to-halo evolution and the cosmic clock in `disk_history`, and their tuning constants. Stored star ages, populations and clusters change (generator epoch, `--affects plan,scatter,fill,population`).
