**Headline:** Lighter discs now form their stars later and heavier ones earlier, because the disc's inflow carries the mass term of the Milky Way infall model.

### Changed

- The thin disc's inflow (`galaxy/disk_history.py`) is multiplied by `exp(0.4 t (1 - V / 220 km/s))` per gigayear ([R5103]: equation 21), with the disc speed V = 220 (M / 5e10)^(1/3) km/s. A disc slower than the Milky Way's is fed more and more into the present and a faster one the other way, matching the rising, flat and falling histories of [R5104] (GEN.1246).
- The mass term brings back the downsizing the mass-blind infall of GEN.1230 had flattened: the share of young stars falls about 2.3 times from a 1e10 to a 5e10 solar mass spiral. The reference galaxy forms 2.7 solar masses a year (it was 2.8). Stored star ages, populations and clusters change (generator epoch, `--affects plan,scatter,fill,population`).
