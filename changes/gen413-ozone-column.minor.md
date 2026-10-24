**Headline:** A world's ozone layer now follows its oxygen and its star, and the UV index reads it.

### Changed

- **Ozone column (GEN.413).** Every rocky planet and moon draws an ozone column in Dobson units from its oxygen and its star's temperature (`physics/ozone.py`, stored as `planets.ozone_column_du` and `moons.ozone_column_du`, migration 0135, schema v135): the fifth-degree fit of Kozakis et al. 2022 for five hosts (G0V to M5V, valid for 0.01 to 150 percent of the present oxygen), mixed in temperature between them, scaled to the F-star columns of Rugheimer et al. 2015 above 6000 K, times a structure factor of 0.85 to 1.71 that spans the published models (Earth air gives 163 to 328 DU). Hot hosts peak at a quarter of the present oxygen and cool ones rise to the highest. Below 0.01 percent of the present oxygen the column is 0.01 to 0.1 percent of the host's full column (15 percent for an F host); a world with no oxygen has none.
- **The UV index and the ozone-loss flag read the column.** The square-root proxy on the oxygen pressure is gone; the index is the star's ultraviolet against the Sun's times the column over 300 DU to the power -1.6, as before, and the flag needs a column above zero. A row made before this revision has no column (NULL).
- **Not yet.** The ozone is not a gas of the mixture (it joins when the mixture table lands, GEN.411), the Chappuis absorption table waits for the sky-colour items, and the galactic events do not yet thin the column (the flag stays).

generation-output: changed
