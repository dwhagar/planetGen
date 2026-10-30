**Headline:** Ammonia, sulphur dioxide, hydrogen sulphide, carbon monoxide and hydrogen cyanide now absorb infrared from line-list tables instead of rough band guesses, so airs holding them warm a planet by their real amount.

### Changed

- **Five more infrared gas tables (GEN.852, part 2).** Setup now also downloads the ExoMol correlated-k tables of NH3 (CoYuTe), SO2 (ExoAmes), H2S (AYT2), CO (Li2015) and HCN (Harris), 340 to 470 MB each, reduced like the first three (so the setup download of the infrared data is about 3.4 GB before it is reduced and deleted). The greenhouse model reads them in place of the plateau-and-wing stand-ins, so a world whose air holds any of these five gets a different surface temperature (a trace of ammonia or sulphur dioxide warms more than the guess did). Worlds without these gases do not change. `./update.sh` re-reduces the infrared tables once. Registry rows IR1 and IR2 updated; the CO far-wing entry is gone.
