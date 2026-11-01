# Volatile inventory against the totals of the Earth, Venus, Mars and Titan (GEN.813)

Status: research note and calibration, written on the night of 2026-10-24 (as the lane reported) by the Citations lane. It answers finding 2 of `reports/air-fixed-point-findings-2026-10-24.md`. Two brackets of `mixture_rocky.INVENTORY` change (the dry-class carbon and the nitrogen per carbon). Nothing in the generator draws the mixture yet, so no stored value changes and there is no epoch; the check on the fixed point itself waits on GEN.411 and GEN.422 (section 6).

## 1. The answer in short

- The finding was right and its cause is one number. The delivery weights put a planet at the Earth's distance in the dry class (b1) 85 percent of the time (x = a / a_snow = 0.37 gives s = 0.076, weight (1 - s)^2), and b1 carried 0.6 weight percent CO2: 1640 ppm of carbon, sixteen times the 100 ppm of the bulk silicate Earth (BSE). Nitrogen follows carbon by N/C, so 37 percent of Earth-distance draws held more than ten times the BSE's nitrogen (one seed: 87 bar of N2 on a 1.8 Earth-mass planet).
- The 0.6 is Ortenzi's DELIVERED dry value. The same table's STORED dry value, the carbon that stays in the mantle, is 0.018 weight percent (49 ppm C). The inventory is read as the carbon in the mantle, so the delivered value was the wrong column for it. The dry-class CO2 is now a log bracket from the stored 0.018 to the wet High case's delivered 0.1 weight percent. Its log-middle, 0.042 weight percent, is 115 ppm of carbon against the BSE's 100.
- The nitrogen per carbon bracket was 1e-3 to 0.1. The measured ratios are the BSE 0.020, the CI chondrites 0.059, the air of Venus 0.074 and the air of Mars 0.055 (Bergin et al. 2015 Table S1). The bracket is now 0.005 to 0.1, whose log-middle (0.022) is the Earth's.
- With these, an Earth, Venus or Mars twin draws a median of 80 to 110 ppm of carbon (BSE 100), the BSE's carbon and nitrogen sit inside the 5 to 95 percent range of the draw for all three, the air of Venus is what a degassed share of 0.43 of the draw brings up, Mars's source magma range (5 to 1000 ppm CO2) holds the median, and no Earth-distance draw holds ten times the BSE's nitrogen (it was 37 percent).
- Titan is not in the rocky inventory. Its air is set by the `titan` family from the pressure it is given, and the nitrogen of an icy moon needs an ammonia-ice inventory that does not exist (section 5); the check there is the pressure.

## 2. Sources, and what was read

| Source | Read | Used for |
|---|---|---|
| Bergin, Blake, Ciesla, Hirschmann and Li 2015, PNAS 112, 8965, "Tracing the ingredients for a habitable Earth from interstellar space through planet formation" (arXiv:1507.04756) with its supporting information | full text, Table S1 and the Earth section | BSE C/Si 1.1e-3 (plus or minus 2.0e-4) and N/Si 2.5e-5 (plus or minus 5e-6), C/N 49 plus or minus 9; Earth surface carbon 9.3e22 g and nitrogen 6.4e21 g; mantle carbon 20 and 165 ppm for the depleted and enriched sources; the ratios for the CI chondrites (17), the air of Venus (13.5) and the air of Mars (18) |
| NASA planetary fact sheets for Venus and Mars (nssdc.gsfc.nasa.gov) | read | Venus 92 bar, 96.5 percent CO2, 3.5 percent N2, g 8.87 m/s2, radius 6051.8 km; Mars 6.36 mbar, 95.1 percent CO2, 2.59 percent N2 |
| arXiv:2302.04241, "Constraints on the Size and Composition of the Ancient Martian Atmosphere from Coupled CO2-N2-Ar Isotopic Evolution Models" | full text | the CO2 of the Martian source magma, 5 ppm at IW-1 to 1000 ppm at IW+1 (it cites Hirschmann and Withers 2008 and Grott et al. 2011, which were not opened: second-hand); 1.9 ppm of nitrogen (same); the atmosphere of 0.3 to 1.5 bar CO2 and 0.1 to 0.5 bar N2 at 3.8 Ga; 300 to 1400 mbar of CO2 in carbonates |
| Erkaev, Scherf, Kislyakova, Lammer et al. 2020, MNRAS, "Escape and evolution of Titan's N2 atmosphere constrained by 14N/15N isotope ratios" (arXiv:2010.12342) | the parts on the ice and nitrogen budget | Titan's mean surface pressure 1.467 bar, the origin of the nitrogen in ammonia ice, and the CI-like nitrogen content of the rock |
| Ortenzi et al. 2020 Table 1 (already read, part 05 section 3.1) | already read | delivered and stored CO2 of the Low, High and Dry cases |

The Earth's silicon fraction, 21 weight percent of the BSE (45 weight percent SiO2), turns the Si ratios into ppm; it is recalled, not read (tier 3).

## 3. The totals

All as ppm by weight of the mantle (the mass the inventory is per).

| Body | Carbon | Nitrogen | Basis |
|---|---|---|---|
| Earth (BSE, mantle plus surface) | 99 (82 to 117) | 2.6 (2.1 to 3.1) | Bergin Table S1 ratios times 21 weight percent Si |
| Earth surface only | 23 | 1.6 | 9.3e22 g C and 6.4e21 g N over 4.0e27 g |
| Venus air | 39 | 3.3 | 92 bar is 4.77e20 kg of air; 3.3e24 kg of mantle; a floor, the interior is unknown |
| Mars source magma | 1.4 to 273 | 1.9 | 5 to 1000 ppm CO2, redox dependent; nitrogen the silicate Earth's, adopted |
| Mars air | 0.013 | 0.001 | 6.36 mbar over 4.8e23 kg of mantle (a floor, mostly lost) |

Titan: 1.467 bar is 9.0e18 kg of air, almost all N2, so about 8.9e18 kg of nitrogen, 65 ppm of the whole 1.35e23 kg body, and its nitrogen per carbon in the air is above 30 (5.65 percent methane in N2). The CI chondrites carry about 2000 ppm of nitrogen (N/Si 0.042, silicon 10.6 weight percent recalled, tier 3), so the nitrogen of Titan is a small share of an ice-rock inventory but 25 times the BSE's: the two cannot share one bracket.

## 4. The change and its effect

| Draw (median, 5 to 95 percent of Earth-distance draws, ppm) | Before | After | Earth BSE |
|---|---|---|---|
| Carbon | 1637 (27 to 1637) | 101 (27 to 250) | 99 |
| Nitrogen | 12.9 (0.2 to 145) | 1.1 (0.1 to 12) | 2.6 |
| Share of draws with more than 26 ppm N (ten times the BSE) | 0.37 | 0.00 | |

The remaining spread is wide on purpose: the 27 ppm low end is Ortenzi's Low case (6 times under the BSE) and the 250 ppm top is the High case. Venus and Mars draws are the same inventory (the delivery weights differ: Venus 93 percent dry, Mars 65 percent dry, 31 percent Low).

Registry (GEN.725): AM13 keeps the wet classes (tier 1, tweak level 0); AM14 is the nitrogen per carbon bracket (tier 1, tweak level 2: read, then tuned so that its log-middle is the Earth's); AM50 is the dry-class carbon (tier 1, tweak level 2: both ends read from Ortenzi's table, the choice of columns tuned to the BSE carbon). Tests: `src/tests/test_volatile_inventory_totals.py`.

## 5. Limits and what is not done

- **Mars** is checked only at the source magma, which spans 200 times; the 50 percent coverage of the draw is not a test of Mars. The outgassed and lost CO2 (the 70 mbar to 13 bar of Hirschmann and Withers through the Mars paper) depends on the degassed share of Mars, which the interior heat model supplies in the fixed point.
- **Venus** is a one-sided check (the air is a floor); the mantle carbon of Venus is unknown. The degassed share 0.43 is the one the supply law note takes from the same air.
- **Titan** has no inventory: the `titan` family takes its pressure from the retention and supply terms and the nitrogen of an icy moon is ammonia ice (Erkaev et al.), which the rocky inventory does not carry. An ice-moon nitrogen inventory (ammonia over water, with the outgassing efficiency of Erkaev et al. Table 2) is a new item; the Titan check will be the 1.467 bar.
- **Sulfur** per carbon (0.1 to 3) is untouched; the BSE ratio of 1.7 and the CI 0.6 are inside it. The sulfur sink is finding 4 (a separate item).
- **Correlation.** Carbon and nitrogen are drawn independently inside their brackets. The BSE ratio C/N = 49 plus or minus 9 would allow a tighter pairing; the independent draw keeps the spread of the measured ratios (0.013 to 0.074) rather than the Earth's alone.
- **The water-rich class** (b4) is GEN.814's.

## 6. What waits on the fixed point

GEN.411 and GEN.422 (the air fixed point) are not in the repository yet. The "Done" line of GEN.813 asks that the four bodies reproduce within the stated spread. At the level of the inventory this is done and tested here. The full check (an Earth twin at about 1 bar of N2, a Venus twin at about 3 bar, a Mars twin thin and a Titan at 1.467 bar) is a run of the fixed point; it is the same run GEN.819 needs, and it will be done there.

- TEST.547 (2026-10-31, as the lane reported): GEN.1005 replaced the bracket tuned to the bulk silicate Earth with the one derived from Blanchard et al. 2022 (0.029 to 0.108 weight percent CO2, 80 to 294 ppm C of the body, log middle 153 ppm C). The test of the middle now holds it to the paper's own mantle carbon, 140 +/- 40 ppm (100 to 180), and still asks that the Earth's 99 ppm lie inside the bracket; the old bound of 20 percent about 99 was the tuning, not a measurement.
