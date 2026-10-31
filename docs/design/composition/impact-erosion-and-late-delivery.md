# Impact erosion and late delivery (GEN.739)

Status: research note and module, written on the night of 2026-10-24 (as the lane reported) by the Citations lane. The module is `src/planetgen/physics/impact_history.py` and its supply-law hook is `surface_pressure.impact_inventory`. Nothing in the generator calls them yet: the surface-pressure supply law is still being wired (GEN.411, GEN.422), so no stored value changes and there is no epoch. The numbers in section 4 are COMPUTED here from published formulas; they are a model result, not a measurement.

## 1. The answer in short

- After the giant impacts a rocky planet goes on accreting leftover planetesimals. The mass is about 0.7 weight percent of an Earth mass for the Earth, and goes as (M / M_earth)^(-1/3) (cross section over mass): Mercury 1.9, Venus 0.75, Earth 0.7, Mars 1.5 percent, against 2.7, 0.8, 0.7 and 1.4 in the N-body runs of Brasser et al. 2016 (Table 1; scatter 0.5 to 4.5 percent). The Moon is off the law (0.025 weight percent from the highly siderophile elements, as quoted in the sources above); nothing here is fitted to it.
- The same impacts do two opposite things to the air. They deliver volatiles (water, CO2, N2) and they blow the air away. The balance for each species is the one of Sakuraba, Kurokawa and Genda 2019: dm/dS = (1 - zeta) x - eta m / m_air, with S the cumulative impactor mass, x the species' mass fraction in the impactors, zeta the share of the impactor's own vapour that escapes and eta the air mass one unit of impactor mass removes.
- With the carbonaceous reference composition of Sakuraba et al. (10 percent water, 2 percent CO2, 0.2 percent N2) scaled by a drawn factor of 0.01 to 1 (0.1 is the enstatite-chondrite-like veneer that the Earth's isotopes ask for), the Earth twin ends with 5.9e17 to 7.3e18 kg of N2 in the air; the real 3.9e18 kg lies inside. The water the impactors bring and the planet keeps (above what the air can hold) is 3.7e20 to 3.7e21 kg; the ocean's 1.4e21 kg lies inside.
- Venus, with no ocean and no carbonate sink, holds everything the impacts do not blow away: 1.4e17 to 9.3e18 kg of N2 against the observed 1.1e19 kg.
- Mars comes out with much more N2 than it has (1.8e18 kg at scale 1, 8.8e16 kg at 0.1, against 4.4e14 kg observed). This is the result of Sakuraba et al. 2019 as well: late impacts alone do not make Mars thin. The thin air of Mars is the escape of the retention rows of the pressure supply law, not erosion by impacts. A planet the size of Mars needs those rows (GEN.411/422), and the module's test says so.
- About a tenth to a fifth of the impactor vapour escapes (zeta 0.10 for the Earth, 0.13 for Venus, 0.17 for Mars).

## 2. Sources, and what was read

| Source | Read | Used for |
|---|---|---|
| Sakuraba, Kurokawa and Genda 2019, Icarus 317, 48, "Numerous chondritic impactors and oceanic composition" (arXiv:1805.07094) | full text | the balance Eq. 1 and 2; Svetsov's erosion as Eq. 5 to 7; Shuvalov's vapour escape Eq. 8; the carbonaceous reference composition (section 2.4); the size distribution (section 2.6); the surface limits of an ocean and carbonate cycle (section 2.3); the Earth, Venus and Mars check |
| Brasser, Mojzsis, Werner, Matsumura and Ida 2016, Earth Planet. Sci. Lett. 455, 85 (arXiv:1609.01785) | full text | Table 1: late-accreted mass, impact speeds and escape speeds for Mercury, Venus, Earth and Mars; the 0.5 to 0.8 weight percent of Walker 2009 |
| Schlichting, Sari and Yalinewich 2015, Icarus 247, 81 (arXiv:1406.6435) | full text | the stronger erosion by the largest impactors in the cap regime: used only as the upper bracket on the erosion (section 5) |
| Raymond, Schlichting, Hersant and Selsis 2013, Icarus 226, 671 (arXiv:1306.4325) | full text | the stochastic nature of the late veneer (the planet-to-planet scatter that the log-normal width stands for) |

Shuvalov 2009 were not opened; their formulas are used as printed in Sakuraba et al. 2019 (second-hand). The registry rows AM43 and AM44 say so.

## 3. The model

- **How much accretes** (`late_accreted_fraction`): 0.7 weight percent times (M / M_earth)^(-1/3), log-normal with sigma 0.65 in ln (the scatter of the Earth value, 0.7 plus or minus 0.5), capped at 3 percent, drawn once per planet and recorded as `late_accreted_percent`.
- **What it brings** (`veneer_volatiles`): H2O 0.10, CO2 0.02, N2 0.002 of the impactor mass (carbonaceous reference of Sakuraba et al. section 2.4), times a scale drawn log-uniform between 0.01 and 1, recorded as `veneer_volatile_scale`. The scale is the one free number of the model; the Earth's N2 and the ocean bracket it (section 4).
- **How hard it hits** (`impact_speed_ms`): the speed at the surface is sqrt(v_esc^2 + v_inf^2) with v_inf = 0.6 of the orbital speed. Brasser et al. give 24.3, 21.1 and 13.4 km/s for Venus, Earth and Mars; the module gives 23.4, 21.1 and 15.3 (within 16 percent, tested).
- **What it removes** (`erosion_per_impactor_mass`): Svetsov's air mass removed per unit impactor mass for an impactor of diameter D, with the plume fraction that outruns the escape speed, averaged over a size distribution whose mass is flat in the diameter from 10^-1.5 to 10^3 km (Bottke et al. 2005 as used by Sakuraba et al.). A thin air loses less per impactor mass, a big impactor less than a small one, and an impact too slow to launch the plume past escape erodes nothing (all tested).
- **The caps** (`caps_pa`): the gas a species can hold is limited by the surface. An ocean holds the water above its vapour pressure, a carbonate cycle the CO2 above its balance, an ice cap the rest. The excess is kept in the ocean, rock or ice, not eroded. The caller gives the caps (the module knows no climate), and the test uses the surface limits of Sakuraba et al. section 2.3.
- **To the supply law** (`surface_pressure.impact_inventory(balance, mantle_kg)`): the air plus the kept mass, as moles of H, C and N atoms per 100 g of mantle (H = 2 H2O, C = CO2, N = 2 N2).

## 4. Checks against the terrestrial planets

Late-accreted fraction at the mean (0.70, 0.75, 1.47 percent for Earth, Venus, Mars). Masses in kg, scale 0.1 and 1.

| Planet | Speed km/s | zeta | N2 delivered | N2 in the air | Observed N2 in the air | H2O kept | CO2 kept |
|---|---|---|---|---|---|---|---|
| Earth | 21.1 | 0.10 | 8.4e18 to 8.4e19 | 5.9e17 to 7.3e18 | 3.9e18 | 3.7e20 to 3.7e21 (ocean 1.4e21) | 7.3e18 to 6.9e20 |
| Venus | 23.4 | 0.13 | 7.3e18 to 7.3e19 | 1.4e17 to 9.3e18 | 1.1e19 | none (no ocean) | none (all in the air) |
| Mars | 15.3 | 0.17 | 1.9e18 to 1.9e19 | 8.8e16 to 1.8e18 | 4.4e14 | 7.8e19 to 7.8e20 | 0 to 4.9e19 |

Tests (`src/tests/test_impact_history.py`): the Brasser means inside their scatter, the Earth N2 and ocean inside the bracket of the scale, Venus within a factor of five of its N2, Mars above 100 times its N2 (so that escape, not erosion, has to explain the thin air), zeta between 0.05 and 0.25, mass conservation and the caps held.

## 5. Limits and what is not done

- **Erosion strength.** Svetsov's formula (used by Sakuraba et al.) is the lower bracket. Schlichting et al. 2015 find a stronger erosion by the largest impactors (about 10 times that of Svetsov's law at 15 km/s, as recorded in registry row AM43). The upper bracket would lower the Earth's and Venus's N2 by the same order; with the veneer scale free over two dex, the Earth's N2 alone cannot decide between the two erosion laws. A later hydrocode check (Svetsov 2007, new) would.
- **Single-species caps.** The caps are supplied by the caller one species at a time. Water and CO2 compete for the same ocean and the same climate; the real coupling is part of the pressure supply law (GEN.411/422).
- **Solar-System remnant mass.** The leftover planetesimal mass of 1e-3 Earth masses is the Solar System's. For other systems it is a stand-in (registry row AM36).
- **Not wired.** `impact_inventory` is the hook for the supply law; the law does not call it until 4c lands. A follow-up (GEN.739b, to be filed by TODO) is the wiring and the epoch that goes with it.
- **The Moon** is off the mass law by a factor of about 30 and is not covered.

## 6. Registry

Rows AM36 to AM44 of `docs/design/research/gap-filling-plan.md` (and the generated `stand-in-values.md`) list each constant with its source tier and tweak level, and the test that would falsify it.
