# Atmospheric escape: a thermosphere model, Jeans and diffusion-limited loss, and the XUV history (GEN.738)

Status: model and module, written on the night of 2026-10-24 (as the lane reported) by the Citations lane. The module is `src/planetgen/physics/thermosphere.py`; `mixture_escape.apply` (GEN.693) now calls it in place of the logistic filter. Nothing in the generator draws the mixture yet (it waits on GEN.411 and GEN.422), so no stored value changes and there is no epoch. The numbers in section 4 are COMPUTED here; they are model results, not measurements, except where a measured value is named.

## 1. The answer in short

- What a body loses from its air is fixed by three things: how hot the top of the air is (the exobase temperature), how thin it is there, and how much XUV the star has put out over the body's life. The old filter had none of them: a fixed ratio of the exobase to the equilibrium temperature (1.3 or 4, chosen by a threshold at half CO2), a logistic in the Jeans parameter and a rule of thumb for the age. Those three (rows AM19 to AM21) are gone.
- The exobase temperature is the mesopause temperature (0.78 of the equilibrium temperature) plus a rise that goes as the XUV flux to the power 0.96 and falls as the CO2 of the air grows, `1 / (1 + (x_CO2 / 0.0126)^0.70)`. The two shapes and four constants are a fit to the 1-D upper-atmosphere model of Johnstone, Guedel, Lammer and Kislyakova 2018: the exobase temperatures of the modern Earth over a hundred-fold range of CO2 (1310 K at 4e-6 down to 290 K at 0.4) and the 4,600 K of the young Sun's 4.9 times XUV. A CO2 air is cold (Venus 274 K, Mars 186 K), an N2 air is hot (the Earth 1,204 K).
- From that temperature the module finds the exobase (the height where the mean free path equals the scale height) on a temperature profile that rises from the mesopause, the Jeans flux of every species through it, and caps each species lighter than the air by the diffusion-limited flux of Hunten (1973): a light gas cannot leave faster than the air below can carry it up.
- The history is followed. The star's XUV falls as the star ages (`activity.coronal_xuv_at`); the flux, the temperature and the flux of every species are recomputed on 48 log-spaced steps from 30 Myr to the body's age, and each species' column is lowered by what leaves. The share of a species that stays is its column at the end over the one at the start. The air that is left is what the next step sees, so a hydrogen-rich air carries the heavier gases only while it still has its hydrogen.
- The checks: the Earth keeps its N2, O2, Ar and CO2 and loses H2 and He; Venus keeps all of its air (its exobase is cold); Mars keeps CO2 and N2 and loses H2; Titan's diffusion-limited hydrogen flux (1.0e28 per second) is within 15 percent of the 9.2e27 that Cassini measured, and its Jeans flux is a factor of 4.5 under it (section 4).
- The exobase search is bounded (GEN.841): at most `MAX_EXOBASE_STEPS` (200,000) steps, where the longest search over the whole range of bodies swept (radius 1e6 to 7e7 m, gravity 0.5 to 300 m/s2, mesopause 5 to 400 K, exobase to 5000 K) takes about 15,000, so the bound never acts on a real air; a NaN, infinite or non-positive input is an error and no longer a search that never ends. The binary diffusion parameters and molar masses are cached, which changes no value (the shares of 120 random airs are bit-identical before and after).
- The mass loss of a hydrogen envelope (a sub-Neptune) is the analytic fit of Kubyshkina et al. 2018 to 7,000 hydrodynamic models (`hydro_based_mass_loss_g_s`); it reproduces the nine planets of its Table 2 to a factor of 2.3 or better (seven of the nine within 5 percent, the other two within the 2.5 of the test). It is not wired into the envelope rows yet (section 5).

## 2. Sources, and what was read

| Source | Read | Used for |
|---|---|---|
| Johnstone, Guedel, Lammer and Kislyakova 2018, A&A 617, A107, "Upper atmospheres of terrestrial planets: Carbon dioxide cooling and the Earth's thermospheric evolution" (arXiv:1806.06897) | full text and Figures 11, 12 and 13 | the Earth exobase temperature against the CO2 mixing ratio (Fig. 12: 1310, 1260, 1155, 985, 540 and 290 K, read from the figure), the 4,600 K at 4.9 times the XUV (their model 3 Gyr ago), the exobase altitudes (475 km now, 2,000 to 2,500 km at the hottest), the mesopause of 200 K at 1e-2 mbar, the statement that the Earth model is never hydrodynamic and that its Jeans rates of N and O are orders below that of H |
| Kubyshkina, Fossati, Erkaev et al. 2018, Frontiers in Astronomy and Space Sciences 5, 19, "Overcoming the limitations of the energy-limited approximation for planet atmospheric escape" (arXiv:1810.06920) | Eq. 9 and 10, Table 1, Table 2 | the hydro-based mass loss of a hydrogen envelope and its nine test planets |
| Zhang 2024, "Atmospheres of Solar System Moons and Pluto" (arXiv:2410.04595) | the Titan section | Titan's hydrogen escape (limited by the maximum diffusion rate, about 9.2e27 per second), the exobase at 1,500 to 1,600 km |
| Erkaev et al. 2020 (arXiv:2010.12342), already read for GEN.813 | the Titan parts | Titan's surface pressure 1.467 bar |
| Hunten 1973 (J. Atmos. Sci. 30, 1481); Chapman and Cowling's hard-sphere diffusion; Watson, Donahue and Walker 1981; Erkaev et al. 2007; Chamberlain 1963; Breck 1974; Bates 1959; Marrero and Mason 1972 | not opened | the diffusion-limited flux, the diffusion coefficient, the energy-limited rate, the exobase condition and cross section, the kinetic diameters, the temperature profile and the H2 in N2 diffusion coefficient of the test, as recalled (source tier 3; the registry rows say so) |

## 3. The model

- **Exobase temperature** (`exobase_temperature_k`): `T = 0.78 T_eq + 1095 K (F / F_E)^0.96 / (1 + (x_CO2 / 0.0126)^0.70)`, `F_E = 4.6e-3 W/m2` the Earth's XUV (`constants.EARTH_XUV_FLUX_W_M2`). `F` is held to 4.9 Earth units, the top of the range the fit was checked in. The mesopause of a body whose lower air is not in radiative balance with the star (Titan's haze: 150 K) can be given.
- **Exobase** (`exobase`): the air rises from 1 Pa at the mesopause temperature to the exobase temperature on `T(z) = T_exo - (T_exo - T_meso) exp(-z / L)` with `L` 2 scale heights of the mesopause (12 km for the Earth), in hydrostatic balance with gravity falling off as the inverse square. The exobase is where `n sigma H = 1`, `sigma = 2e-19 m2`. This profile gives 441 km at 1,155 K and 2,160 km at 4,600 K against the 475 and the 2,000 to 2,500 km of Johnstone et al. (an isothermal air from the mesopause, the first try, put the 4,600 K exobase at 14,800 km and made every Earth-mass draw with no CO2 airless; it is not used). The altitude is above the 1 Pa level, not above the surface or the 1 bar level.
- **Jeans flux** (`jeans_flux`): `n v0 (1 + lambda) exp(-lambda) / (2 sqrt(pi))` with the density of the species `n = f n_exo`, `v0` its most probable speed and `lambda = G M m / (k T r)` at the exobase.
- **Diffusion limit** (`diffusion_limited_flux`): `b f g (m_air - m) / (k T_meso)` for a species lighter than the air, with the binary diffusion parameter `b = 3 / (8 sigma_ij^2) sqrt(k T / (2 pi mu_ij))` (hard spheres, kinetic diameters in `KINETIC_DIAMETER_PM`, Blanc's law for the mixture). The test checks the H2-in-N2 coefficient against 0.68 cm2/s at 273 K within 20 percent.
- **Energy limit** (`energy_limited_mass_loss_kg_s`): the total mass loss of all species is at most `0.15 pi F R_exo^3 / (G M)`; if the species demand more, all are scaled down together.
- **History** (`kept_shares`): described in section 1. The flux at an age is the present flux times `activity.coronal_xuv_at(age) / coronal_xuv_at(now)`, held to 4.9 Earth units. The present flux is `subject["xuv_flux_earth"]` if the subject has it, else the star's corona over the body's distance (from the luminosity and the insolation), else zero. The loss of the first 30 Myr and of the hotter phases above 4.9 times is not integrated here: it belongs to the retention rows d1 to d5 (GEN.511), which `surface_pressure` and the fixed point apply, so it is not counted twice.
- **What `mixture_escape.apply` does with it**: the share each species keeps multiplies its mole fraction; a species of which less than 1e-3 is left is gone; what stays is renormalised; if less than 1e-2 of the gas is left the body is airless (both are numerical cutoffs, rows AM22 and AM23). The body's mass, radius (from `radius_km`, else from the escape speed and the mass), pressure, temperature and age come from the subject; no new required input was added.
- **Envelope** (`hydro_based_mass_loss_g_s`): `Mdot = exp(beta) F^alpha1 d^alpha2 R^alpha3 Lambda^(zeta + theta ln d)`, `Lambda = G M m_H / (k T_eq R)`, with the constants of the boil-off branch below `exp(Sigma)` and of the XUV-driven branch above it (Table 1). The grid of the fit is 1 to 40 Earth masses and 300 to 2,000 K (`HYDRO_BASED_MASS_RANGE_EARTH`, `HYDRO_BASED_TEQ_RANGE_K`); outside it the rate is an extrapolation.

## 4. Checks

Tests: `src/tests/test_thermosphere.py` and `src/tests/test_atmosphere_mixture_escape.py`.

| Body | Model | Measured or published | Test |
|---|---|---|---|
| Earth exobase temperature, CO2 4e-6 to 0.4 | 1,290 to 288 K | 1,310 to 290 K (Johnstone Fig. 12) | within 8 percent at each of six points |
| Earth at 4.9 times the XUV | 4,820 K | 4,600 K | within 6 percent |
| Earth exobase altitude | 441 km now, 2,160 km at 4,600 K | 475 km; 2,000 to 2,500 km | 15 percent now; 1,800 to 2,600 km |
| Earth air history | N2, O2, Ar, CO2 kept; H2 and He lost | the Earth has them and has lost its primordial H2 and He | tested |
| Venus exobase temperature | 274 K (XUV 1.9 Earth units, CO2 0.965) | dayside about 285 to 300 K, global mean lower | 200 to 320 K; the whole air kept |
| Mars exobase temperature | 186 K (XUV 0.43), 11 K per mW/m2 of global-mean EUV | 180 to 260 K, about 9.5 K per mW/m2 of global mean (MAVEN, as recalled, tier 3) | 180 to 260 K; slope 4 to 19 |
| Mars air history | CO2, N2 and Ar kept; H2 lost | thin CO2 air | tested |
| Titan exobase temperature | 164 K (mesopause 150 K, XUV 0.011) | 150 to 190 K | tested |
| Titan hydrogen | diffusion limit 1.0e28 per second; Jeans 1.3e27 per second at the 723 km exobase | 9.2e27 per second, limited by the diffusion rate (Zhang 2024) | within a factor of 5 of each |
| Sub-Neptune and hot Jupiters, Table 2 of Kubyshkina et al. | ratio to the published rate 0.95 to 2.3 | nine planets, 1.9e9 to 4.5e10 g/s | within a factor of 2.5 for each |

A pure N2, O2 and Ar air with no CO2 on the Earth has a thermosphere of 1,300 K and more; it still keeps its nitrogen (the same finding as Johnstone et al.). A hydrogen-rich draw (40 to 70 percent H2 on an Earth-mass body) loses the hydrogen first; the first version of the history kept the initial mean mass throughout and made every such draw airless, which is why the history follows the air that is left.

## 5. Limits and what is not done

- **Non-thermal escape** (sputtering, photochemical and ion escape) dominates Mars's and Venus's loss of hydrogen and oxygen and is not modelled; the Earth model of Johnstone et al. says the same of the Earth. The thermal loss here is a lower bound for light species.
- **Cooling by gases other than CO2.** CH4 and HCN cool Titan's thermosphere, water cools a wet air; only CO2 is in the cooling term. Titan is therefore checked at the present XUV and with the observed mesopause, and its history is not integrated (the `titan` family is not filtered). A water- or methane-rich air is hotter in the model than in nature.
- **Day and night.** The temperature is the dayside of a fast rotator, as in the Earth model. Venus's night side and a tidally locked body are not separated.
- **Magnetic shielding** is GEN.804's (it wires the stored field to the loss).
- **Envelope escape is not wired.** `hydro_based_mass_loss_g_s` is checked and ready, but the sub-Neptune rows (`radius_valley`, `runaway_history`) still use the energy-limited rate of GEN.639 with its own efficiency. Wiring it is a follow-up item (GEN.738b, to be filed by TODO); it changes stored radii and needs an epoch.
- **The 4.9 times cap.** The flux is held at the top of the range the temperature fit was checked in. A star more active than that for longer (an M dwarf) is a longer time at the cap, not a hotter one.
- **Diameters and the cross section** are recalled values (tier 3); the hydrogen check on Titan is the test of the product.
- **The surface mixture is the exobase mixture.** The model does not separate the species by diffusion above the homopause and does not cold-trap water.

## 6. Registry

Rows AM63 to AM77 of `docs/design/research/gap-filling-plan.md` (and the generated `stand-in-values.md`) list each constant of the module with its source tier and tweak level. Rows AM19 to AM21 are replaced; AM22 and AM23 stay as numerical cutoffs.
