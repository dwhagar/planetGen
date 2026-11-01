# Part C: volatile-rich and gas bodies (ocean and ice worlds, sub-Neptunes, ice and gas giants, brown dwarfs, rogue planets)

Agent C, 2026-10-17. Read-only research; nothing in the repo was touched. Scripts, parsed tables and results are in `/mnt/project-files/design-drafts/parts/C-volatile-and-giants-files/` (see its README; the scripts point at scratch paths and need editing to rerun).

Labels: READ, 2ND-HAND, COMPUTED, SNIPPET, RECALLED, NOT REACHED as in the brief. Network note: after the widening, arXiv, PubChem, Wikipedia, Nature, IOP landing pages worked; IOP full-text PDFs, ADS, Science, A&A, MNRAS and ScienceDirect did not ([R2726] gave only its abstract; the PDF is behind a bot check, which I did not try to get around).

The GEN.425 report (`/mnt/project-files/reports/bulk-composition-sources-check-2026-10-17.md`, cited below as **G425 s.N**) already READ and verified [R3992, R4881] (method), [R4880, R4882] `f(x)`, [R2741] (Eq. 1, 3, tables, CSV), [R1429: Eq. 7 and 8], [R2865] (ice/rock 1 +/- 0.5), [R3193] snow line, [R0709, R3932]. I did not re-verify those; I re-used their numbers where noted and built on them. Where my own work agrees or differs it is stated in section 2.8.

---

## 1. Scope and interfaces

**Covers.** Where the volatile (ice, water) mass fraction of a body comes from; the water and ice phases and their equations of state; mass, radius and ocean structure of water worlds and sub-Neptunes; composition by mass of ice giants, gas giants and hot Jupiters; brown dwarfs; rogue planets as bodies. Two generator routes (keep drawn bulk figures and fill a fitting composition; draw composition first and derive radius) with a recommendation, error bars measured against a reference I built.

**Takes from other parts / the repo (interface values I need).**

| Need | From | Name used here |
|---|---|---|
| Rock+iron core radius, core mass fraction, rock constants and EOS | Part B (file did not exist when I wrote; I took [R3992: Table 1] constants, which G425 s.2.7 verified) | `CMF` (iron mass share of the *dry* part), `R_dry(M, CMF)` |
| Stellar [Fe/H], mass, luminosity, age | stored on stars (schema v95) | `feh`, `M_star`, `L_star`, `age` |
| Snow line | `formation.snow_line_au(L)` = 2.7 AU * sqrt(L/Lsun) (170 K blackbody, G425 s.2.12) | `r_snow` |
| Surface gravity, temperature, pressure | existing planet fields | `g`, `T_surf`, `P_surf` |
| Hydrosphere state | `physics/hydrosphere.py` stores `water_mass_fraction`, `ocean_fraction`, `ocean_depth_km`, `ice_shell_km`, `hp_ice_km`, `ocean_class`, `ocean_ph`, `water_activity`, `phosphorus`; giants (`body_type == "g"`) have all fields `None` | `x_w` = `water_mass_fraction` (whole-body mass fraction) |

**Hands over.**
- To the atmosphere lane: for giants the envelope heavy-element mass fraction `Z_env`, helium mass fraction `Y_env`, interior (1-bar) temperature and surface gravity; for water worlds the water inventory (`x_w`), the pressure at the base of the volatile layer `P_base` and whether the surface is steam, liquid or ice; for sub-Neptunes the envelope mass fraction `f_env`. No gas mixes, hazes or colours here.
- To the radiation lane: nothing from this part is a radioactive heat source except the rock core (Part B); water and H/He are radiogenically empty. Giants and brown dwarfs have no surface; surface radiation is atmosphere only.
- To the hydrosphere module: a better rule for `water_mass_fraction` (section 3.2), compressed layer depth and base pressure (2.6), and constants it currently guesses (2.5).

**Does not cover.** Atmosphere chemistry (other lane); rock/iron EOS and radiogenic heat (B); stellar and star-remnant composition (other parts).

---

## 2. Findings

### 2.1 Condensation sequence and snow lines

Saturation curves for the ices are given by [R2865] ([R2865: Tables 3 to 5]; the CO, CO2, CH4, H2S, N2, NH3, CH3OH curves are the [R1445] compilation, `ln P(bar) = sum_{i=0..5} A_i / T^i`; water uses an Eq. 12 and 13 form with seven coefficients). I recomputed the condensation temperature of each species from those curves for the nebula abundances of their Table 2 (CO:CO2 = 1:1 column, molar X/H2 in the second column). Each row is **COMPUTED** (inputs READ from Marboeuf Tables 2 to 4; method: solve `P_sat(T) = (X/H2) * P_H2`).

| Species | X/H2 (Marboeuf Table 2) | T_cond at P_H2 = 1e-12 | 1e-10 | 1e-8 | 1e-6 | 1e-4 bar | Source | Label |
|---|---|---|---|---|---|---|---|---|
| H2O | 5.95e-4 | 116.7 K | 128.1 | 141.9 | 158.9 | 180.5 | [R2865: Eq. 12], 13, Table 4 | COMPUTED |
| CH3OH | 8.92e-5 | 100.6 | 109.6 | 120.4 | 133.5 | 149.8 | Table 3 | COMPUTED |
| NH3 (monohydrate not modelled) | 4.16e-5 | 69.7 | 76.2 | 84.0 | 93.5 | 105.5 | Table 3 | COMPUTED |
| CO2 | 1.19e-4 | 60.7 | 66.3 | 73.1 | 81.5 | 92.2 | Table 3 | COMPUTED |
| H2S | 1.19e-5 | 52.1 | 57.2 | 63.3 | 71.0 | 80.7 | Table 3 | COMPUTED |
| CH4 | 3.57e-5 | 24.2 | 26.7 | 29.8 | 33.8 | 38.9 | Table 3 (valid from 20.6 K) | COMPUTED |
| CO | 1.19e-4 | 19.4 | 21.4 | 23.9 | 27.1 | 31.2 | Table 3 | COMPUTED |
| N2 | 4.16e-5 | 17.9 | 19.7 | 21.8 | 24.5 | 28.0 | Table 3 (valid from 10 K) | COMPUTED |

Local H2 pressure for orientation (COMPUTED, `P = n k T`): n_H2 = 1e10, 1e12, 1e14 cm^-3 at 100 K gives 1.4e-10, 1.4e-8, 1.4e-6 bar.

| Statement | Value | Source | Label |
|---|---|---|---|
| Evaporation (snow-line) temperatures (density-dependent, binding energies 5800, 2000, 850 K) | CO 18 to 22 K (adopt 20), CO2 42 to 52 (47), H2O 120 to 150 (135); carbon grains over 150 K; silicate about 1500 K | [R3322: Table 1] | READ |
| Average disk midplane temperature law | `T = T0 (r / 1 AU)^-q`, T0 = 200 K, q = 0.62 ([R0115] average) | same, Eq. 1 | READ |
| Gas-phase C/O rises outward through the H2O line, max about 1 between the CO2 and CO lines at 10 to 40 AU; interior to both lines the gas has stellar C/O | qualitative | same, section 2 and 4 | READ |
| H2O line at 135 K in that disk | `(200/135)^(1/0.62)` = 1.9 AU | Oeberg Eq. 1, my arithmetic | COMPUTED |
| Condensation temperatures from computed curves vs Oeberg | H2O agrees (128 to 142 K at n 1e10 to 1e12 vs 120 to 150). CO within 2 K (21 to 24 vs 18 to 22). CO2 differs: 66 to 73 K vs 42 to 52 K (Oeberg's CO2 binding energy 2000 K is lower than the Fray and Schmitt curve implies). Treat CO2 line temperature as uncertain by 20 K | my comparison | COMPUTED |
| Water line vs disk heating | `r_ice = 9.2 AU (alpha/0.01)^0.61 (Sigma_g,0 / 1000 g cm^-2)^0.8 (f_DG/0.01)^0.37` for a 1 Msun star ([R3935] fit; reaches inside 1 AU within the first Myr in most models) | [R0405: Eq. 1] | 2ND-HAND |
| Water line position in a population synthesis | median 2.81 AU; planetesimal density 3.2 g/cm3 inside, 1.0 g/cm3 outside | [R1268, R1269] (NGPPS II), [R1269] | READ |
| Snow line, viscous disk (G425 s.2.12): 3.3, 2.6, 1.9, 1.2, 0.5 AU for 1.0, 0.8, 0.6, 0.4, 0.2 Msun (small grains), 1-sigma about 0.2 dex | | [R3193] | READ via G425 |
| Water line temperature range quoted | 140 to 180 K typical disks; Ciesla uses 140 K | [R0803] | READ |
| Repo constant | 2.7 AU sqrt(L/Lsun) = the 170 K blackbody value (G425 s.2.12) | `physics/constants.py:92` | READ via G425 |

Consequence: the repo's single 170 K line is the warm end. For pressures typical of a mid-disk (n_H2 1e10 to 1e12 cm^-3) water condenses at 128 to 142 K, which puts the optically thin line at 2.7 * (170/135)^2 = 4.3 AU, while accretion-heated viscous disks give 3.3 AU (Mulders) and the synthesis median is 2.8 AU. A 0.2 dex scatter around 2.7 AU covers all of these (a parametrisation).

### 2.2 Water and ice mass fraction from formation

| Quantity | Value | Source | Label |
|---|---|---|---|
| Ice mass fraction of condensed solids beyond the snow line, all ices, range over models and 1 to 30 AU | 10 to 60 wt% (ice/rock 0.1 to 1.5); irradiated disks 10 to 50 wt% (0.1 to 1.1); refractory organics 20 to 30 wt%, minerals 40 to 65 (45 to 70) wt%; without organics ices are about 1.5 times higher; strong gradient near the water line | [R2865] section 4.2 | READ |
| Ice/rock, model-averaged | 1 +/- 0.5 (ice 45 +/- 15 wt%) without organics; 0.35 +/- 0.25 (30 +/- 10 wt%) with | Marboeuf abstract (via G425 s.2.12) | READ via G425 |
| Nebula ice species, molar vs H2O: CO 20 %, CO2 20 %, CH4 6, H2S 2, N2 7, NH3 7, CH3OH 15 (CO:CO2 = 1:1) | | Marboeuf Table 2 | READ |
| Water share of the ice mass if every ice is condensed | 0.429 by mass (0.565 by moles) | Marboeuf Table 2 times molar masses | COMPUTED |
| Water-only, ice/rock calibrated to 1.00 when all ices present: water share of solids by cumulative set (P_H2 1e-10 bar order H2O 128 K, CH3OH 110, NH3 76, CO2 66, H2S 57, CH4 27, CO 21, N2 20) | H2O only 0.300; +CH3OH 0.352 (ice 0.35); +NH3 0.364; +CO2 0.438; +H2S 0.444; +CH4 0.451; +CO 0.488; +N2 0.500 are the **ice fraction of solids**; the water mass fraction of solids is 0.300 / 0.304 / 0.294 / 0.245 / 0.241 ... / 0.214 | my cumulative sum of Marboeuf Table 2 | COMPUTED (assumption: ice/rock = 1 when all eight ices are present, the Marboeuf mean; ice/rock 0.5 and 1.5 give water 0.143 and 0.257 of solids) |
| Convention: solids beyond the water line are about 50 % water by mass ([R2726] chemical equilibrium); comets and icy satellites 40 to 50 %; carbonaceous chondrites 5 to 10 % (hydrated, after water loss); Ceres at least 20 to 30 % | | [R0803] section 2 | READ (Lodders cited, original NOT REACHED beyond abstract) |
| Stoichiometric estimate | water / (water + refractories) = 0.566 (all O not in silicate counted as water) | G425 s.2.12 | READ via G425 (COMPUTED there) |
| Rock:ice assumed for icy cores | x between 1/2 (solar rock:H2O) and 2/3 (rock:all ices); rock to ice ratio "does not strongly depend on host metallicity" | [R4880] | READ |
| Pebble-accretion cores that start beyond the ice line stay water-rich, `f_ice` about 0.5; inside they are dry (0); rocky core masses peak near 3 Mearth, max about 5; icy cores 1 to 36 Mearth. With collisions: about 33 % (20 %) of planets end with `0.05 <= f_ice < 0.45` (mean, max collisional model), 30 % stay 0, 37 % (50 %) about 0.5 | | [R4528], section 3, Appendix C and D | READ |
| Water-rich sub-Neptune scenarios: formed beyond the ice line gives "up to 50 %" water; formed inside, with a giant core blocking the pebbles, gives dry; giant forming inside the ice line gives a water-rich atmosphere and a few percent total water | | [R0405] section 2 and 4 | READ |
| Water inside the snow line, N-body: innermost planets around low-mass stars 1.2 to 3.3 % when planetesimals carry 50 % water, versus almost none in the [R3683] comparison | | [R0803] section 4 | READ |
| Water delivered to terrestrial planets inside the snow line by a forming giant's scattered planetesimals: inner-belt ordinary chondrites at most 0.1 % water; C-type asteroids 5 to 20 % water; Earth's water consistent with 2.5e-3 to 25e-3 Mearth of C-type material; "whenever a giant planet forms it invariably pollutes its inner system" | | [R3681], abstract and section 8.1 | READ |
| Delivered water fraction implied for a 1 Mearth planet (my arithmetic: (2.5e-3 to 25e-3) * (0.05 to 0.20)) | 1.2e-4 to 5e-3 | | COMPUTED |
| Synthetic planets in the TRAPPIST-1 zone around a 0.1 Msun star (M over 0.1 Mearth, a under 0.1 AU, 2222 planets): most purely rocky, 14 % have ice fraction over 0.1; ice-rich close-in planets reach down to about 1 Mearth for 0.1 Msun and about 3 Mearth for 1 Msun; ice-rich composition about 50 % ice outside the ice line | | [R0602] (NGPPS IV), [R0602], sections 3 and 4.6 | READ |
| Planets of the NGPPS population: ice-poor/ice-rich boundary set at 1 % volatiles by mass (set so that water gives high-pressure ice at the bottom of a 1 Mearth ocean) | | [R1268, R1269] (NGPPS II) | READ |
| Observed: three populations around M dwarfs (34 planets with mass error under 25 % and radius error under 8 %): rocky mean density 0.94 +/- 0.13 and radius 1.21 +/- 0.28 Rearth; water-rich 0.47 +/- 0.05 and 1.97 +/- 0.28; puffy 0.24 +/- 0.04 and 2.85 +/- 0.63 (densities relative to a scaled Earth-like model, CMF 0.325); density gap at 0.65; no intermediate water fractions in the sample; rocky planets up to about 10 Mearth, water worlds from 2 to 3 Mearth (2 Mearth for 0.3 to 0.5 Msun hosts), puffy from 6 Mearth; innermost planets of multis are rocky | | [R2774] | READ |
| Rocky:icy core number ratio about 1:2 for 1 to 4 Rearth planets; water worlds 2 to 4 Rearth with more than 1/4 and possibly over 1/2 water by mass | | [R4880] | READ |

Cross-checks and disagreements:
- The water fraction of the solids beyond the line is **0.14 to 0.39** (Marboeuf-based, my emulation, mean about 0.21 to 0.30) in one family of sources and **about 0.5** in another (Lodders via Ciesla, Zeng, Venturini, Bitsch, Burn). The factor of two is not resolved by any source read. I carry both as a range in section 3.
- "Dry inside the snow line" is a scenario, not a law: Venturini has exactly zero by construction (pebble accretion), Bitsch scenario A dry, scenario B a few percent, Ciesla 1 to 3 %, Raymond 1e-4 to 5e-3. The statement "dry inside, up to about 50 percent outside" is supported only in the narrow form: **cores that grow beyond the ice line have about 50 % ice in pebble models (Venturini, Bitsch), and cores inside are dry unless water is delivered**.
- Metallicity: Zeng says the rock-to-ice ratio does not depend strongly on [Fe/H]; metallicity changes how much solid there is (the repo already scales disk solids with `10^[Fe/H]` in `formation.disk_surface_density_scale`), not the water share. 
- Stellar mass: the minimum mass of an ice-rich close-in planet follows the two NGPPS IV anchors (1 Mearth at 0.1 Msun, 3 at 1 Msun), which a power law `M_ice,min = 3 (M_star/Msun)^0.48 Mearth` joins (COMPUTED two-point fit) and which gives 1.7 to 2.2 Mearth at 0.3 to 0.5 Msun, matching Luque and Palle's 2 Mearth.

### 2.3 Ice phases and the water equation of state

Densities, bulk moduli and their pressure derivatives from the [R2204] Gibbs-energy B-splines as distributed in SeaFreeze 1.1.3 (PyPI, GPL-3; **I used the spline data only and wrote my own evaluator, `sf_eval.py`; no package code was run**). `rho = 1 / (dG/dP)`, `K_T = -(dG/dP)/(d2G/dP2)`, `K' = dK/dP` by finite difference. Spline domains: Ih 0 to 400 MPa and 1 to 301 K; II 0 to 900 MPa; III 0 to 500 MPa; V 0 to 1000 MPa; VI 0 to 3000 MPa and 0 to 400 K; VII/X (French) 1.7 GPa to 1 TPa and 20 to 1800 K; liquid (Bollengier) 0 to 2300 MPa and 239 to 501 K. Original paper ([R2204], JGR Planets) not read. Every row COMPUTED from the spline data (2ND-HAND for the underlying measurements).

| Phase, state | P (MPa) | T (K) | rho (kg/m3) | K_T (GPa) | dK/dP | alpha (1/K) |
|---|---|---|---|---|---|---|
| Ih, melting point, 1 bar | 0.1 | 273.15 | 916.72 | 8.49 | 4.89 | 1.60e-4 |
| Ih, 200 K | 0.1 | 200 | 926.13 | 9.46 | 6.38 | 1.16e-4 |
| Ih, 100 K | 0.1 | 100 | 933.05 | 10.45 | 8.07 | 3.0e-5 |
| Ih at Ih-III-L point | 200 | 251.15 | 939.82 | 9.86 | 5.39 | 1.13e-4 |
| III at Ih-III-L | 209.5 | 251.15 | 1154.5 | 9.19 | 6.25 | 2.2e-4 |
| III at III-V-L | 355 | 256.43 | 1170.8 | 10.01 | 6.00 | 2.05e-4 |
| II (metastable region of its field) | 300 | 200 | 1196.0 | 14.89 | 5.95 | 1.14e-4 |
| V at III-V-L | 355 | 256.43 | 1244.9 | 12.92 | 6.15 | 1.8e-4 |
| V at V-VI-L | 618.4 | 272.73 | 1265.7 | 14.24 | 5.89 | 1.68e-4 |
| VI at V-VI-L | 618.4 | 272.73 | 1327.2 | 14.87 | 6.73 | 2.44e-4 |
| VI | 1000 | 300 | 1351.1 | 16.64 | 6.34 | 2.26e-4 |
| VI at VI-VII-L | 2216 | 355 | 1422.5 | 22.44 | 5.53 | 1.78e-4 |
| VII at VI-VII-L | 2216 | 355 | 1476.4 | 26.57 | 2.11 | 2.6e-4 |
| VII | 3000 | 300 | 1537.6 | 29.29 | 2.37 | 2.0e-4 |
| VII | 10000 | 300 | 1845.3 | 50.85 | 3.53 | 1.0e-4 |
| VII/X | 50000 | 300 | 2674.6 | 204.8 | 5.11 | 2.0e-5 |
| X | 100000 | 300 | 3107.1 | 503.8 | 4.62 | 7.4e-6 |
| X | 500000 | 300 | 4669.7 | 1708 | 2.79 | 1.6e-6 |
| X | 1000000 | 300 | 5814.4 | 2880 | 1.91 | 7.8e-7 |

Checks (COMPUTED): ice Ih at 273.15 K 916.72 kg/m3 matches the 0.9167 to 0.9168 g/cm3 on the Wikipedia ice page (SNIPPET: Wikipedia, 0.9340 at 93 K, spline 0.9330 at 100 K); liquid water 996.56 at 300 K and 1 bar. Frank, Fei and Hu 2004 formula for ice VII (READ in [R4881]) `rho(mol/cm3) = 0.0805 + 0.0229 (1 - e^(-0.0743 p)) + 0.1573 (1 - e^(-0.0061 p))`, p in GPa, gives 1.584 g/cm3 at 3 GPa against 1.538 here (3 %, their 300 K and mine). Seager's DFT Table 2 (cold, READ) gives 2.576 at 44.3 GPa and 3.326 at 126.8 GPa, 3.485 at 155.9; the spline 2.67 at 50 GPa, 3.11 at 100 GPa: agree within about 3 %. Seager's ice VII zero-pressure BME (K0 23.7, K0' 4.15, rho0 1460, READ) vs the spline's VII at 3 GPa (K 29.3, K' 2.4): the BME is a different representation of the same compression curve; use the spline values above or Seager's polytrope, not both.

**Phase boundaries, recomputed (COMPUTED from the same splines, `G(ice) = G(liquid)`):** Ih melting 273.15 K at 0.1 MPa, 264.2 K at 100 MPa, 252.2 K at 200 MPa; III 251.2 K at 210 MPa, 256.2 K at 350 MPa; V 256.2 K at 350, 271.7 K at 600 MPa; VI 272.1 K at 620, 286.8 K at 800, 300.4 K at 1000, 327 K at 1500, 347 K at 2000, 353.8 K at 2200 MPa.

| Triple point | [R4881] (READ; [R0777]) | `hydrosphere.py` `LIQUIDUS_K_MPA` [R4585] | Recomputed here |
|---|---|---|---|
| Ih-III-liquid | 209.5 MPa, 251.15 K | 207.5 MPa, 251 K | about 210 MPa, 251.2 K |
| III-V-liquid | 355.0 MPa, 256.43 K | 346 MPa, 256 K | about 350 MPa, 256.2 K |
| V-VI-liquid | 618.4 MPa, 272.73 K | 626 MPa, 273.3 K | about 620 MPa, 272.1 K |
| VI-VII-liquid | 2.216 GPa, 355 K | 2200 MPa, 355 K | about 2.2 GPa, 353.8 K |
| VII-X-liquid | 47 GPa, 1000 K | not in the module | not modelled |

The three sets agree within 1 to 3 %. The module's linear interpolation between triple points puts the 280 K liquidus at 755 MPa; the spline puts it near 715 MPa (between 620 MPa/272.1 K and 800 MPa/286.8 K), a 5 % difference that is smaller than the uncertainty in the bottom temperature.

**Constants in `hydrosphere.py` vs these numbers (COMPUTED):** `ICE_DENSITY_KG_M3 = 917` fine for Ih; `HP_ICE_DENSITY_KG_M3 = 1300` is below ice VI (1327 at its lowest pressure to 1422 at 2.2 GPa) and ice VII (1476 and up), so `hp_ice_km` is overestimated by about 4 to 10 %. `DEEP_WATER_DENSITY_KG_M3 = 1050` is the mean for about 0.2 GPa; the liquid column down to the 280 K liquidus is denser (1044 at 100 MPa, 1110 at 300, 1180 at 600, 1218 at 800 MPa at 280 K, so a mean of about 1110 to the 715 MPa liquidus), 6 % above the constant, which makes `max_liquid_depth_km` (73 km at Earth gravity) about 4 % too large.

**Superionic and X:** ice VII transforms to ice X near 47 GPa and 1000 K (as cited by [R4881], READ); [R1457] DFT phase diagram (superionic phase, used by [R4881] from 37.4 GPa to 8.9 TPa, READ as described there; original paper not read, 2ND-HAND). Experiment: laser-shock water solidifies to a bcc (probably ice X) then face-centred-cubic superionic **ice XVIII** at **100 to 400 GPa and 2000 to 3000 K**; superionic conduction (protons mobile in an oxygen lattice) is predicted over 100 GPa and over 2000 K ([R3081]:251, abstract READ via nature.com landing page; paper body paywalled). Superionic ice is a dense solid water layer in the deep interior of water worlds and ice giants; the generator needs only its density (the spline "VII_X_French" above, valid to 1 TPa and 1800 K) and the fact it is solid.

**Low-pressure ices of other species (small icy bodies, comets, moons):**

| Material | Value | Source | Label |
|---|---|---|---|
| CO2 solid (dry ice) | 1.56 g/cm3 at -79 C (about 194 K) and 1.512 at -56.6 C; 1600 kg/m3 used for cometary CO2 ice | PubChem density annotation (CRC and USCG entries); Blum and Bischoff, [R0430] | READ (PubChem, arXiv text) |
| H2O ice for comets | 934 kg/m3 (at low T) | Blum and Bischoff | READ |
| N2 solid | 1.0265 g/cm3 | PubChem N2 density entry | READ (temperature not given there) |
| CO solid, hexagonal | 929 kg/m3 at 65 K; liquid 788.6 at 81.63 K | PubChem CO | READ |
| CH3OH liquid | 0.792 g/cm3 (liquid at 20 C); the solid is not given | PubChem | READ (liquid only) |
| CH4, NH3, H2S solid densities | no reachable primary value found. Crystallographic estimates from recalled lattice constants would be about 0.52 (CH4 fcc phase I, a about 5.9 A), about 0.82 to 0.85 (NH3, a about 5.1 A), but these are RECALLED and not verified | RECALLED | RECALLED |
| CH4 clathrate hydrate | about 0.9 g/cm3 (908 to 917 kg/m3 depending on structure and cage filling; 0.92 for natural hydrate extrapolated to 273 K) | web search snippets (Wikipedia, MarineInfo) | SNIPPET |
| Clathrate cage and composition | hydrate number `n_hyd = 6` used for all (ideal 5.75 for structure I, 5.67 for structure II); CO, CH4, CO2, H2S form structure I, N2 structure II; NH3 forms a stoichiometric monohydrate, not a clathrate | [R2865] section 3.4 | READ |
| Clathrate stability, `ln P_cl(Pa) = A0 + A1 / T`: CO 23.3818, -1890.39; CO2 23.294, -2565; CH4 22.6864, -2176.2; H2S 22.8931, -3111.02; N2 22.719, -1696.05; NH3 19.515, -2878.23 | | Marboeuf Table 5 (from [R1925]) | READ |
| Refractory material in comets | density 1800 to 3500 kg/m3 (organics 1000 to 2000, minerals 2600 to 5000); bulk packing density of 67P Phi_bulk 0.15 to 0.37 | Blum and Bischoff section on porosity | READ |
| Planetesimal densities in a planet-formation model | 3.2 g/cm3 inside the ice line, 1.0 g/cm3 outside | NGPPS II | READ |

For bodies below about 0.01 Mearth, G425 s.3.2 recommends a volume-additive density; with the numbers above, an ice-rock small body of ice fraction `w` has `rho = 1 / (w/rho_ice + (1-w)/rho_rock)` times a packing density (0.15 to 0.37 measured for a comet; no value was found for other small bodies).

### 2.4 Mass-radius of solid and water-rich bodies

**Reference I built** (`structure.py`, COMPUTED): hydrostatic integration, `dm/dr = 4 pi r^2 rho`, `dP/dr = -G m rho / r^2`, surface pressure 1 bar, shooting on central pressure; layers iron ([R3992: Table 1], epsilon-Fe Vinet, 8.30 g/cm3, K0 156.2, K0' 6.08), silicate (Seager enstatite BME below 23.83 GPa, then MgSiO3 perovskite BME4 247 GPa, 3.97, K0'' -0.016 per GPa, 4.10 g/cm3), water (the 300 K isotherm of the Journaux splines: liquid to 1 GPa, ice VI to 2.2 GPa, VII/X to 1 TPa, then Seager's Table 3 polytrope matched in density). Above 1 TPa the Fe and silicate branches switch to the Seager Table 3 polytropes (`rho0 + c P^n`, scaled to be continuous) in place of Thomas-Fermi-Dirac. No thermal expansion, no mixing.

Validation against the six two-layer curves tabulated in [R4881: Table 1] (READ; parsed to `zs13_table1.json`; their curves use TFD above the Vinet/BME range and perovskite throughout the silicate):

| Composition (mass) | M = 0.1 | 0.5 | 1 | 2 | 5 | 10 | 20 | 50 Mearth |
|---|---|---|---|---|---|---|---|---|
| Pure H2O | +0.4 % | -0.3 | -0.7 | -1.2 | -1.2 | -0.7 | +0.2 | +1.3 |
| 25 % MgSiO3 + 75 % H2O | +1.4 | -0.2 | -0.5 | -0.8 | -0.8 | -0.5 | +0.4 | +1.6 |
| 50 % MgSiO3 + 50 % H2O | +2.6 | +0.2 | -0.2 | -0.4 | -0.4 | -0.2 | +0.6 | +1.8 |
| Pure Fe | -0.1 | -0.1 | -0.1 | -0.1 | 0.0 | -0.2 | -0.8 | -2.9 |
| 50 % Fe + 50 % MgSiO3 | +5.5 | +2.8 | +1.9 | +1.4 | +1.1 | +1.0 | +1.0 | +0.5 |
| Pure MgSiO3 | +7.7 | +3.5 | +2.4 | +1.9 | +1.5 | +1.6 | +2.1 | +2.8 |

(My radius over theirs; positive means mine is larger.) Water-bearing cases agree within 1.5 % from 0.5 to 50 Mearth; the silicate-bearing dry cases are 2 to 3 % large at 1 Mearth and up to 8 % at 0.1 Mearth because I use enstatite (3.22 g/cm3) in the upper mantle where Zeng and Sasselov use perovskite; the Earth-like case then gives 0.989 Rearth at 1 Mearth (30 % Fe, 1 bar) against [R4882]'s 1.005 (1.6 % low), the same few-percent spread G425 s.3.2 reports among published closed forms. Pure water at 1 Mearth: 1.401 Rearth against [R4882]'s 1.41 (G425 s.2.4) and 1.42 interpolated from Zeng and Sasselov Table 1.

**Three-layer grid** (`grid3.py`, `grid3.json`, 175 integrations, COMPUTED): M = 0.3, 0.5, 1, 2, 4, 8, 16 Mearth; `CMF` of the dry part = 0, 0.1, 0.325, 0.5; water mass fraction of the whole body `x_w` = 0, 0.05, 0.1, 0.25, 0.5, 0.75, 1 (CMF irrelevant at 1). Earth-like dry radii (CMF 0.325): 0.699, 0.810, 0.989, 1.203, 1.452, 1.736, 2.049 Rearth. Ratio `R(x)/R(0)` at CMF 0.325 compared with the published `f(x) = 1 + 0.55 x - 0.14 x^2` ([R4880], READ):

| x | f(x) | reference ratio at M = 0.3 / 1 / 4 / 16 Mearth |
|---|---|---|
| 0.05 | 1.027 | 1.043 / 1.041 / 1.038 / 1.039 |
| 0.10 | 1.054 | 1.074 / 1.070 / 1.066 / 1.066 |
| 0.25 | 1.129 | 1.151 / 1.140 / 1.134 / 1.135 |
| 0.50 | 1.240 | 1.256 / 1.243 / 1.230 / 1.235 |
| 0.75 | 1.334 | 1.348 / 1.332 / 1.316 / 1.325 |
| 1.00 (CMF 0) | 1.410 | 1.437 vs the CMF 0.325 dry radius / 1.417 / 1.397 / 1.412 |

So the published `f(x)` is good to 0.5 to 2 % from x = 0.25 to 1 and low by 1.5 to 2.5 % for small x (0.05 to 0.1), against my reference. `f(x)` was defined on a CMF 0.325 base only; used on top of [R4882]'s CMF-dependent dry radius it errs by up to +10 % (CMF 0.5 or x near 1), because the iron factor and the water factor are not independent. A refit (below) fixes that.

**Recommended closed form** (my fit to the 175-point grid, COMPUTED; same structure as [R4880, R4882], which it reproduces in the limits):

`R / Rearth = (1.0688 - 0.2317 CMF (1 - x)) * M^(1/3.7) * (1 + 0.4492 x - 0.1597 x^2)`

M in Mearth, `x` water mass fraction of the body, `CMF` iron mass share of the dry part. Limits: dry Earth-like (CMF 0.325, x 0) gives 0.993 at 1 Mearth; the dry coefficients 1.0688 and 0.2317 equ[R4882]'s 1.07 and 0.21 within 10 % (an independent refit landing on nearly the same numbers); pure water gives `1.0688 * 1.2893 = 1.378` at 1 Mearth against the integrated 1.401 and Zeng's 1.41.
Errors against my reference: rms 1.2 %, maximum 3.6 % (M 0.3 to 16 Mearth, CMF 0 to 0.5, x 0 to 1); by mass band rms 1.06 % (0.3 to 1), 1.06 % (1 to 8), 1.42 % (8 to 16). Against the [R4881: Table 1] water mixtures (CMF 0): +2.7 / +0.1 / -0.3 / +1.2 % for 50 % water at 0.3 / 1 / 3 / 10 Mearth, -0.6 to +1.4 % for 75 % water, -2.4 to -0.9 % for pure water, and **growing beyond 20 Mearth** (+6 to +7 % at 30, +17 to +20 % at 100 Mearth, the power law bends because compression grows), so **valid 0.3 to 20 Mearth**. Dry cases (x = 0, vs Zeng and Sasselov, CMF 0.5): +3 % at 1 Mearth, +7 % at 0.3. Cold reference (300 K isotherm, 1-bar surface): hot steam or supercritical envelopes change this (below). Inverse (CMF, M and R known, solve for `x`) is a monotone one-dimensional root-find.

**Irradiated water worlds.** [R4448] show that a water-rich rocky planet more irradiated than the runaway-greenhouse limit (about 1.1 times Earth's insolation for a Sun-like star) holds its water as a thick steam atmosphere, which inflates the radius so that a given water mass fraction gives a much lower bulk density; they give an empirical steam-atmosphere thickness as a function of core gravity and radius, water mass fraction and irradiation (coefficient tables, not copied here). Example in the abstract: TRAPPIST-1 b, c, d can hold at most 2, 0.3 and 0.08 % water (READ abstract). Consequently the cold fit above is for planets *below* the runaway limit; above it the thickness law is needed. Also from the abstract of Luque and Palle (READ): rocky planets close to their stars may hold supercritical water that raises the radius. Both are inputs to a later flag, not implemented here.

**Ocean-world structure** (`ocean_table.py`, COMPUTED from the integrated reference, CMF 0.325, 300 K isotherm, surface 1 bar; the two x = 0.2 rows were recomputed on 2026-10-25 from `physics/water_layer.py` at the table's own radii (DOC.79), the earlier 1635 km / 26.0 GPa and 2519 km / 114 GPa did not follow from those radii, and `tests/test_water_layer.py` now reads this table and holds every row to the code):

| M (Mearth) | x_w | R (Rearth) | water layer thickness (km) | pressure at base (GPa) | surface g (m/s2) | `g x M / (4 pi R^2)` estimate (GPa) |
|---|---|---|---|---|---|---|
| 1 | 1e-4 | 0.9895 | 1.2 | 0.012 | 10.01 | 0.012 |
| 1 | 1e-3 | 0.9907 | 11.7 | 0.119 | 9.98 | 0.119 |
| 1 | 1e-2 | 1.0004 | 102 | 1.18 | 9.79 | 1.14 |
| 1 | 0.1 | 1.0583 | 746 | 11.1 | 8.75 | 9.1 |
| 1 | 0.2 | 1.1281 | 1294 | 19.2 | 7.70 | 14.2 |
| 1 | 0.5 | 1.2301 | 3005 | 45.4 | 6.48 | 25.0 |
| 5 | 1e-4 | 1.5406 | 2.4 | 0.051 | 20.64 | 0.051 |
| 5 | 1e-3 | 1.5428 | 22.8 | 0.51 | 20.58 | 0.505 |
| 5 | 1e-2 | 1.5571 | 173 | 5.03 | 20.21 | 4.87 |
| 5 | 0.1 | 1.6419 | 1189 | 48.4 | 18.17 | 39.4 |
| 5 | 0.2 | 1.7463 | 2024 | 84.7 | 16.07 | 61.7 |
| 5 | 0.5 | 1.8942 | 4623 | 205 | 13.65 | 111 |

The naive layer depth in `hydrosphere.water_layer_depth_km` (incompressible, 1000 kg/m3, planet radius fixed) is right to 2 % for `x_w <= 0.01` (104 km vs 102 for 1 Mearth) and **too deep by 40 % at x_w = 0.1** (1047 km vs 746 for 1 Mearth). `P_base ~ g x M / (4 pi R^2)` holds to 3 % for `x_w <= 0.01`, 20 % low at 0.1. For Earth gravity a liquid column deeper than about 70 km (about 6e-3 of an Earth mass, 26 Earth oceans) reaches the ice VI liquidus (715 to 755 MPa at 280 K): above `x_w` of order 6e-3 the ocean is sealed from the rock by high-pressure ice (the module's `hp_ice_km` rule), and the water layer thickness at `x_w = 0.1` already holds ice VII/X, which is solid water at 300 K (the integration uses 300 K, so deep layers are really hotter and may be fluid or superionic; Seager's note that temperature changes the density by under about 4 %, and radius by 1.2 % per 3.5 % in density, bounds this at the few-percent level for radius).

### 2.5 Sub-Neptunes: envelope on a rock or water core

All three rows are READ via G425 s.2.5 and s.2.6, plus two checks of my own.

| Item | Value | Source | Label |
|---|---|---|---|
| Rock/iron core | `R_core = (M_core/Mearth)^0.25 Rearth` (2:1 rock:iron, within about 2 %) | [R2741: Eq. 1] | READ via G425 |
| Envelope | `R_env = 2.06 Rearth (M/Mearth)^-0.21 (f_env/5 %)^0.59 (F/Fearth)^0.044 (age/5 Gyr)^-0.18`, total `R = R_core + R_env + R_atm`, `R_atm ~ 9 H`; M 1 to 20 Mearth, f_env 0.01 to 20 %, F 0.1 to 1000 Fearth, 100 Myr to 10 Gyr; radius at 20 mbar; rock/gas divide suggested 1.75 Rearth | Eq. 3 | READ via G425 |
| Fit accuracy against the paper's own tables (f_env at least 0.5 %) | mean +3.1 %, rms 6.3 % (1 Gyr, enhanced opacity); rms 8.3 % at 10 Gyr | G425 s.2.5 | COMPUTED there |
| **Inverse**: given M, R, F, age, find `f_env` from Eq. 1 and 3 (`M_core = M (1 - f_env)`), checked on the 432 table entries with enhanced opacity and f_env at least 0.5 % (CSV from G425) | rms **0.106 dex**, median absolute 0.069 dex, worst 0.39 dex (by age: 100 Myr 0.118, 1 Gyr 0.079, 10 Gyr 0.117 dex rms); same size as the paper's own "within about 0.1 dex" | my test of `check_lopez_fortney_2014.py` data | COMPUTED |
| Water-rich cores | no envelope-on-water-core fit exists in what was read. [R4880]: "a gas envelope of at most a few percent by mass" is enough to inflate any planet to 2 to 4 Rearth; envelope isotherms 300 to 1000 K; surface 1 mbar | [R4880] | READ via G425 |
| Transition | most 1.6 Rearth planets are not rocky ([R3800], not read); the observed M-dwarf gap is in density (0.65 of Earth-like), not radius | [R2774] | READ |

### 2.6 Gas and ice giants: composition by mass

| Quantity | Value | Source | Label |
|---|---|---|---|
| Heavy-element mass vs planet mass, 47 transiting giants with `20 Mearth < M < 20 MJ`, insolation under 2e8 erg s^-1 cm^-2 (no hot-Jupiter inflation) | `M_Z = (57.9 +/- 7.03) (M/MJ)^(0.61 +/- 0.08) Mearth` (units from Fig. 7: M in Jupiter masses, M_Z in Earth masses; M_Z roughly proportional to sqrt(M)); intrinsic spread: 1 sigma within a factor `10^sigma = 1.82 +/- 0.09` (0.26 dex); "using planet mass alone ... appears accurate to a factor of a few" | [R4352] (not 1511.07442), section 5.1 | READ |
| Enrichment relative to the star | `Z_planet / Z_star = (9.7 +/- 1.28) (M/MJ)^(-0.45 +/- 0.09)`, `Z_star = 0.014 * 10^[Fe/H]` (their simple approximation); weaker fit to [Fe/H] than earlier work: Kendall tau 0.088, p 0.38 for M_Z against [Fe/H] (no correlation in the full sample); planets over 100 Mearth of metals orbit metal-rich stars; few planets over 1 to 3 MJ around metal-poor stars; no relation with orbital distance, stellar mass or eccentricity | same, sections 5.2, 5.3 | READ |
| Model used for the above | [R3932] H/He with Y = 0.27, additive volumes; metals as a 50:50 rock-ice mixture (ANEOS), **distributed through the envelope** (not all in a core: an all-core model needs implausibly large cores); more than 50 Mearth of metals in many planets implies metal-enriched envelopes | same, section 3.1, 3.2, abstract | READ |
| Same relations tabulated (my arithmetic) | M (MJ) 0.05 / 0.1 / 0.3 / 1 / 3 / 10: M_Z = 9.3 / 14.2 / 27.8 / 57.9 / 113 / 236 Mearth; Z = 0.59 / 0.45 / 0.29 / 0.18 / 0.12 / 0.074; second relation at Z_star 0.014: Z = 0.52 / 0.38 / 0.23 / 0.14 / 0.083 / 0.048. The two relations differ by up to 30 % in Z; `Z = 1` at about 4.5 Mearth (M_Z = M) so the fit cannot be extended below that | | COMPUTED |
| Heavy-element mass in giants from earlier work | extended M_Z vs M trends and [Fe/H] link [R3077] | not re-read | NOT READ (fetched, 469 lines, not used) |
| Atmospheric metals vs mass | H2O abundance rises with decreasing mass: generally sub-stellar in gas giants, stellar to super-stellar in Neptunes and mini-Neptunes; the overall water trend is far below the solar-system carbon mass-metallicity relation; Na and K in gas giants are stellar or super-stellar and follow the solar-system trend; elements are not all enriched together | [R4663], abstract | READ (abstract) |
| Ice giants | H/He envelope about 10 to 20 % of the total mass; metals mass fraction 0.75 to 0.92 (Uranus), 0.76 to 0.90 (Neptune) in empirical density-profile models; interiors need not be water-rich; the 3-layer rock-water-gas models need rock-to-water ratios far from solar (water-to-rock 19 to 35 by one model), which is why they are doubted; possible composition gradients and non-adiabatic profiles | [R1892] | READ |
| Transitional planets 4 to 10 Rearth | ice-rich with gas envelopes of at least 5 to 10 % by mass, a few tens of Mearth; gas giants over 10 Rearth dominated by H2-He; cosmic H2-He mix 75 % H2, 25 % He by mass | [R4880] | READ |
| Protosolar mass fractions | X 0.7110, Y 0.2741, Z 0.0149 (photosphere 0.7491 / 0.2377 / 0.0133) | [R2726] abstract (IOP landing page) | READ (abstract); star lane owns Y(Z) |
| Giant radius, age and irradiation | four-dimensional grid: 17 to 3594 Mearth (0.052 to 11.3 MJ), core 0 to 100 Mearth, 0.02 to 9.5 AU around the Sun (Teq about 1960, 1300, 875, 260, 78 K), ages 0.3, 1, 4.5 Gyr; radii in RJ at 1 bar; SCvH H/He with Y = 0.28, olivine/ice/iron cores; 4.5 Gyr example values are in the correction paragraph below the table | [R1429], astro-ph/0612671, Tables 2 to 4 (parsed to `fortney2007_tables.json`) | READ |
| Radius vs mass shape: Terran `R = 1.008 Rearth M^0.279` up to 2.04 Mearth (rms 4.0 %), Neptunian `R ~ M^0.589` (+0.044 -0.031) to 0.414 MJ (rms 14.6 %), Jovian `R ~ M^-0.044` (+/-0.018) to 0.080 Msun (rms 7.4 %), stellar `M^0.881` beyond | | [R0758: Table 1] | READ |
| The repo's giant relation: Neptunian exponent 0.519 anchored on Saturn (95.16 Mearth, 9.449 Rearth), Jovian -0.04 anchored on Jupiter (317.8, 11.21); the code comment says the shape is Chen and Kipping's. The paper's own Neptunian exponent is 0.589 | | `physics/constants.py:182-183` | READ (code) |
| Hot-Jupiter inflation | heating efficiency (anomalous power as a fraction of incident flux) rises with `T_eq` to a maximum of about 2.5 % near 1500 K and falls to about 0.2 % at 2500 K (Ohmic dissipation predicted this); the Gaussian model's peak is about 2 % (favoured 75:25 over a 3.5 % peak); planets below about 0.5 MJ behave differently and were excluded | [R4351], abstract and Fig. 8 text | READ (abstract, text) |

Test of the table interpolation (COMPUTED, `fortney_fit.py`; interior grid points left out and predicted by linear interpolation, M at least 0.15 MJ): in log mass rms 2.0 %, worst 11.5 %; in log `T_eq` (distance, 3 interior of 5 rows) rms 3.7 %, worst 19.6 %; in log age (1 Gyr from 0.3 and 4.5) rms 1.0 %, worst 4.6 %. A single global regression (ln M, ln T_eq, ln age, core/M, cross terms) is worse: rms 4.1 %, worst 32 % (M at least 0.15 MJ), 2.6 % and 13.6 % for M at least 0.41 MJ. **The table is the better representation.**

Correction to the Fortney row above: for a 1 MJ planet at 4.5 Gyr the table gives, for core 0, 10, 25, 50, 100 Mearth: 1.081, 1.058, 1.026, 0.974, 0.878 RJ at 1 AU, and 1.019, 0.999, 0.972, 0.926, 0.842 at 9.5 AU (READ). Interpolating the Thorngren value `M_Z = 58 Mearth` as a core at 1 AU gives about 0.96 RJ. The 1-bar radius is smaller than the roughly 1.0 to 1.1 RJ seen for cool Jupiter-mass planets because Thorngren spreads the metals through the envelope, which gives a larger radius than the same metals in a core (his section 3.1): using his M_Z as a Fortney *core* therefore understates the radius by a few percent. I did not quantify it.

### 2.7 Brown dwarfs, deuterium and hydrogen burning, degenerate radius limit

| Quantity | Value | Source | Label |
|---|---|---|---|
| Deuterium burning mass (50 % of initial D burned), solar Z, Y = 0.25, D/H 2e-5 | 13.55 MJ; range 12.7 to 14.35 MJ for 10 to 90 % burned; for "most proto-brown dwarf conditions" `13.0 +/- 0.8 MJ` | [R4170: Table 1] and abstract | READ |
| Dependence | on Y: He22 14.08, He25 13.55, He28 12.83, He30 12.42, He32 12.05 (50 %): about `-0.2 MJ` per 0.01 in Y; on initial D: D1 (1e-5) 13.74, D2 13.55, D4 (4e-5) 13.09: `-0.8 MJ` per dex of D; on Z (Z/Zsun 0, 0.01, 0.1, 0.3, 1): 15.40, 14.56, 13.82, 13.33, 12.79 MJ (50 %): `-0.9 MJ` per dex of Z; with COOLTLUSTY atmospheres (0.3, 1, 3 Zsun): 13.77, 13.48, 13.13 (`-0.7 MJ` per dex). Whole range 11.0 MJ (3 Zsun, 10 %) to 16.3 MJ (Z = 0, 90 %). Their Zsun is 0.0189 (Anders and Grevesse) | same Table 1, 2 | READ |
| Hydrogen-burning minimum mass | 0.07 to 0.074 Msun at solar Z and Y = 0.25 (Teff 1700 to 1750 K, L about 6e-5 Lsun); 0.092 Msun at Z = 0 (Teff about 3600 K); `dM_edge/dY` about -0.1 Msun; weak function of metallicity, `M_edge ~ 0.0865 (kappa_R / 1e-2 cm2 g^-1)^(-1/9) Msun`; lithium burning below 0.06 to 0.065 Msun | [R0606], astro-ph/0103383, section III | READ |
| Cooling power laws for old objects (solar Z): `L ~ 4e-5 Lsun (t/Gyr)^-1.3 (M/0.05 Msun)^2.64 (kappa/1e-2)^0.35`, `Teff ~ 1550 K (t/Gyr)^-0.32 (M/0.05 Msun)^0.83 (kappa)^0.088`, `M ~ 35 MJ (g/1e5)^0.64 (Teff/1000 K)^0.23`, `t ~ 1.0 Gyr (g/1e5)^1.7 (1000 K/Teff)^2.8`, `R ~ 6.7e4 km (1e5/g)^0.18 (Teff/1000 K)^0.11`; for Z = 0 divide L by about 4 and use indices 1.25, 2.4 and Teff 1140 K with 0.31 and 0.77 | [R0606: Eq. 1 to 5] | READ |
| Degenerate-matter radius | on the very low-mass main sequence `R ~ M^0.6`; for a cold substellar object the radius peaks near 4 MJ [R4867] and is mass-independent within 30 % from 0.3 to 70 MJ because Coulomb (`R ~ M^(1/3)`) and degeneracy (`R ~ M^(-1/3)`) compete, "a feature of a polytrope of index 1.0" | [R0606] section II | READ ([R4867] itself NOT REACHED, not on arXiv) |
| Radius table (solar Z, Y about 0.275, dust-free COND models) | R in RJ at ages 0.1 / 0.5 / 1 / 5 / 10 Gyr (3 decimals; R/Rsun is quoted to 3 digits so these carry about 0.5 % rounding): 1 MJ 1.139 / 1.061 / 1.032 / - / -; 2.1 MJ 1.168 / 1.090 / 1.061 / 1.022 / -; 5.2 MJ 1.187 / 1.100 / 1.070 / 1.022 / 1.002; 10.5 MJ 1.168 / 1.070 / 1.041 / 0.983 / 0.963; 21 MJ 1.187 / 1.012 / 0.973 / 0.924 / 0.905; 42 MJ 1.285 / 0.973 / 0.905 / 0.827 / 0.808; 63 MJ 1.460 / 0.993 / 0.895 / 0.769 / 0.740; 75 MJ 1.576 / 1.041 / 0.924 / 0.808 / 0.788; 84 MJ 1.654 / 1.080 / 0.993 / 0.963 / 0.963; 105 MJ 1.839 / 1.246 / 1.216 / 1.216 / 1.216. Teff at 5 Gyr: 129 K (2.1 MJ), 322 K (10.5), 760 K (42), 1120 K (63), 1524 K (73), 2320 K (84), 2785 K (105 MJ). Minimum radius about 0.077 Rsun (0.74 RJ) near 0.06 Msun (63 MJ) at 10 Gyr | [R0253], astro-ph/0302293, Tables 1 to 5 (parsed to `baraffe2003_cond.json`) | READ |
| Mass cut of the giant classes | The repo's J class goes to 13 x 317.8 Mearth (13 MJ) and a giant is at most 1 % of the star's mass (`tuning.GIANT_*`); the 13 MJ cut moves by -0.9 MJ per dex Z and -0.2 MJ per 0.01 Y | | READ (code) with Spiegel |
| Rogue brown dwarf cooling used by the repo | `L ~ t^-1.3` from Burrows; Teff capped at 2800 K | `docs/design/systems/rogue-planet-surface.md` | READ |

### 2.8 Rogue planets as bodies

| Item | Finding | Source | Label |
|---|---|---|---|
| Composition | a rogue keeps the composition it formed with; ejection does not change the heavy-element budget of a rock-ice core. What changes: no stellar irradiation, so surface temperature is set by internal heat alone (`docs/design/systems/rogue-planet-surface.md`, steps 1 to 4: radiogenic heat for rock, Kelvin-Helmholtz cooling for giants, 2.725 K background, hydrogen envelopes kept with chance 10 % terrestrial and 50 % sub-Neptune, water 1e-4 to 10 % of mass in 35 % (terrestrial) or 80 % (sub-Neptune) of rogues) | the design note | READ |
| Population synthesis ejection statistics (model, 1000 systems per stellar mass): fraction of systems that eject planets 17.5 % (1 Msun), 12 % (0.3), 0.8 % (0.1); mean ejected mass 7.2, 5.2, 1.6 Mearth; median 0.8, 0.3, 1.4 Mearth; mean total ejected per ejecting system 44.2, 24.5, 3.3 Mearth; population = Jupiter-mass planets, then more numerous super-Earths, with a tail down to 1e-5 Mearth; the paper gives no composition of the ejected planets | Ma et al. 2016, [R2783] | READ |
| Stripping of volatiles by ejection | no source read states how much of a planet's water or envelope is lost on ejection. Ice-rich bodies come from beyond the snow line where planets and embryos are ejected preferentially (outer planets are the ones scattered), so the water fraction distribution of ejected planets should be taken from the *outer* branch (2.2) unless a source says otherwise; this is a parametrisation | none | NOT REACHED |
| Radius at fixed mass | a rogue giant keeps cooling: radius of a gas giant is nearly flat in mass and shrinks with age (Baraffe, Fortney tables) and no irradiation inflation applies, so use the cold (9.5 AU, Teq 78 K) rows of the Fortney tables | [R1429] | READ |

### 2.9 Agreement and differences with G425 and the existing repo

- Agree: Seager constants, [R4880] `f(x)` (I confirm it to 0.5 to 2 % on a base of CMF 0.325 and show it errs by up to 10 % when stacked on [R4882]'s CMF term; section 2.4), Lopez and Fortney Eq. 1 and 3 (added an inverse test), Marboeuf ice/rock 1 +/- 0.5 (I extend it to a per-species cumulative sequence), Mulders snow line.
- New here and in G425's "gaps that are yours": low-pressure ice densities and bulk moduli (2.3, via the Journaux splines), triple points recomputed, water inside the snow line (2.2), ocean-world structure (2.4), Thorngren heavy-element masses (2.6), brown dwarf and rogue bodies (2.7, 2.8). Not obtained: [R3932] numerical tables (CDS refused in G425; I did not retry); I did not need them because [R1429] and Lopez and Fortney are already built on SCvH.
- The repo `hydrosphere.py` constants `HP_ICE_DENSITY_KG_M3`, `DEEP_WATER_DENSITY_KG_M3` and the liquidus line are within 1 to 10 % of the recomputed numbers (2.3); the water layer depth formula needs a compressed version above `x_w` about 0.03 (2.4).

---

## 3. How the generator should compute or draw it

### 3.1 Stored quantities

Per body (all new unless noted):

| Name | Unit | Range | Scope | Notes |
|---|---|---|---|---|
| `x_w` (existing `water_mass_fraction`) | mass fraction of body | 1e-6 to 0.5 (0.6 max sourced) | per body | rocky and icy planets and moons |
| `x_ice_other` | mass fraction of body | 0 to 0.3 | per body | CH3OH, NH3, CO2, CO, CH4, H2S, N2 ices summed (needed for comets, icy moons and outer ice worlds); optional, set 0 if not tracked |
| `cmf` | fraction of the dry part | 0 to 0.5 | per body | Part B owns it; used here as the iron share |
| `f_env` | mass fraction of body | 0 to 0.9 | per body | H/He envelope of sub-Neptunes, ice giants |
| `Z_env`, `Y_env` | mass fractions of the envelope | 0 to 0.6; 0.22 to 0.32 | per giant | `Z_env` from `M_Z` below; `Y_env` = initial helium mass fraction of the star's natal gas (star lane) |
| `M_Z` | Mearth | 3 to 1e3 | per giant | total heavy-element mass |
| `core_mass` | Mearth | 0 to 100 | per giant | the part of `M_Z` in a core; rest in the envelope |
| `born_beyond_snow_line` | boolean | | per body | needed to keep the drawn water fraction reproducible if radius is later recomputed |
| `ice_phase_surface` | enum Ih, III, V, VI, VII, X, liquid, none | | per body | derived from `T_surf`, `P_surf`; the module has `hydrosphere` state already |
| `radius` reference level | | | global | one choice (G425 s.3.2 item 4: 1 bar), documented; Lopez and Fortney is 20 mbar; Fortney 1 bar |

Already stored and used: mass, radius, density, `feh`, age, `T_surf`, `P_surf`, orbital distance, snow line from luminosity.

### 3.2 Route 1: keep the drawn bulk figures, fill a composition that fits

Inputs: drawn mass `M`, radius `R` (from the existing class draws). Output: composition consistent with `(M, R)`.

1. **Solid planets and moons (class by existing rock/ice classes; M below about 0.01 Mearth: use density).** Compute the bulk density; subtract compression with the closed form (section 2.4) by writing it for density `rho / rho_dry-Earth-like(M)` and solving for `x`. Algorithm: for given `(M, R, cmf)`, find `x_w` in [0, 1] with `(1.0688 - 0.2317 cmf (1-x)) M^(1/3.7) (1 + 0.4492 x - 0.1597 x^2) = R`. If the root is below 0 (the radius is below the dry model), the planet is iron-richer: lower `cmf` ≥ 0.5 is out of range; report the maximum `cmf` (Fe-only ≈ `R_Fe(M)`, Part B). If the root is above 1 (radius above pure water), the planet needs a gas envelope: go to step 2. The ratio `rho / rho_Earth-like(M)` of 0.65 maps to `x_w` about 0.3 (COMPUTED from the fit; Luque and Palle's gap at 0.65 corresponds to about 30 % water, and their 0.47 mean for the water-world population to about 55 to 60 %).
2. **Sub-Neptunes (R above the water-world curve for x = 0.5).** Fix a core from step 1 at a water fraction drawn from the section 3.3 branch, then solve `f_env` from the Lopez and Fortney inverse (rms 0.106 dex, section 2.5). Use the model radius reference level 20 mbar and `R_atm`.
3. **Ice giants (10 to 60 Mearth, current class I) and gas giants (J) with drawn M and R.** Heavy-element mass from Thorngren's mean relation with its scatter, `log10 M_Z = log10(57.9) + 0.61 log10(M/MJ) + N(0, 0.26)`, capped at `M` (Z at most 1); put the H/He share `1 - Z` in the envelope; the radius stays the drawn one. For M below about 30 Mearth use the ice giant branch: `f_env` log-uniform 0.05 to 0.20 (Helled and Fortney 10 to 20 %; the Thorngren fit is outside its sampled range there and gives 0.4 to 0.6 envelope fraction at 16 Mearth, which disagrees). Consistency check afterwards: `R_model(M, Z)` from the Fortney table (3.4 step 1) must be within the Fortney tabulated scatter (about 5 %); if not, redraw `M_Z` (or keep the drawn radius and flag it).
4. **Brown dwarfs (13 MJ to 0.075 Msun)**: no composition draw beyond `Y`, `Z` of the host gas; radius from the age/mass table (2.7).

Disadvantages: nothing ties `x_w` to formation, so the same planet may sit at an unphysical distance (a 50 % water planet at 0.05 AU around an F star would be a steam world, 2.4); the radius draw must stay within what composition can produce, or step 1 returns no root.

### 3.3 Route 2: draw composition from formation, then derive the bulk figures (recommended)

Recommended because: the bulk density becomes a consequence of composition and mass (what Boss asked), the number of free draws stays small, and every draw has a sourced range. All distribution shapes are parametrisations of the quoted ranges.

**Step A. Snow line.** `r_snow = 2.7 AU * sqrt(L/Lsun) * 10^(N(0, 0.2))` (centre = repo constant; 0.2 dex = [R3193, R3196], covers 2.8 to 4.3 AU for 1 Lsun, section 2.1). Draw it once per star so every planet of a system shares it (store `snow_line_au` per system). Born-beyond flag `B`:
- orbital distance now `a >= r_snow` → `B = true`;
- `a < r_snow` → a migrant with probability `p_mig`; `p_mig = 0` if `M < M_ice,min = 3 (M_star/Msun)^0.48 Mearth` (two NGPPS IV anchors; 2 Mearth for 0.3 to 0.5 Msun per Luque and Palle), else `p_mig` in [0.14, 0.67] (lower end NGPPS IV 0.1 Msun, upper end [R4880]'s 2:1 icy:rocky core ratio fitted to Kepler occurrence). **PARAMETRISATION; weakest link, see gaps.**
- if a giant planet orbits outside the snow line in the same system and formed earlier, set `p_mig` to 0 (blocks the pebbles; Bitsch scenario A) — optional, only if the system already stores the giant.

**Step B. Water mass fraction of the solids `W`.**
- **If `B`** (formed outside): draw `r_ir` (ice-to-rock mass ratio when every ice is condensed) from a normal with mean 1.0, sigma 0.5 truncated to [0.1, 1.5] (Marboeuf); multiply by 0.67 with probability 0.5 if the system is drawn "with refractory organics" (Marboeuf 30 +/- 10 % vs 45 +/- 15 %). Draw the formation temperature `T_form` of the solids as the temperature at the formation distance (store or draw `r_form = r_snow * 10^U(0, 1.3)`, up to ~20 times the snow line), then use the cumulative sequence of 2.2: the set of ices condensed at `T_form` (using the T_cond table at the local pressure; default P_H2 1e-10 bar) and `W = m_H2O / (m_set + 25.0 / r_ir)` with masses in mg per mol of H2: H2O 10.72, CH3OH 2.858, NH3 0.708, CO2 5.237, H2S 0.406, CH4 0.573, CO 3.333, N2 1.165 (25.0 for the total). Result: `W` between 0.14 and 0.39, centre about 0.21 to 0.30. A second parameter `convention_high_water = true` replaces this by `W ~ N(0.5, 0.1)` truncated to [0.3, 0.6] (Zeng, Venturini, Bitsch, Ciesla 50 %). **Boss must choose; section 5.** The remaining ice mass (`1 - m_H2O/m_set` of the ice) goes to `x_ice_other`.
- **If not `B`**: `W = 0` with probability `1 - p_deliv` (dry), else log-uniform `[1.2e-4, 5e-3]` (Raymond and Izidoro, delivered by planetesimals scattered by a forming giant; give `p_deliv = 1` when the system has a giant planet outside the snow line, otherwise a draw of 0.5 — no source for giantless systems; for stars under 0.5 Msun use Ciesla's 0.012 to 0.033 for the innermost planets if the draw says "water-rich planetesimals" — no source gives that probability either).
- Metallicity: no [Fe/H] term on `W` [R4880]. [Fe/H] sets how much solid mass the disk has (existing `disk_surface_density_scale`), which moves mass, not `W`.

**Step C. Body water fraction.** `x_w = W * (1 - f_env)` when the body has an envelope; `x_w = W` otherwise. For rocky planets and icy moons `x_w` replaces the class range draw in `hydrosphere.py` (so `WATER_MASS_FRACTION` classes become the *output* of this rule, not an input); the class label is then chosen from `x_w`, `T_surf` and `P_surf`.

**Step D. Core composition.** `cmf` from Part B (stellar Fe/Si and Mg/Si; G425 s.3.3 gives CMF 0.2 to 0.47 from Thiabaud corners and 0.320 for the Lodders Sun).

**Step E. Radius.** Rocky/icy planets 0.3 to 20 Mearth: the closed form of 2.4 with `(M, cmf, x_w)`, rms 1.2 %, worst 3.6 % against my integration (valid for `T_surf` below the runaway greenhouse; above it add Turbet's steam atmosphere, not implemented). Below 0.3 Mearth the fit extrapolates (the Fe and rock parts get +5 to +8 % too large against Zeng and Sasselov at 0.1 Mearth): use `rho = 1 / sum(w_i / rho_i)` with the Table of 2.3 and a Seager polytrope correction; the error is not measured. Sub-Neptunes with an envelope: `R = R_core(M(1-f_env), cmf, x) + R_env(f_env, M, F, t) + R_atm` with Lopez and Fortney Eq. 3 and 20 mbar; *the water-core + envelope combination is not tested* (G425 and Lopez use rock/iron cores only).

**Step F. Giants.** `M_Z` from Thorngren with 0.26 dex scatter (3.2 step 3), cap `Z = 1`; `core_mass` = `min(M_Z, 10 to 50 Mearth)`? **No source gives a core-mass distribution for modern giants** (the core mass is the least constrained quantity; Thorngren and Helled and Fortney put metals in the envelope, in a diluted core), so store `Z_env = (M_Z - core_mass)/(M - core_mass)` with `core_mass` log-uniform 5 to 20 Mearth (the minimum core for runaway accretion, G425 and [R4880]'s "typically a few tens of Mearth" for transitional planets) as a placeholder. Radius from the [R1429] table (mass, `T_eq` via a bilinear or table interpolation in log M and log T_eq, age in log t, core mass), rms 2 to 4 % and worst 11 to 20 % interpolation error, plus an inflation term for `T_eq > 1000 K` and `M > 0.5 MJ`: heating efficiency `eps(T_eq)` Gaussian peaking at about 2 % near 1500 to 1600 K (falling to 0.2 % at 2500 K; the full Gaussian parameters are in Thorngren and Fortney's text/figure, not extracted, see gaps). Scatter: Chen and Kipping 7.4 % (Jovian), 14.6 % (Neptunian), 4 % (Terran) rms in radius at fixed mass; the project currently uses 10 % and 7 %. Fallback if tables are not wanted: the repo's anchored power laws with Chen and Kipping scatter.

**Step G. Brown dwarfs and the giant mass cut.** Deuterium limit `M_D = 13.55 MJ - 0.9 MJ * log10(Z/Zsun) - 0.2 MJ * (Y - 0.25)/0.01` (linear fits of Spiegel Table 1 around the Z = Zsun, Y = 0.25 centre; Z = Zsun means 0.0189 in their units; 50 % burned; +/-0.8 MJ spread). Hydrogen-burning limit `M_H = 0.074 Msun` at solar Z, Y 0.25, rising toward 0.092 Msun as Z → 0 (Burrows; log-linear in Z is a parametrisation), `-0.1 Msun` per unit Y. Radius: [R0253] COND table interpolated in log age and log mass (solar Z only); for other Z: not available. Cooling luminosity and Teff from Burrows' power laws only for old objects.

**Step H. Rogue planets.** Keep all composition fields; after ejection set irradiation to zero, keep `x_w`, `f_env`, `Z_env`; the existing `rogue_surface` module then computes heat flow and surface state. Preserve the formed branch (outer-branch `W` if the planet was ejected from beyond the snow line; ejection from within is a minority since ejection is driven from the outer planets, a parametrisation).

### 3.4 Consistency invariants a test can assert

1. Mass budget: `x_w + x_ice_other + (1 - x_w - x_ice_other) = 1` and layer masses sum to M: `M_Fe + M_Si + M_H2O + M_ice_other + M_env = M` (relative error 1e-9).
2. `R_model(M, cmf, x_w)` is strictly increasing in `x_w` and decreasing in `cmf` at fixed M (monotone), and `R(x_w = 0, cmf = 0.325, M = 1 Mearth)` within 1.6 % of 1.00 Rearth (a calibration, not a preset).
3. Round trip: from `(M, cmf, x)`, compute R, invert for `x` and recover it within 0.01 (absolute) for x in [0.02, 1] (the fit is monotone).
4. Hydrosphere: `P_base ≈ g x M / (4 pi R^2)` within 3 % for `x_w <= 0.01` against the integrated profile; `x_w * M / (rho_w * 4 pi R^2)` equals the layer depth within 40 % for `x_w <= 0.1`.
5. Giants: `0 <= Z <= 1`, `M_Z <= M`, and `M_Z` median over many draws at 1 MJ within 15 % of 57.9 Mearth.
6. Snow line: median `r_snow / (2.7 sqrt(L))` over many draws is 1.0 +/- 0.05; log scatter 0.2 dex.
7. A body with `born_beyond_snow_line = false` and no giant in the system has `x_w <= 5e-3`.
8. Brown dwarf limits: `M_D(Z = Zsun, Y = 0.25)` = 13.55 MJ within rounding; radius from the table at (age, mass) is within 1 % of the table entry at grid points.

---

## 4. Calibration and test targets (never presets)

| Target | Check | Source |
|---|---|---|
| Earth-like dry planet, CMF 0.325, 1 Mearth | R_model within 1.6 % of 1.00 Rearth | [R4882] (via G425) |
| Pure water, 1 Mearth | 1.401 integrated, fit 1.378 (-1.6 %) vs Zeng 1.41 | [R4882]/2019 |
| Water mixtures, 50 and 75 % water, 0.3 to 10 Mearth | within 3 % of [R4881: Table 1] | READ |
| Luque and Palle water-world population density 0.47 +/- 0.05 of Earth-like, radius 1.97 +/- 0.28 Rearth | 50 % water on an Earth-like rock gives 0.54 in the cold fit (so the observed value needs slightly more water, about 55 to 60 %, or hotter/puffier envelopes) | Luque and Palle |
| Earth's oceans, 2.34e-4 of mass | `P_base` 0.027 GPa (1 Mearth, 2.34e-4: about 0.028 by interpolation of the table) | derived |
| Ice Ih density 916.7 kg/m3 at 273.15 K, 1 bar | 916.72 spline | NIST/Wikipedia SNIPPET and spline |
| Jupiter-mass cool giant (not a preset): 4.5 Gyr, 1 MJ, 1 AU, core 10 to 50 Mearth | Fortney table 1.058 to 0.974 RJ | READ |
| Neptune/Uranus class | H/He envelope fraction 0.10 to 0.20; Z = 0.75 to 0.92 | Helled and Fortney |
| Deuterium limit | 13 +/- 0.8 MJ | [R4170] |
| TRAPPIST-1 (b, c, d) | water at most 2, 0.3, 0.08 % when inside the runaway limit | [R4448] abstract |
| TRAPPIST-1 interior ice | below 1 % (low-iron core), below 7 % (CMF 0.325), below 14 % (CMF 0.5); ice fraction rises with orbital period | [R0039] as quoted by [R0602] (2ND-HAND) |

---

## 5. Gaps: Boss must research, and to verify

### 5.1 Boss must research (ranked)

1. **Water fraction of solids beyond the snow line**: the two conventions (0.14 to 0.39 Marboeuf-derived, vs 0.5 Lodders/Zeng/Venturini). Needs [R2726] full text (ApJ 591:1220, Table 8 condensation temperatures at 1e-4 bar and the 50 % water claim; IOP PDF behind a bot check from here) and [R2865: Figure 8 and 9] values by distance; this controls the whole ice-world branch. Number it gives: water mass share of condensed solids at the H2O line and at larger radii.
2. **Probability that a planet inside the snow line is a migrant** (`p_mig`) as a function of mass and stellar mass: [R4528]/2022 and Emsenhuber (NGPPS) give planet-by-planet data (Zenodo / DACE tables) not extracted here. Gives the fraction of 2 to 20 Mearth close-in planets that are water-rich (0.14 to 0.67 is a range, not a distribution).
3. **Hot-Jupiter inflation efficiency function**: the Gaussian hyper-parameter posterior in [R4351] (Table of medians; text read gives only peak about 2 % at 1500 to 1600 K, 0.2 % at 2500 K).
4. **Core mass and diluted core of giants**: modern interior models ([R4586] [R1038, R1896] review "Giant planet interiors") give core mass and metal gradient; not read; core mass distribution is a placeholder.
5. **Solid CH4, NH3, H2S, CH3OH ice densities and bulk moduli** ([R1445] review, ScienceDirect, refused; [R2392] tables). Gives the low-pressure density of small icy bodies and the outer-Solar-System-type composition rows.
6. **[R3932] numerical tables ** (CDS refused in G425): only needed if a lane integrates the H/He envelope itself; the Fortney tables replace them.
7. **[R1457, R1458] superionic phase diagram** (PRB, APS, not on arXiv by title search); the superionic boundary in P-T and the density along it; for deep layers of water worlds and ice giants.
8. **[R1754] AQUA** water EOS tables (title search found nothing on arXiv; Zenodo is reachable but I did not find the record); gives a single tabulated water EOS from 0.1 bar to 400 TPa and 300 to 5000 K, which would replace my 300 K isotherm plus polytrope.
9. **Ejection and stripping of volatiles from rogue planets**: no source found on how ejection changes water or envelope content.
10. **[R4352] Z/Z_star scatter and [R3077]** mass-metallicity details (fetched but not read).

### 5.2 To verify (RECALLED or weak)

1. Solid CH4 0.52 and NH3 0.82 to 0.85 g/cm3 crystallographic estimates (RECALLED lattice constants).
2. CH4 clathrate density about 0.9 g/cm3 (SNIPPET).
3. [R3935] snow-line fit as quoted in [R0405] (2ND-HAND).
4. [R0039] TRAPPIST-1 ice limits (2ND-HAND via [R0602]).
5. The SeaFreeze splines are [R2204]; I used them as data and checked them against NIST-class values only through Ih (916.7, 0.934 at 93 K SNIPPET) and the triple points (agree with Zeng and Sasselov to 1 %); the original paper and its stated uncertainties were not read.
6. Units of Thorngren's `M_Z` fit (Earth masses for M in Jupiter masses) are inferred from Fig. 7 and the abstract text; the printed formula has no explicit units.
7. The transition formula `M_ice,min = 3 (M_star/Msun)^0.48` is a two-point fit to NGPPS IV anchors (about 1 and about 3 Mearth), read off figure descriptions, not a fitted law.

---

## 6. Proposed work items for the TODO tree

**Research and design**
1. *Water-fraction convention decision*; settle 0.14 to 0.39 vs 0.5 for solids beyond the line, with Lodders full text; depends on: Boss; size small.
2. *Migrant probability study*; extract planet-by-planet ice fractions from NGPPS/Venturini tables to replace `p_mig`; depends on: 1; size medium.
3. *Giant interior design note*; core mass, diluted core, `Z_env` placeholder replaced by modern interior models; depends on: none; size medium.
4. *Design note: volatile-rich and gas bodies* (this part into `docs/design/`); depends on: Part B and E done; size small.

**Tables**
5. *Ice and water materials table*: Journaux spline outputs (rho, K, K', alpha per phase) as a compact table plus the VII/X curve to 1 TPa; also a 300 K isotherm for the integrator; depends on: none; size small.
6. *Volatile ices table*: T_cond(P) for the eight species (Marboeuf curves), molar masses, abundance shares, densities (CO2, N2, CO measured; CH4, NH3, H2S pending); depends on: 5; size small.
7. *[R1429] giant radius table + [R0253] brown dwarf table* as data files with interpolators; depends on: none; size small.

**Drawing rules**
8. *Snow line draw per system* with 0.2 dex scatter and `born_beyond_snow_line` flag; depends on: 1; size small.
9. *Water fraction draw* (route 2 steps A to C), replacing class ranges in `hydrosphere.py`; depends on: 1, 6, 8; size medium.
10. *Giant heavy-element draw* (`M_Z`, `Z_env`, `Y_env`, `core_mass`); depends on: 3; size small.
11. *Deuterium and hydrogen burning limits with Z and Y*; depends on: none; size small.

**Derivations**
12. *Closed-form radius for rock, iron and water* (the 1.0688 fit), incl. inverse; depends on: B's CMF; size small.
13. *Sub-Neptune envelope* (Lopez inverse and forward, 20 mbar convention); depends on: 12; size small.
14. *Giant and brown dwarf radius lookup* with inflation; depends on: 7, 10; size medium.
15. *Hydrosphere upgrade*: compressed water-layer depth, `P_base`, constants (`HP_ICE_DENSITY`, `DEEP_WATER_DENSITY`), liquidus from recomputed boundaries; depends on: 5, 12; size medium.
16. *Steam-atmosphere radius for irradiated water worlds* ([R4448] formula); depends on: 12; size medium.
17. *Rogue planet composition carry-over rules*; depends on: 9, 10; size small.

**Storage**
18. *Schema*: `x_ice_other`, `f_env`, `Z_env`, `Y_env`, `M_Z`, `core_mass`, `born_beyond_snow_line`, per-system `snow_line_au`; depends on: 9, 10; size medium.

**Tests**
19. *Invariant tests* of 3.4; depends on: 12 to 15; size medium.
20. *Calibration tests* of section 4 (as tolerances, no presets); depends on: 12, 14; size small.

**Docs**
21. *Docs and sources section* for this part, including the correction that [R1319] is not [R4352] (it is 1511.07854); depends on: 4; size small.

---

## 7. Verification summary

Counts are mentions of each label word in this file (grep, 2026-10-17; they include the legend and gap lists, so they are upper bounds on rows):

| Label | Mentions |
|---|---|
| READ | 78 |
| 2ND-HAND | 8 |
| COMPUTED | 33 |
| SNIPPET | 7 |
| RECALLED | 7 |
| NOT REACHED | 5 |

Caveats: the 4.5 Gyr Fortney 1-bar radii to use are those in the correction paragraph below the 2.6 table. The Earth-ocean base pressure in section 4 is interpolated from the ocean table (COMPUTED), not a measured target. The Fortney row's first draft contained a self-correction; it has been removed.

READ items rest on 20 full-text arXiv papers ([R3992] [via G425 and own read of Tables 1 to 4], [R0253, R0405, R0602, R0606, R0758, R1429, R1892, R2865, R4170, R4352, R4528, R4880, R4881], [R1268, R1269], [R0803, R2774, R3681], Ma 2016, [R4448] [abstract], [R4351] [abstract and text], [R4663] [abstract], [R0430] [density lines]) plus PubChem, Nature (abstract) and IOP (abstract) pages. The COMPUTED items are four own calculations: the Journaux-spline ice properties and phase boundaries, the hydrostatic reference and its validation against Zeng and Sasselov Table 1, the 175-point three-layer fit and its errors, and the condensation and cumulative-ice computations, plus the Lopez inverse test and Fortney interpolation tests.

Three corrections to the brief's pointers: [R4352] is [R4352] (1511.07442 is an unrelated spintronics paper); [R2865] is [R2865] (0804.0406 is the 2008 predecessor); Bitsch's water-content paper is 2021 [R0405].
