# Part E: surface ground, radiogenic heat, radiation shielding, material properties

Draft of 2026-10-17 for planetGen (`dwhagar/planetGen`). Agent E of five research agents. Read-only research: nothing in the repo, `docs/` or `TODO.md` was touched. Sources were fetched after the network widening noted in the brief (12:03 AM EDT). Downloads sit in their own directories and all Python on them ran with `-I`.

Label key (from the brief): READ = paper or data table read and quoted; 2ND-HAND = read in a reachable paper or open-source package that cites the original; COMPUTED = re-derived here with Python from the stated inputs; SNIPPET = search summary only; RECALLED = memory, not checked; NOT REACHED = source exists, nothing reachable gave the number. A repo file read for what the code does is labelled READ and the source column says "repo". Every table row below ends with its label.

Rule 1 of the brief (nothing fixed, nothing copied from the Solar System): every Solar System number below is a calibration or test target or a worked example. None is a preset. Where the literature gives only a range, the range and a distribution shape are given and called a parametrisation.

## 1. Scope and interfaces

This part covers four things.

1. Ground: which crust and surface types a body can have, their oxide chemistry, regolith structure, thermal properties, albedo in 0.4 to 0.7 micrometres, space weathering, and melting, latent heat and Debye data for rock, iron, water ice, CO2, CH4, NH3, N2, SO2 and S.
2. Radiogenic heat: K, Th and U in chondritic, bulk silicate and crustal reservoirs, isotope heat with time, early 26Al and 60Fe heating, a heat-flow formula with the Urey ratio, and the link to the star's abundances.
3. Radiation at and under the ground: gamma and radon from the rock, shielding by regolith, basalt, ice, iron, CO2 and N2, the atmosphere-plus-ground column rule, magnetic field and stellar activity, and what `physics/radiation.py` already does.
4. One machine-loadable materials table (section 2.10), with the ten atmosphere gases as pointer rows.

| Interface | Direction | Quantity | Unit | Status |
|---|---|---|---|---|
| Part B (rocky structure) | in | core mass fraction, mantle oxide fractions (at least Mg, Si, Fe, Al, Ca, Na, K), mantle FeO | mass fraction | new (B to store) |
| Parts C and D (ices, small bodies, moons) | in | surface ice species and mass fractions, S and C inventory, tidal heating | mass fraction, W/m2 | new |
| Part D | out | shared heat-per-kg formula, decay correction, and the materials table | n/a | this part; D keeps asteroid taxon albedo and chondrite chemistry |
| Part A (stars) | in | [Th/H] or [Eu/Mg] proxy, [Mg/H], [Al/H], [Fe/H], age | dex, Gyr | [Fe/H] and age stored; Th, Al, Mg need A |
| Atmosphere lane | in | surface pressure, gravity, gas column, O2 partial pressure | Pa, g, g/cm2, kPa | stored (`atmospheric_pressure`, `gravity`) |
| Atmosphere lane | out | ground type, band albedo, outgassing-relevant K, ground gamma | various | this part |
| Hydrosphere and magnetism code | in | `water_mass_fraction`, `mobile_lid_chance`, magnetic moment | various | stored |

## 2. Findings

### 2.1 What the repo already computes

| Item | What the code does | Source | Label |
|---|---|---|---|
| Cosmic-ray (GCR) dose | Log-interpolation in column X between `(0, 500)`, `(16.4, 237)`, `(100, 110)` mSv/yr, then `exp(-(X-100)/164)` with a floor of 0.1 percent of the sea-level value; column from `X = P/g` in g/cm2 | repo `physics/radiation.py` lines 64 to 118, 273 to 277 | READ |
| Magnetic cut | Dipole cutoff `14.9 GV * (m/7.768e22) * (R_E/R)^3`; dose cut up to 40 percent, full at 15 GV | repo `radiation.py` lines 71 to 74, 127 to 139 | READ |
| Heliosphere multiplier | `1 + 1.5 exp(-X/100) * compression` (1 to 2.5) | repo `radiation.py` lines 76 to 77, 150 to 152 | READ |
| Stellar particle (SEP) dose | Free dose `20 mSv/yr * sqrt(flare rate / 10^-2.5) / d^2`; shielded by `exp(-p_cut / 0.2 GV)` where `p_cut` comes from the range law `R = 0.0022 E^1.77` g/cm2 (E in MeV) | repo `radiation.py` lines 79 to 85, 155 to 176 | READ |
| Ground dose | `a_rad * (0.48 * heat_ratio + 0.5 if radon)` mSv/yr; `heat_ratio` is the rogue radiogenic flux at the body's age over that at 4.5 Gyr; radon only above 1000 Pa | repo `radiation.py` lines 87 to 90, 179 to 184 | READ |
| `a_rad` | One log-uniform draw in (0.5, 2.0) per body in `generate_dose`, independent of the heat calculation | repo `radiation.py` (docstring of `generate_dose`), `tuning.py` `ROGUE_RADIOGENIC_ABUNDANCE_RANGE` | READ |
| Radiogenic heat flux | Four isotopes at [R4449] mantle values (U 31 ppb, Th 124 ppb, K 310 ppm), scaled back in time by `exp(lambda (4.5 - t))`, mantle fraction 0.67, abundance 1.0 | repo `rogue_surface.py` (the former `radiogenic_flux_w_m2` is gone: `radiogenic` and `interior_heat` give the heat); `tuning.py` `ROGUE_RADIOGENIC_ISOTOPES` | READ |
| Internal heat flux | Radiogenic flux at the star's age divided by Urey ratio 0.67 | repo `planets.py` `internal_heat_flux_w_m2`; `tuning.py` `ROGUE_UREY_RATIO` | READ |
| Ice | density 917 kg/m3, vacuum ice stable below 220 K, conductivity 567/T | repo `hydrosphere.py` | READ |
| Mobile lid | chance 0.6 (mass up to 5 Earth masses) or 0.3, only with a surface ocean and mass of 0.5 Earth masses or more | repo `magnetism.py` `mobile_lid_chance` | READ |

Problems found in the existing code (all COMPUTED from the tables below unless noted):

1. Two separate abundance paths. Heat flux uses abundance 1 at the star's age; the dose draws its own `a_rad`. A body cannot be rich in Th for its dose and average for its heat. Fix: one stored set of K, Th, U feeds both.
2. The ground dose is scaled by total radiogenic heat ratio, but gamma comes from specific isotopes at the surface rock, not from the mantle heat. Heat is dominated by U-238 and Th-232 today and by K-40 early (section 2.8); gamma weights differ.
3. No crust type. Surface gamma differs by a factor of about 650 between an ultramafic and a felsic surface (table 2.9a), against the code's factor of 4 range.
4. Star age is used as planet age (the repo comments say so for the heat call); fine to first order, but a rogue planet's age is drawn separately.
5. The code's mantle concentrations give 29.8 TW for a 4.03e24 kg bulk silicate Earth against 19.8 to 21.8 TW from current bulk-silicate models (table 2.8b): 37 to 50 percent too high.
6. Urey ratio 0.67 against about 0.47 to 0.52 from the totals in table 2.8e.
7. `ROGUE_RADIOGENIC_ABUNDANCE_RANGE` (0.5, 2.0) is narrower than the star-to-star scatter in the literature (0.3 to 3, table 2.8d).

### 2.2 Crust types and how they form

Source for all rows: Hu, [R2024: Table 2] (mineral assemblages) and Table 4 (near-infrared geometric albedo of an airless planet covered by one type, G star at 1 AU). Both READ. [R1719], section 4 (READ) gives the melting rule.

| Type | Assemblage (Hu Table 2) | Where it comes from | J / H / K band geometric albedo (Hu Table 4) | Label |
|---|---|---|---|---|
| Ultramafic (primary, mantle overturn) | 60% olivine, 40% enstatite | magma-ocean crust after mantle overturn, or hot early lavas | 0.32 / 0.55 / 0.61 | READ |
| Feldspathic or anorthositic (primary, no overturn) | 97% Fe-plagioclase, 3% augite | magma-ocean flotation crust | 0.56 / 0.61 / 0.63 | READ |
| Basaltic (secondary) | 76% plagioclase, 8% augite, 6% enstatite, 5% glass, 1% olivine | volcanism, decompression melting | 0.24 / 0.25 / 0.24 | READ |
| Granitoid (tertiary) | 40% K-feldspar, 35% quartz, 20% plagioclase, 5% biotite | tectonic reprocessing with water | 0.60 / 0.65 / 0.71 | READ |
| Clay (aqueous alteration) | 50% Mg-smectite, 50% serpentine | liquid water on basalt | 0.68 / 0.73 / 0.51 | READ |
| Ice-rich silicate | 50% water ice, 50% basalt | ice-rich mantle exposed | 0.34 / 0.23 / 0.21 | READ |
| Fe-oxidised | 50% nanophase hematite, 50% basalt | oxidative weathering | 0.21 / 0.21 / 0.21 | READ |
| Metal-rich | pyrite | mantle stripped | 0.12 / 0.18 / 0.19 | READ |

Melting facts for the genealogy ([R1719], section 4, READ): decompression melting works in any convecting mantle that is not buried under a thick volatile envelope, so a secondary basaltic crust is the default for any rocky body that is still or was once convecting; the mantle solidus moves by more than 100 degrees C with alkali content; thermodynamic melt models are most accurate for Earth-like compositions and are untested for strange ones. Crust is also the second silicate reservoir: it holds most of the incompatible heat producers K, Th, U (same section, READ).

Types the brief lists that Hu does not: komatiitic (high-Mg lava), carbonaceous, sulfur, evaporite, carbon-rich. Their oxide or stoichiometric tables are in 2.3 and the albedo in 2.6; komatiite oxide numbers and a carbon-rich crust rule were not reached (section 5).

### 2.3 Oxide tables

Source: the reference-composition data files shipped in the `pyrolite` package (version 0.3.8.post1 wheel from PyPI, files `refcomp/*.csv`). Those files carry the numbers of [R1494] (N-MORB), [R3842] (continental crust), [R4775] (depleted mantle) and [R3417] (bulk silicate Earth and CI). The originals sit behind Elsevier and AGU and were not reached, so each number is 2ND-HAND. Oxide percentages for the bulk silicate Earth and CI were computed from element weight percent with `periodictable` masses (COMPUTED). FeO is total iron as FeO. These are the compositions of Earth reservoirs used as chemistry models for crust types; nothing here is a preset (a draw shifts them, see section 3).

| Reservoir (ground type it models) | SiO2 | TiO2 | Al2O3 | FeO | MgO | CaO | Na2O | K2O | Source | Label |
|---|---|---|---|---|---|---|---|---|---|---|
| N-MORB (basaltic secondary crust) | 50.42 | 1.53 | 15.13 | 9.81 | 7.76 | 11.35 | 2.83 | 0.14 | [R1494] via pyrolite `NMORB_Gale2013` | 2ND-HAND |
| Upper continental crust (felsic tertiary) | 66.6 | 0.64 | 15.4 | 5.04 | 2.48 | 3.59 | 3.27 | 2.8 | [R3842] via pyrolite `UCC_RudnickGao2014` | 2ND-HAND |
| Middle continental crust | 63.5 | 0.69 | 15.0 | 6.02 | 3.59 | 5.25 | 3.39 | 2.3 | same, `MCC_` | 2ND-HAND |
| Lower continental crust (mafic granulite) | 53.4 | 0.82 | 16.9 | 8.57 | 7.24 | 9.59 | 2.65 | 0.61 | same, `LCC_` | 2ND-HAND |
| Bulk continental crust | 60.6 | 0.7 | 15.9 | 6.7 | 4.7 | 6.4 | 3.1 | 1.8 | same, `BCC_` | 2ND-HAND |
| Depleted mantle (ultramafic primary crust and mantle) | 44.71 | 0.13 | 3.98 | 8.18 | 38.73 | 3.17 | 0.13 | 0.01 | [R4775] via pyrolite `DMM_WorkmanHart2005` | 2ND-HAND |
| Bulk silicate Earth (undifferentiated rocky mantle plus crust) | 45.4 | 0.21 | 4.5 | 8.1 | 36.76 | 3.65 | 0.35 | 0.03 | [R3417] via pyrolite `PM_`, oxides computed from elements | COMPUTED |
| CI chondrite (carbonaceous; oxide equivalents of the elements, S and C not included) | 22.89 | 0.07 | 1.59 | 24.01 | 15.82 | 1.27 | 0.67 | 0.07 | [R3417] via pyrolite `CH_`, oxides computed from elements | COMPUTED |
| Anorthosite end-member (anorthite CaAl2Si2O8, stoichiometry) | 43.2 | 0 | 36.6 | 0 | 0 | 20.2 | 0 | 0 | stoichiometry with `periodictable` masses | COMPUTED |

Not found in reachable sources (gap list, section 5): komatiite oxide table (RECALLED only: MgO about 18 to 32 weight percent, source needed), Mars, Venus and lunar crust oxide tables (the Wikipedia pages fetched had no table; the primary sources are Taylor and McLennan and Wanke style compilations behind publishers). Chemistry of the other types: ice is H2O (formula, COMPUTED mass fractions H 0.1119, O 0.8881), sulfur is S, evaporites are halite NaCl, gypsum CaSO4.2H2O and calcite CaCO3 (formulas), carbonaceous is the CI row above and agent D's [R2729] table (READ in part D, section 2.5).

### 2.4 Regolith structure and thermal properties

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| Surface-layer density | 1100 | kg/m3 | [R1861], Appendix table (lunar fit) | READ |
| Deep-layer density | 1800 | kg/m3 | same | READ |
| Density with depth | `rho(z) = rho_d - (rho_d - rho_s) exp(-z/H)` | n/a | same, equation 1 | READ |
| Scale depth H | 0 to above 0.2, average 0.06 | m | same | READ |
| Phonon conductivity, surface and deep | 7.4e-4 and 3.4e-3 | W/m/K | same | READ |
| Conductivity law | `K = Kc [1 + chi (T/350)^3]`, chi = 2.7; Kc linear in density between the two layers | n/a | same | READ |
| Heat capacity polynomial | `cp = -3.6125 + 2.7431 T + 2.3616e-3 T^2 - 1.2340e-5 T^3 + 8.9093e-9 T^4` | J/kg/K | same) | READ |
| cp from that polynomial | 138 (50 K), 254 (90 K), 424 (150 K), 555 (200 K), 720 (273 K), 771 (300 K), 850 (350 K), 910 (400 K) | J/kg/K | evaluated from the polynomial | COMPUTED |
| Thermal inertia at 273 K, surface / deep | 36.6 / 100 | J/m2/K/s^1/2 | evaluated from the Hayne laws, `I = sqrt(K rho cp)` | COMPUTED |
| Granular conductivity model | `lambda = lambda_solid H(r,T,phi) + 8 sigma eps T^3 Lambda`; `Lambda = e1 (1-phi)/phi r`, e1 = 1.34; `H = [(9 pi/4)(1-mu^2) gamma(T)/(E r)]^(1/3) f1 exp(f2 phi) chi`; f1 = 5.18e-2, f2 = 5.26, chi = 0.41; basalt grain values mu = 0.25, E = 7.8e10 Pa, gamma = 6.67e-5 T J/m2; lambda_solid 1.5 to 4 W/m/K | n/a | [R1349: Eq. 6 to 8] | 2ND-HAND |
| Example from that model, filling factor 0.3 / 0.5, T = 250 K, lambda_solid = 3 | r = 10 um: 1.8e-3 / 4.7e-3; 100 um: 2.7e-3 / 3.0e-3; 1 mm: 2.0e-2 / 9.5e-3; 1 cm: 0.20 / 0.086 | W/m/K | my evaluation; implementation not yet checked against the original paper | COMPUTED |
| Asteroid thermal inertia range used for grain-size inversion | 150 to 300 | J/m2/K/s^1/2 | [R4577], section on thermal inertia and albedo | READ |
| Water ice conductivity | `k = 567/T` (crystalline) | W/m/K | [R2392] as quoted in [R1055] | 2ND-HAND |
| Water ice heat capacity | `cp = 7.49 T + 90` | J/kg/K | same | 2ND-HAND |
| Limit of the ice laws | good above about 100 K; amorphous ice conducts far less | n/a | search summary of the same paper | SNIPPET |
| Lunar interior heat flow in that model | 0.018 | W/m2 | [R1861] Appendix (cites [R2560]); calibration target only | READ |

Failure modes: the Hayne laws are fitted to lunar basaltic and anorthositic regolith in vacuum and carry no gas conduction; with an atmosphere the pore-gas term dominates for coarse grains and must be added by the atmosphere lane. The Gundlach and Blum form needs the grain radius and filling factor, which no part yet stores.

### 2.5 Heat capacity, Debye data, melting, latent heat

Phase data for pure materials are in the CSV (section 2.10). Mineral grain data come from the [R4232] database parameters read from the `burnman` 2.1.0 wheel source (`SLB_2011.py`), so they are 2ND-HAND. Debye heat capacity per kg was computed from those parameters (the Debye function with the Debye temperature and atom count).

| Mineral | Molar mass g/mol | Grain density g/cm3 | Debye temperature K | Debye Cv at 150 / 300 / 600 / 1000 K (J/kg/K) | Dulong-Petit limit J/kg/K | Label |
|---|---|---|---|---|---|---|
| Forsterite | 140.691 | 3.2266 | 809.2 | 402 / 885 / 1135 / 1201 | 1241 | COMPUTED |
| Fayalite | 203.771 | 4.4021 | 618.7 | 415 / 699 / 813 / 841 | 857 | COMPUTED |
| Enstatite | 200.774 | 3.2034 | 812.2 | 400 / 884 / 1136 / 1202 | 1242 | COMPUTED |
| Ferrosilite | 263.854 | 4.0014 | 674.5 | 409 / 744 / 888 / 924 | 945 | COMPUTED |
| Diopside | 216.547 | 3.2791 | 781.6 | 397 / 839 / 1060 / 1117 | 1152 | COMPUTED |
| Anorthite | 278.203 | 2.7652 | 752.4 | 427 / 868 / 1079 / 1133 | 1166 | COMPUTED |
| Albite | 262.218 | 2.6104 | 713.8 | 492 / 947 / 1153 / 1206 | 1237 | COMPUTED |
| Quartz | 60.083 | 2.5384 | 816.3 | 398 / 883 / 1137 / 1205 | 1245 | COMPUTED |

Debye data for water ice, CO2, CH4, NH3, N2, SO2, S and iron were not reached (section 5). The Debye model gives Cv from one Debye temperature per mineral. Forsterite Cv at 300 K is 885 J/kg/K from the table; its measured Cp is roughly 840 J/kg/K (that calorimetric value is RECALLED), so a single Debye temperature overestimates by about 5 percent there, and Cp is a further few percent above Cv at high temperature. For regolith use the Hayne polynomial above (READ) between 50 and 400 K.

Latent heat and melting in the CSV: pure-material values are from NIST WebBook pages (READ) and from the CRC Handbook tables shipped in the `chemicals` 1.5.2 package (2ND-HAND). Rock and iron melting and vaporisation values are mostly RECALLED and flagged in the CSV labels.

### 2.6 Albedo versus wavelength, 0.4 to 0.7 micrometres

All band values are means of laboratory reflectance spectra, averaged over the band by me (COMPUTED from READ data). Two open datasets:

- USGS Spectral Library version 7 (splib07a, sciencebase.gov, USGS Data Series 1035). Beckman (BECK) spectra are absolute reflectance of powders or slabs; ASD spectra are field or lab reflectance. File names are given so the spectrum can be found.
- [R3425] (ApJ 981, 130; data at Zenodo record 17986129, `hemi_refls.csv`): hemispherical reflectance of slab, crushed and powder samples.

Caveat for every row: a lab reflectance factor is neither a geometric albedo nor a Bond albedo. Use it as a normal reflectance proxy. Real regolith is darkened by space weathering (2.7). The relation to Bond albedo is a gap.

| Material (ground type) | State and file | 0.4-0.5 | 0.5-0.6 | 0.6-0.7 | 0.4-0.7 mean | Label |
|---|---|---|---|---|---|---|
| Basalt, powder (basaltic) | Paragas `K1919_basalt_powder` | 0.207 | 0.214 | 0.224 | 0.213 | COMPUTED |
| Basalt, slab | Paragas `new_basalt_slab` | 0.125 | 0.121 | 0.116 | 0.122 | COMPUTED |
| Basalt, crushed | Paragas `new_basalt_crushed` | 0.094 | 0.099 | 0.102 | 0.097 | COMPUTED |
| Basaltic andesite, powder | Paragas | 0.149 | 0.153 | 0.166 | 0.154 | COMPUTED |
| Fresh basalt (dark) | USGS BR93-46B (BECK) | 0.077 | 0.092 | 0.107 | 0.092 | COMPUTED |
| Weathered basalt | USGS BR93-43 (BECK) | 0.105 | 0.128 | 0.145 | 0.126 | COMPUTED |
| Dunite xenolith, powder (ultramafic) | Paragas | 0.629 | 0.690 | 0.706 | 0.665 | COMPUTED |
| Olivine pyroxenite, powder (ultramafic) | Paragas `EG_19_70` | 0.348 | 0.371 | 0.370 | 0.360 | COMPUTED |
| Olivine clinopyroxenite, powder | Paragas `EG_19_63` | 0.393 | 0.441 | 0.435 | 0.417 | COMPUTED |
| Olivine Fo89 below 60 um | USGS GDS70.d (BECK) | 0.466 | 0.571 | 0.589 | 0.544 | COMPUTED |
| Olivine Fo89 165 um | USGS GDS70.a (BECK) | 0.235 | 0.384 | 0.362 | 0.330 | COMPUTED |
| Granite, powder (felsic) | Paragas Orlando Gold | 0.688 | 0.760 | 0.815 | 0.739 | COMPUTED |
| Granite, powder | Paragas Dalmatian | 0.395 | 0.436 | 0.451 | 0.420 | COMPUTED |
| Granite, slab | Paragas Orlando Gold | 0.269 | 0.314 | 0.366 | 0.305 | COMPUTED |
| Andesite, powder | Paragas STM 101 | 0.247 | 0.280 | 0.297 | 0.268 | COMPUTED |
| Anorthite, below 74 um (anorthositic, pure) | USGS GDS28 (BECK) | 0.710 | 0.734 | 0.750 | 0.732 | COMPUTED |
| Anorthite, mineral specimen | USGS HS349.3B (BECK) | 0.595 | 0.653 | 0.663 | 0.638 | COMPUTED |
| Kaolinite (clay) | USGS CM9 (BECK) | 0.687 | 0.751 | 0.768 | 0.737 | COMPUTED |
| Serpentine (clay) | USGS HS8.3B (BECK) | 0.250 | 0.267 | 0.241 | 0.253 | COMPUTED |
| Hematite 60-104 um (oxidised) | USGS GDS69.c (BECK) | 0.092 | 0.086 | 0.085 | 0.087 | COMPUTED |
| Hematite below 10 um | USGS GDS69.g (BECK) | 0.041 | 0.042 | 0.061 | 0.048 | COMPUTED |
| Goethite (oxidised, red) | USGS WS222 (BECK) | 0.030 | 0.122 | 0.206 | 0.122 | COMPUTED |
| Ilmenite (metal-bearing) | USGS HS231.3B (BECK) | 0.055 | 0.054 | 0.056 | 0.055 | COMPUTED |
| Magnetite | USGS HS195.3B (BECK) | 0.059 | 0.058 | 0.058 | 0.058 | COMPUTED |
| Desert varnish | USGS GDS141 (BECK) | 0.046 | 0.061 | 0.084 | 0.064 | COMPUTED |
| Carbon black (carbon-rich, dark end) | USGS GDS68 (ASD) | 0.019 | 0.019 | 0.019 | 0.019 | COMPUTED |
| Kerogen (carbon-rich organic) | USGS BK-Cornell (BECK) | 0.013 | 0.016 | 0.021 | 0.016 | COMPUTED |
| Sulfur powder (yellow, absorbs blue) | USGS GDS94 (BECK) | 0.362 | 0.777 | 0.824 | 0.664 | COMPUTED |
| Halite (evaporite) | USGS HS433.3B (BECK) | 0.817 | 0.855 | 0.876 | 0.850 | COMPUTED |
| Gypsum (evaporite) | USGS HS333.3B (BECK) | 0.901 | 0.917 | 0.925 | 0.915 | COMPUTED |
| Calcite (evaporite or carbonate) | USGS WS272 (BECK) | 0.911 | 0.941 | 0.946 | 0.934 | COMPUTED |
| Melting snow (fine, clean) | USGS mSnw01a (ASD) | 0.832 | 0.832 | 0.822 | 0.829 | COMPUTED |
| Melting snow (older) | USGS mSnw05 (ASD) | 0.622 | 0.627 | 0.622 | 0.624 | COMPUTED |
| Slush | USGS mSnw16 (ASD) | 0.194 | 0.212 | 0.204 | 0.203 | COMPUTED |
| Hexane frost at 85 K (cold frost stand-in) | USGS (ASD) | 0.746 | 0.751 | 0.755 | 0.751 | COMPUTED |
| Water ice at 77 K | USGS GDS136 (BECK) starts at 0.86 um, so no visible value | n/a | n/a | n/a | n/a | NOT REACHED |
| Red organic-rich ice, Arrokoth surface (calibration target) | [R1970]: geometric albedo 0.21 (+0.05, -0.04) at 550 nm, 0.24 at 610 nm | n/a | 0.21 | 0.24 | n/a | READ |
| Near-infrared types (J 1.1-1.4, H 1.5-1.8, K 2.0-2.4 um) | see table 2.2 | n/a | n/a | n/a | n/a | READ |

Not reached: visible albedo of CO2, CH4, NH3, N2 and SO2 frosts; Mars bright dust (search snippets only: Pancam single-scattering albedo of outcrop about 0.9 and spherule soil about 0.6 at 753 nm, SNIPPET); tholin-like organics other than Arrokoth; pyrite and metal in the visible (Hu gives only the near infrared, table 2.2).

Grain-state effect (COMPUTED from the Paragas columns above, powder over slab, mean 0.4 to 0.7): Orlando granite 2.4, basaltic andesite 1.9, Dalmatian granite 2.9, STM 101 andesite 2.5, K1919 basalt 2.1, olivine pyroxenite 2.8 (powder over crushed 4.1), dunite 3.4 (powder over crushed). Powder is 2 to 4 times brighter than slab or crushed rock for most igneous types. The generator must therefore pick a grain state, not just a rock type.

### 2.7 Space weathering

| Item | Value | Source | Label |
|---|---|---|---|
| Model form | the unweathered fraction decays exponentially with time since the surface was refreshed (single-tau and dual-tau models) | [R4713] | READ |
| Dual-tau characteristic times, S-complex asteroids | weathering 960 +/- 160 Myr, gardening 2000 +/- 290 Myr | same, quoting the earlier fit | READ |
| Enhanced dual-tau fit | weathering 2050 +/- 80 Myr, gardening 4400 (+700) Myr | same | READ |
| Spread of published timescales | 1 Myr to 1 Gyr for significant reddening | Marchi et al. as summarised in a search result | SNIPPET |
| Direction of the effect | darkening and reddening of the spectral slope, from nanophase iron | search summary | SNIPPET |
| Lab-to-regolith gap | lunar mean Bond albedo in the Hayne model is 0.12 [R1861] against 0.15 to 0.74 for the lab powder means of 2.6 (basaltic andesite to felsic): weathered regolith is darker than fresh lab powder by a factor of about 1.3 to 6, depending on type. A Bond albedo and a lab reflectance factor are not the same quantity, so this is an order-of-magnitude guide only | [R1861] and table 2.6 | COMPUTED |

Parametrisation (mine, not from a paper): `A_band(t) = A_fresh * [1 - d * (1 - exp(-t_surface / tau_w))]` with `tau_w` drawn log-uniform in 1 to 2 Gyr and `d` drawn uniform in 0.3 to 0.7 for airless silicate regolith. The values of `d` and the exponent form are placeholders to fit against lunar and asteroid reflectance data (section 5).

### 2.8 Radiogenic heat

#### 2.8a Isotope constants and heat per kilogram

Method ([R1203], equation 6, READ): heat per kg of isotope `h = N_A * lambda * Qh / mu`, with `Qh` the heat released per decay after removing the neutrino energy. `Qh` and `mu` are from Dye's Table 2; half-lives are from the `radioactivedecay` 0.6.1 package (ICRP-107 data); my values reproduce Dye's table within 0.1 percent for U-238, U-235 and Th-232 and within 2.2 percent for K-40.

| Isotope | Half-life Gyr | Qh pJ (Dye Table 2) | h W/kg of isotope (Dye) | h W/kg (this work) | Label |
|---|---|---|---|---|---|
| U-238 | 4.468 | 7.648 | 9.513e-5 | 9.511e-5 | COMPUTED |
| U-235 | 0.704 | 7.108 | 5.6847e-4 | 5.682e-4 | COMPUTED |
| Th-232 | 14.05 | 6.475 | 2.628e-5 | 2.627e-5 | COMPUTED |
| K-40 | 1.251 | 0.110 | 2.847e-5 | 2.910e-5 | COMPUTED |

Dye's inputs (Qh, mu, Dye's h values) are READ. Per element at present isotopic composition (Dye Table 3, READ): U 98.5 microW/kg, Th 26.3, K 3.33e-3; mine (COMPUTED): U 98.48, Th 26.27, K 3.48e-3 microW/kg (the K difference, 4.5 percent, was not traced to a cause).

Per ppm by mass of the element, today's isotopic mix: U 9.85e-11 W/kg, Th 2.63e-11, K 3.48e-15 (COMPUTED). So at the same ppm, U produces 3.7 times Th.

#### 2.8b Reservoir concentrations and heat

Concentrations are ppm by mass. Heat is in pW/kg at the present isotopic mix (COMPUTED from the concentrations and 2.8a). Concentration sources as in table 2.3 (2ND-HAND), except CI K, Th, U which are Palme and O'Neill via pyrolite (2ND-HAND).

| Reservoir | K ppm | Th ppm | U ppm | Heat pW/kg | Th/U | Label |
|---|---|---|---|---|---|---|
| CI chondrite | 546 | 0.030 | 0.0081 | 3.5 | 3.70 | COMPUTED |
| Bulk silicate [R3417] | 260 | 0.0849 | 0.0229 | 5.4 | 3.71 | COMPUTED |
| Bulk silicate Earth (pyrolite), McDonough and Sun 1995, via pyrolite `PY_McDonoughSun1995` | 240 | 0.0795 | 0.0203 | 4.9 | 3.92 | COMPUTED |
| Bulk silicate Earth in the repo (Turcotte and Schubert) | 310 | 0.124 | 0.031 | 7.4 | 4.00 | COMPUTED |
| Depleted mantle [R4775] | 51 | 0.0079 | 0.0032 | 0.7 | 2.47 | COMPUTED |
| N-MORB | 1162 | 0.252 | 0.083 | 18.8 | 3.04 | COMPUTED |
| Upper continental crust | 23244 | 10.5 | 2.7 | 622.6 | 3.89 | COMPUTED |
| Middle continental crust | 19093 | 6.5 | 1.3 | 365.2 | 5.00 | COMPUTED |
| Lower continental crust | 5064 | 1.2 | 0.2 | 68.8 | 6.00 | COMPUTED |
| Bulk continental crust | 14943 | 5.6 | 1.3 | 327.1 | 4.31 | COMPUTED |

Bulk silicate Earth totals for a bulk silicate mass of 4.03e24 kg (the mass Dye uses): Palme and O'Neill 21.7 TW, McDonough and Sun 19.8 TW, the repo values 29.8 TW (all COMPUTED). Dye's own ranges (Table 7, READ): chondritic-composition model 13 TW; cosmochemical model 11 to 21 TW (U 12 to 22 ng/g, Th 45 to 83, K 140 to 260 ug/g); geophysical model 19 to 38 TW (U 20 to 35, Th 77 to 135, K 240 to 420). [R3293] ([R3293], READ) assume 260 ppm K, 23 ppb U, 85 ppb Th, giving 22 TW, of which about 15 TW sits in the mantle.

Refractory and volatile behaviour (COMPUTED from the table above): Th/U is the same in CI (3.70) and the bulk silicate Earth (3.71), so Th and U track each other. K/U is 67 407 in CI and 11 354 in the bulk silicate Earth, a depletion factor of 0.168: K is volatile and must be drawn separately from Th and U. Th and U are lithophile and refractory, so they follow Al, Ca and Ti; the bulk silicate Earth is enriched over CI by a factor of about 2.8 in Th (0.0849 over 0.030), the same factor that applies to Al2O3 (4.5 over 1.59, 2.8). Use that: Th and U in the silicate part scale with the refractory (Al) enrichment of the body.

Condensation temperatures (50 percent, nebular pressure 1e-4 bar, from [R2733: Table 1], READ): Th 1659 K, Al 1653 K, U 1610 K, Ti 1582 K (all refractory, host phase hibonite or titanite); Na 958 K and Rb 800 K (moderately volatile). Potassium is not in the part of that table that could be read; its value of 1006 K is now READ in [R2726] (feldspar, solar-system composition; see part H section 3.4). The ordering Th, U above Mg and Si (1354 K), and Na and K well below them, is the physical reason for the volatile K and the refractory Th and U. The solar and CI abundance tables of [R0169] were also downloaded; the paper states that its proto-solar values for radioactive nuclides come from [R2733] (READ, text).

#### 2.8c Heat with time

Concentrations fixed to the Palme and O'Neill bulk silicate Earth at 4.567 Gyr; "age" is time since formation; each isotope's amount scales as `exp(lambda (4.567 - age))`. Total per kg and share by element. COMPUTED.

| Age Gyr | H / H(4.567) | H W/kg | K share | Th share | U share | Label |
|---|---|---|---|---|---|---|
| 0.0 | 4.98 | 2.69e-11 | 0.42 | 0.10 | 0.47 | COMPUTED |
| 0.5 | 3.80 | 2.05e-11 | 0.42 | 0.13 | 0.45 | COMPUTED |
| 1.0 | 2.98 | 1.61e-11 | 0.41 | 0.17 | 0.43 | COMPUTED |
| 2.0 | 1.98 | 1.07e-11 | 0.35 | 0.24 | 0.41 | COMPUTED |
| 3.0 | 1.44 | 7.76e-12 | 0.28 | 0.31 | 0.41 | COMPUTED |
| 4.567 | 1.00 | 5.39e-12 | 0.17 | 0.41 | 0.42 | COMPUTED |
| 6.0 | 0.79 | 4.24e-12 | 0.10 | 0.49 | 0.41 | COMPUTED |
| 8.0 | 0.61 | 3.29e-12 | 0.04 | 0.57 | 0.39 | COMPUTED |
| 10.0 | 0.50 | 2.68e-12 | 0.02 | 0.64 | 0.35 | COMPUTED |
| 12.0 | 0.42 | 2.24e-12 | 0.01 | 0.69 | 0.30 | COMPUTED |

Isotope fractions at formation (4.567 Gyr ago, COMPUTED): K-40/K 1.47e-3 (12.6 times today), U-235 24.3 percent of U (90 times today's U-235 amount), U-238 75.7 percent of U; Th-232 amount 1.25 times today's.

Cross-checks: Nimmo et al. state heat production was 3.5 times higher about 4 Gyr ago ([R3293], READ); mine gives 3.67 at an age of 0.567 Gyr (COMPUTED), within 5 percent. [R4469] ([R4469], READ) scale back to a formation heat of 51.6 TW against a present 16.4 TW (ratio 3.15); my U-plus-Th-only ratio is 3.45 (K-40 left out) and my total ratio is 4.98 (COMPUTED). The 10 percent gap against Unterborn's ratio is not explained and is listed in section 5; the large K-40 share at early times is real physics, so the generator should include K-40 for young bodies.

#### 2.8d Link to the star (agent A)

| Statement | Value | Source | Label |
|---|---|---|---|
| Th abundance in Sun-like stars, relative to the Sun | 59 to 251 percent (log eps_Th); measured log eps_Th from -0.14 +/- 0.04 to 0.49 +/- 0.05 | [R4469] abstract and results | READ |
| Th and U relative to Mg between Sun-like stars | 30 to 300 percent of solar (Eu/Mg range -0.5 to +0.5 dex, 1 to 99 percent of the Hypatia catalogue; 10 to 90 percent range about -0.3 to +0.3 dex) | [R3293] | READ |
| Th and U are r-process products with star-to-star scatter, K is not | text | same | READ |
| Eu/Mg decreases with stellar age by about 0.1 dex per 10 Gyr | text | same | READ |
| Stellar Mg/Si scatter is only tens of percent | text | same | READ |
| Decay correction from observed to formation Th | Th x 1.051 (1 Gyr), 1.160 (3), 1.280 (5), 1.484 (8), 1.638 (10), 1.808 (12); U-238 x 1.17, 1.59, 2.17, 3.46, 4.72, 6.43 | `exp(lambda age)` with half-lives from `radioactivedecay` | COMPUTED |
| Suggested draw (parametrisation) | `x = [Th/Mg]` offset, normal with sigma 0.22 dex, truncated to +/- 0.5 dex, same offset for U | chosen so the 1 to 99 percent range is about +/- 0.5 dex (2.33 sigma = 0.51 dex) | COMPUTED |

#### 2.8e Heat-flow formula, Urey ratio, targets

Formula (a derivation from the pieces above; all terms COMPUTED or READ):

`H(t) = M_sil * [ C_K * h_K * f_K(t) + C_Th * h_Th * f_Th(t) + C_U * (w_238 * h_238 * f_238(t) + w_235 * h_235 * f_235(t)) ]`

`f_i(t) = exp(lambda_i * (T_ref - t))` where `T_ref` is the age at which the concentrations `C` are quoted; `w_238 = 0.9928, w_235 = 0.0072` at present, so at age t the stored uranium splits by `exp` factors. `M_sil` is the mass of the part that carries the heat (bulk silicate: planet mass times (1 - core mass fraction) from part B). Mean surface flux from radiogenic heat: `q_rad = H / (4 pi R^2)`. Total flux: `q = q_rad / Ur`, with Ur the Urey ratio, so the remainder `q - q_rad` is secular cooling.

Urey ratio targets (COMPUTED from [R3293], READ): total surface heat flow 42 to 47 TW, bulk silicate radiogenic heat 22 TW, so the present Urey ratio is 0.47 to 0.52; they quote 35 to 40 TW when the crustal contribution is removed. Borexino geoneutrinos give a total of 38.2 (+13.6, -12.7) TW with a mantle part 24.6 (+11.1, -10.4) TW (same paper, READ). The repo's 0.67 is above all of these. Ur varies with age: young bodies have Ur well below 1 because secular cooling is large; the generator should draw it, not fix it (section 3).

Worked example, an Earth-mass body with core mass fraction 0.325 (a worked example only, not a default): `M_sil = 0.675 * 5.97e24 kg`, Palme and O'Neill concentrations, `H = 21.7 TW`, `q_rad = 0.0426 W/m2` over `R = 6371 km`, `q = 0.091 W/m2` at Ur 0.47 (COMPUTED). The measured total of 42 to 47 TW over 5.1e14 m2 is 0.082 to 0.092 W/m2 (COMPUTED).

#### 2.8f Short-lived heating

| Item | Value | Source | Label |
|---|---|---|---|
| 26Al/27Al at formation (Solar System) | (5.23 +/- 0.13)e-5 | [R2768: Table 2] | READ |
| Range explored | 0.1 to 10 times that value | [R2678] | READ |
| Decay energy and mean life | E = 3.12 MeV, tau = 1.03 Myr | [R2678: Eq. 1] | READ |
| Heat source law | `H = f_Al * (26Al/27Al)_0 * E / tau * exp(-t/tau)` | [R2678: Eq. 1] | READ |
| Energy per decay: 26Al, 60Fe | 5.00e-6 and 4.16e-6 erg | [R2768: Table 4] | READ |
| Total decay energy of 26Al in CI rock | 5.07e10 erg/g (5.07e6 J/kg) | [R2768: Table 4] | READ |
| Same, recomputed (Al 0.84 weight percent, 5.23e-5) | 4.90e6 J/kg (3.4 percent below Lugaro) | my calculation | COMPUTED |
| Initial power, CI Al, ratio 5.23e-5 | 1.51e-7 W/kg; 5.7e-8 at 1 Myr; 2.2e-8 at 2 Myr; 8.2e-9 at 3 Myr; 1.2e-9 at 5 Myr; 9.2e-12 at 10 Myr | my calculation | COMPUTED |
| Initial power, bulk-silicate Al (2.38 weight percent), 5.23e-5 | 4.27e-7 W/kg (about 80 000 times today's long-lived heat 5.4e-12 W/kg) | my calculation | COMPUTED |
| 26Al power falls below today's long-lived heat after about | 10 to 11 Myr (CI Al) | my calculation | COMPUTED |
| 60Fe/56Fe at formation | (1.01 +/- 0.27)e-8 recommended; values 1e-7 to 1e-6 also reported | [R2768: Table 2] and text | READ |
| Interstellar-medium 60Fe/56Fe from gamma-ray data | 0.8 to 13 e-7 | [R2768] text | READ |
| 60Fe half-life | about 2.6 Myr; the `radioactivedecay` ICRP-107 value is 1.5 Myr and is outdated for this isotope | Lugaro text (2.6 Myr, SNIPPET in the search result for the same review); package value READ from the package | SNIPPET |
| Total 60Fe energy at ratio 1e-8, CI Fe 18.66 weight percent | 8.4e3 J/kg (independent of the half-life); initial power 7.0e-11 W/kg with 2.62 Myr half-life (1.2e-10 with the package value) | my calculation | COMPUTED |

Consequences: the 26Al energy per kg is 600 times the 60Fe energy at the recommended ratios (4.9e6 over 8.4e3), so 60Fe matters only for ratios of 1e-7 or more or for bodies that stay hot for 10 Myr or more. Because 26Al/27Al varies between star-forming regions, use the 0.1 to 10 times range as a log-uniform draw (a parametrisation) and take Al from part B.

### 2.9 Surface radiation

#### 2.9a Ground gamma and radon

Dose-rate conversion (air kerma per unit specific activity, infinite flat source, 1 m height, nGy/h per Bq/kg): K-40 0.042; U series 0.463; Th series 0.604, 0.042 / 0.444 / 0.592 ([R2836] text, which also quotes Gasser's 0.036 / 0.357 / 0.482). Source: [R2836], READ.

Specific activity per ppm of the element (COMPUTED with `radioactivedecay` half-lives, N_A and isotope abundances): U 12.35 Bq/kg per ppm (U-238 only; the series is assumed in equilibrium), Th 4.06 Bq/kg per ppm, K 316 Bq/kg per weight percent K (K-40 abundance 1.17e-4). Dose per ppm and per percent (Saito and Jacob factors, COMPUTED): U 5.72 nGy/h per ppm, Th 2.45 nGy/h per ppm, K 13.3 nGy/h per weight percent. With the Malins factors: U 5.48, Th 2.40.

Dose from the reservoirs (COMPUTED from table 2.8b and the factors above; nGy/h; 8766 hours per year):

| Surface rock | Gamma nGy/h (Saito and Jacob) | Gamma nGy/h (Malins) | mGy/yr air kerma | Label |
|---|---|---|---|---|
| Depleted mantle (ultramafic) | 0.11 | 0.10 | 0.0009 | COMPUTED |
| Bulk silicate Earth | 0.68 | 0.67 | 0.0060 | COMPUTED |
| CI chondrite | 0.85 | 0.84 | 0.0074 | COMPUTED |
| N-MORB (basalt) | 2.6 | 2.6 | 0.023 | COMPUTED |
| Lower crust | 10.8 | 10.7 | 0.095 | COMPUTED |
| Bulk continental crust | 41.0 | 40.4 | 0.360 | COMPUTED |
| Middle crust | 48.7 | 48.1 | 0.427 | COMPUTED |
| Upper continental crust (felsic) | 72.1 | 70.9 | 0.632 | COMPUTED |
| World soil mean check: 33 Bq/kg U, 45 Th, 420 K | 60.1 | n/a | 0.527 | COMPUTED |

The soil activities 33 / 45 / 420 Bq/kg were the UNSCEAR world means as RECALLED; UNSCEAR 2008 Annex B is now READ (part H section 3.3): 33 (238U), 32 (226Ra), 45 (232Th) and 412 (40K) Bq/kg, so 412 replaces 420, and the last row is a sanity check only. Spread from ultramafic to felsic: a factor of about 650. The repo's `GROUND_MSV_YR = 0.48` mSv/yr agrees with UNSCEAR's total terrestrial external effective dose, indoors plus outdoors (0.48 mSv/yr, UNSCEAR 2008 Table 12, READ in part H section 3.3): it includes buildings, which a wild planet does not have. Outdoor, full-time, 0.7 Sv/Gy (RECALLED) from the 60 nGy/h soil mean is 0.37 mSv/yr.

Radon: the repo adds a flat 0.5 mSv/yr above 1000 Pa. No source for radon dose per unit of uranium, emanation fraction or atmospheric mixing was reached (NOT REACHED). Radon in air scales with Ra-226 in the soil, so with U, and needs open soil and an atmosphere that holds radon; the rule in section 3 is a placeholder scaling.

Effective dose quality factor: GCR absorbed dose and the repo's dose-equivalent anchors differ. The Moon anchor 500 mSv/yr against absorbed 80 to 99 mGy/yr in lunar orbit (RADOM 0.227 mGy/day at 100 km, 0.257 mGy/day at 200 km; CRaTER 0.22 to 0.27 mGy/day; all READ from [R4718]; conversion to per year COMPUTED) implies a quality factor near 5.5 (COMPUTED). The Mars anchor 237 mSv/yr against Dartnell's absorbed 62 mGy/yr at 16 g/cm2 implies 3.8 (COMPUTED). Both factors depend on the quality factors used and are noted, not fixed. An orbit dose is not a surface dose (the Moon hides part of the sky from an orbiter), so the lunar ratio is rough.

#### 2.9b Photon shielding

Total mass attenuation coefficient mu/rho (cm2/g) of mixtures, from the NIST XCOM element tables (READ) combined with the mixture rule `mu/rho = sum w_i (mu/rho)_i` (COMPUTED; the rule reproduces the NIST water entry exactly). Oxide-derived compositions use the oxide table 2.3 without Mn and P.

| Material | 10 keV | 100 keV | 1 MeV | 10 MeV | Half-value layer at 100 keV g/cm2 | Label |
|---|---|---|---|---|---|---|
| Water or ice | 5.329 | 0.1707 | 0.07072 | 0.02219 | 4.06 | COMPUTED |
| CO2 | 4.975 | 0.1541 | 0.06369 | 0.02054 | 4.50 | COMPUTED |
| N2 | 3.879 | 0.1529 | 0.06364 | 0.02024 | 4.53 | COMPUTED |
| Iron | 170.6 | 0.3717 | 0.05995 | 0.02994 | 1.86 | COMPUTED |
| Sulfur | 50.12 | 0.2020 | 0.06373 | 0.02589 | 3.43 | COMPUTED |
| Basalt (N-MORB oxides) | 36.07 | 0.1900 | 0.06309 | 0.02347 | 3.65 | COMPUTED |
| Felsic (UCC oxides) | 27.6 | 0.1789 | 0.06323 | 0.02301 | 3.87 | COMPUTED |
| Ultramafic (DMM oxides) | 28.51 | 0.1809 | 0.06322 | 0.02301 | 3.83 | COMPUTED |
| SiO2 | 19.01 | 0.1684 | 0.06367 | 0.02263 | 4.12 | COMPUTED |

At 1 MeV all materials differ by under 6 percent per gram (Compton scattering depends on electron count per gram), so the gamma attenuation can use one mass coefficient above 300 keV. Ground gamma comes from the top metre or so of rock (about 100 to 300 g/cm2); an atmosphere does not change it unless the column is above roughly 30 g/cm2 at the source.

#### 2.9c Proton stopping

CSDA range (g/cm2) from NIST PSTAR (READ; tables retrieved for each material):

| Proton energy MeV | Water | SiO2 | Al2O3 | Iron | CO2 | N2 | Air | Label |
|---|---|---|---|---|---|---|---|---|
| 10 | 0.123 | 0.1573 | 0.1618 | 0.2064 | 0.1404 | 0.1393 | 0.1408 | READ |
| 30 | 0.8853 | 1.099 | 1.128 | 1.377 | 1.004 | 0.9982 | 1.007 | READ |
| 100 | 7.718 | 9.39 | 9.621 | 11.37 | 8.722 | 8.68 | 8.744 | READ |
| 250 | 37.94 | 45.69 | 46.79 | 54.54 | 42.79 | 42.61 | 42.9 | READ |
| 1000 | 325.4 | 388.1 | 397.9 | 459.6 | 366.3 | 365.0 | 367.2 | READ |
| 10000 | 4700 | 5493 | 5663 | 6389 | 5044 | 5030 | 5054 | READ |

Proton energy that is just stopped by a column X (g/cm2, from the same tables, COMPUTED by interpolation):

| X g/cm2 | Water MeV | SiO2 MeV | Al2O3 MeV | Iron MeV | CO2 MeV | N2 MeV | Air MeV | Label |
|---|---|---|---|---|---|---|---|---|
| 1 | 32.09 | 28.45 | 28.04 | 24.99 | 29.93 | 30.03 | 29.88 | COMPUTED |
| 10 | 115.7 | 103.6 | 102.2 | 92.89 | 108.0 | 108.3 | 107.9 | COMPUTED |
| 100 | 452.3 | 403.4 | 397.5 | 362.1 | 419.5 | 420.5 | 418.9 | COMPUTED |
| 1000 | 2393 | 2074 | 2028 | 1814 | 2181 | 2187 | 2177 | COMPUTED |

Check of the repo's range law `R = 0.0022 E^1.77` against PSTAR water (ratio law over PSTAR, COMPUTED): 10 MeV 1.05; 30 MeV 1.02; 60 MeV 1.00; 100 MeV 0.99; 200 MeV 1.00; 400 MeV 1.08; 1000 MeV 1.38. The law is good to 5 percent up to 400 MeV (columns up to about 100 g/cm2 of water). At X = 1000 g/cm2 the law gives a cutoff of 2.33 GV against 3.20 GV from PSTAR water (COMPUTED); that matters only for the stellar-particle dose, which is negligible there. Rock needs 1.22 times the mass of water for the same proton range at 100 MeV (9.39 over 7.718, COMPUTED).

#### 2.9d Hadronic shielding and the column rule

Nuclear interaction length lambda_I in g/cm2 (Particle Data Group atomic and nuclear properties pages, READ). It is a property per gram, independent of phase and density, so one value serves ice, liquid and vapour.

| Material | lambda_I g/cm2 | Relative to air | Relative to standard rock | Label |
|---|---|---|---|---|
| Air (dry, 1 atm) | 90.1 | 1.000 | 1.124 | READ |
| N2 | 89.7 | 1.004 | 1.129 | READ |
| CO2 | 88.9 | 1.013 | 1.139 | READ |
| Water | 83.3 | 1.082 | 1.216 | READ |
| SiO2 (fused quartz) | 97.8 | 0.921 | 1.036 | READ |
| MgO | 98.1 | 0.918 | 1.032 | READ |
| CaO | 109.5 | 0.823 | 0.925 | READ |
| Aluminium | 107.2 | 0.841 | 0.945 | READ |
| Silicon | 108.4 | 0.831 | 0.935 | READ |
| Standard rock (PDG definition, 2.65 g/cm3) | 101.3 | 0.889 | 1.000 | READ |
| Sulfur | 112.4 | 0.802 | 0.901 | READ |
| Iron | 132.1 | 0.682 | 0.767 | READ |
| Basalt (N-MORB oxides, `1/lambda = sum w/lambda` with Al for Al2O3 and Fe for FeO) | 103.3 | 0.872 | 0.981 | COMPUTED |
| Felsic (UCC oxides, same rule) | 101.1 | 0.891 | 0.998 | COMPUTED |
| Ultramafic (DMM oxides, same rule) | 100.8 | 0.894 | 0.995 | COMPUTED |

Two independent measures of how much a gram of material is worth in air-equivalent shielding. Nuclear: `lambda_air / lambda_material`. Ionisation: PSTAR range of air over the material's range at 30 to 1000 MeV (COMPUTED). Rock: 0.889 (nuclear) and 0.92 to 0.95 (ionisation). Water: 1.082 and 1.13 to 1.14. Iron: 0.682 and 0.73 to 0.80. Rock and water agree within 3 to 6 percent; iron within 7 to 17 percent, so iron needs an energy-dependent weight if it matters. [R0996] (Biogeosciences 4, 545; PDF read) found the cosmic-ray dose-depth peak at 40 g/cm2 for both 1 and 3 g/cm3 dry regolith and at 30 g/cm2 for ice: a ratio of 0.75, against 0.82 from both lambda_I and PSTAR range (COMPUTED), a difference within the reading precision of a peak position.

Column rule (a derivation from the table above; the structure is Dartnell's finding that atmosphere and ground columns add as shielding depth, READ, and the weights are COMPUTED):

`X_eff = sum over layers k of X_k * (lambda_ref / lambda_k)`, with `X_k` in g/cm2, `lambda_ref` the interaction length of the medium in which the dose law was calibrated. The repo's law is calibrated on Earth air, so `lambda_ref = 90.1`. Then a gram of N2 or CO2 counts as 1.00 to 1.01, ice 1.08, rock 0.89, iron 0.68. The dose law of the repo is then evaluated at `X_eff`, not at the atmosphere column alone.

#### 2.9e Depth-dose shape (flagged)

| Fact | Value | Source | Label |
|---|---|---|---|
| Mars surface dose under 16 g/cm2 atmosphere, solar minimum, 3 g/cm3 dry dust | 6.2 cGy/yr; peak 6.6 cGy/yr at 40 g/cm2 ground depth (1.06 times) | [R0996] | READ |
| Other models | 20 cGy/yr peak at 25 g/cm2; 19.4 surface and 24.9 peak at 30 g/cm2; factor 3 between them | same | READ |
| Peak depth for dry dust at 1 g/cm3 and regolith at 3 g/cm3 | both 40 g/cm2 | same | READ |
| Peak depth for ice at 1 g/cm3 | 30 g/cm2 | same | READ |
| Protons below 1 GeV | simple exponential decay with depth, no build-up; higher energy bands peak deeper (100 GeV to 1 TeV band peaks at 1.3 m) | same | READ |
| Energy absorbed | 73 percent of GCR energy (about 0.02 J/cm2/yr) in 20 m of dry regolith, 5 percent in the atmosphere; 76 percent in a pure ice column | same | READ |
| GCR energy shares by primary | protons 73 percent, helium 19 percent, heavier 8 percent | same (CREME96) | READ |
| Proton-only model with a 1.37 weight reproduces the full-ion result to 3 percent in the top 500 g/cm2 | 3 percent | same | READ |
| Solar maximum against minimum | surface and peak dose about half; integrated over 20 m, two-thirds | same | READ |

Parametrisation proposed (not fitted): ground dose at depth `z` g/cm2 equals the surface dose times `B(z) * exp(-z_eff/Lambda)`, with a build-up `B` rising from 1 to about 1.06 over the first 40 g/cm2 and `Lambda` equal to the repo's 164 g/cm2 beyond. The build-up amplitude is tied to a single Mars case. It is NOT REACHED how it changes with the primary spectrum (Moon without atmosphere, Mars, thick atmospheres). Until it is fitted to dose-depth data, the generator should apply the exponential only (no build-up), which is on the safe side by under 10 percent in the Dartnell case.

#### 2.9f Magnetic field, stellar activity, existing law

The repo rules (2.1) cover dipole cutoff, heliosphere compression, flare-driven particle dose and UV. What the ground adds, unmodelled today: gamma by crust type (2.9a), radon scaling, and shielding of surface and subsurface by the ground column (2.9d, 2.9e). The Moon and Mars anchors of the repo are measured values (Chang'E-4 and MSL-RAD per the repo docstring); they stay calibration targets in section 4, and the rule here only gives the extra ground column an air-equivalent weight.

### 2.10 Materials table (machine-loadable)

Columns: `id, family, composition, molar_mass_g_mol, rho_grain_g_cm3, T_melt_K, T_boil_or_sub_K, T_triple_K, P_triple_bar, T_crit_K, P_crit_bar, H_fus_kJ_mol, H_vap_kJ_mol, H_sub_kJ_mol, H_fus_kJ_kg, lambda_I_g_cm2, debye_K, alb_040_050, alb_050_060, alb_060_070, owner, labels`.

Notes: an empty cell means not reached or not applicable. `T_boil_or_sub_K` is the boiling point at 1 atm, or the sublimation point for CO2. `H_vap_kJ_mol` is at the normal boiling point except where the labels column says otherwise. `H_sub_kJ_mol` is at the temperature given in the source. `owner` is `ground` (this part) or `atmosphere` (gas rows are pointers: the atmosphere lane owns gas properties and these rows only supply phase data that the ground needs, such as frost and condensation limits). `labels` is `column=code` pairs separated by semicolons, codes: R = READ, H = 2ND-HAND, C = COMPUTED, S = SNIPPET, Q = RECALLED, N = NOT REACHED. A column absent from `labels` is empty. Sources: NIST WebBook pages for phase points and critical data (R); CRC Handbook tables in the `chemicals` 1.5.2 package for latent heats and normal boiling points (H); PDG for lambda_I (R); [R4232] via `burnman` for mineral density and Debye temperature (H); albedo from table 2.6 (C); molar mass from `periodictable` (C).

```csv
id,family,composition,molar_mass_g_mol,rho_grain_g_cm3,T_melt_K,T_boil_or_sub_K,T_triple_K,P_triple_bar,T_crit_K,P_crit_bar,H_fus_kJ_mol,H_vap_kJ_mol,H_sub_kJ_mol,H_fus_kJ_kg,lambda_I_g_cm2,debye_K,alb_040_050,alb_050_060,alb_060_070,owner,labels
water_ice,ice,H2O,18.015,0.917,273.15,373.15,273.16,0.0061,647.1,220.64,6.01,40.65,,333.6,83.3,,0.832,0.832,0.822,ground,M=C;rho=R;Tm=H;Tb=R;Tt=H;Pt=R;Tc=H;Pc=R;Hfus=H;Hvap=H;Hfus_kg=C;lamI=R;alb=C;debye=N
co2_ice,ice,CO2,44.009,1.56,216.58,194.67,216.58,5.185,304.18,73.8,9.02,16.7,26.1,204.9,88.9,,,,,ground,M=C;rho=Q;Tm=R;Tb=H;Tt=R;Pt=R;Tc=R;Pc=R;Hfus=H;Hvap=R;Hsub=R;Hfus_kg=C;lamI=R;alb=N;debye=N
ch4_ice,ice,CH4,16.043,,90.67,111.66,90.67,0.1169,190.56,,0.94,8.19,9.7,58.6,,,,,,ground,M=C;Tm=R;Tb=H;Tt=R;Pt=R;Tc=H;Hfus=H;Hvap=H;Hsub=R;Hfus_kg=C;lamI=N;alb=N;debye=N
nh3_ice,ice,NH3,17.031,,194.95,239.72,195.41,0.0606,405.4,113.0,5.66,23.33,31.2,332.3,,,,,,ground,M=C;Tm=R;Tb=H;Tt=H;Pt=R;Tc=R;Pc=R;Hfus=H;Hvap=H;Hsub=R;Hfus_kg=C;lamI=N;alb=N;debye=N
n2_ice,ice,N2,28.013,,63.3,77.34,63.14,0.1252,126.2,33.98,0.71,5.57,,25.3,89.7,,,,,ground,M=C;Tm=R;Tb=R;Tt=R;Pt=R;Tc=R;Pc=R;Hfus=H;Hvap=H;Hfus_kg=C;lamI=R;alb=N;debye=N
so2_ice,ice,SO2,64.064,,200.75,263.13,197.64,0.0167,430.34,,7.4,24.94,,115.5,,,,,,ground,M=C;Tm=R;Tb=H;Tt=R;Pt=R;Tc=R;Hfus=Q;Hvap=H;Hfus_kg=Q;lamI=N;alb=N;debye=N
sulfur,element,S,32.06,2.0,388.4,717.8,,,1313,182.08,1.721,45.0,,53.7,112.4,,0.362,0.777,0.824,ground,M=C;rho=R;Tm=Q;Tb=H;Tc=R;Pc=R;Hfus=H;Hvap=H;Hfus_kg=C;lamI=R;alb=C;debye=N
iron,metal,Fe,55.845,7.874,1811,3133.35,,,,,13.81,340,,247.3,132.1,,,,,ground,M=C;rho=R;Tm=Q;Tb=H;Hfus=H;Hvap=Q;Hfus_kg=C;lamI=R;alb=N;debye=N
basalt_powder,rock,N-MORB oxides,,2.9,1400,,,,,,,,,400,103.3,,0.207,0.214,0.224,ground,rho=Q;Tm=Q;Hfus_kg=Q;lamI=C;alb=C;debye=N
basalt_slab,rock,N-MORB oxides,,2.9,1400,,,,,,,,,400,103.3,,0.125,0.121,0.116,ground,rho=Q;Tm=Q;Hfus_kg=Q;lamI=C;alb=C;debye=N
felsic_granite_powder_a,rock,UCC oxides,,2.7,,,,,,,,,,,101.1,,0.688,0.760,0.815,ground,rho=Q;lamI=C;alb=C;debye=N
felsic_granite_powder_b,rock,UCC oxides,,2.7,,,,,,,,,,,101.1,,0.395,0.436,0.451,ground,rho=Q;lamI=C;alb=C;debye=N
ultramafic_dunite_powder,rock,DMM oxides,,3.3,,,,,,,,,,,100.8,,0.629,0.690,0.706,ground,rho=Q;lamI=C;alb=C;debye=N
ultramafic_pyroxenite_powder,rock,DMM oxides,,3.3,,,,,,,,,,,100.8,,0.348,0.371,0.370,ground,rho=Q;lamI=C;alb=C;debye=N
anorthosite_powder,rock,CaAl2Si2O8,278.203,2.7652,,,,,,,,,,,103.3,752.39,0.710,0.734,0.750,ground,M=H;rho=H;debye=H;lamI=C;alb=C
komatiite,rock,ultramafic lava,,,,,,,,,,,,,,,,,,ground,all=N
carbon_black,carbon,C,12.011,,,,,,,,,,,,,,0.019,0.019,0.019,ground,M=C;alb=C;rho=N
kerogen,carbon,CxHy,,,,,,,,,,,,,,,0.013,0.016,0.021,ground,alb=C;rho=N
halite,evaporite,NaCl,58.44,2.17,,,,,,,,,,,,,0.817,0.855,0.876,ground,M=C;rho=Q;alb=C
gypsum,evaporite,CaSO4.2H2O,172.17,2.31,,,,,,,,,,,,,0.901,0.917,0.925,ground,M=C;rho=Q;alb=C
calcite,evaporite,CaCO3,100.09,2.71,,,,,,,,,,,,,0.911,0.941,0.946,ground,M=C;rho=Q;alb=C
kaolinite,clay,Al2Si2O5(OH)4,258.16,2.6,,,,,,,,,,,,,0.687,0.751,0.768,ground,M=C;rho=Q;alb=C
serpentine,clay,Mg3Si2O5(OH)4,277.11,2.55,,,,,,,,,,,,,0.250,0.267,0.241,ground,M=C;rho=Q;alb=C
hematite_fine,oxidised,Fe2O3,159.69,5.26,,,,,,,,,,,,,0.092,0.086,0.085,ground,M=C;rho=Q;alb=C
hematite_lt10um,oxidised,Fe2O3,159.69,5.26,,,,,,,,,,,,,0.041,0.042,0.061,ground,M=C;rho=Q;alb=C
goethite,oxidised,FeOOH,88.85,4.3,,,,,,,,,,,,,0.030,0.122,0.206,ground,M=C;rho=Q;alb=C
ilmenite,metal-bearing,FeTiO3,151.71,4.7,,,,,,,,,,,,,0.055,0.054,0.056,ground,M=C;rho=Q;alb=C
forsterite,mineral,Mg2SiO4,140.691,3.2266,,,,,,,,,,,,809.17,,,,ground,M=H;rho=H;debye=H
fayalite,mineral,Fe2SiO4,203.771,4.4021,,,,,,,,,,,,618.70,,,,ground,M=H;rho=H;debye=H
enstatite,mineral,Mg2Si2O6,200.774,3.2034,,,,,,,,,,,,812.18,,,,ground,M=H;rho=H;debye=H
ferrosilite,mineral,Fe2Si2O6,263.854,4.0014,,,,,,,,,,,,674.48,,,,ground,M=H;rho=H;debye=H
diopside,mineral,CaMgSi2O6,216.547,3.2791,,,,,,,,,,,,781.61,,,,ground,M=H;rho=H;debye=H
albite,mineral,NaAlSi3O8,262.218,2.6104,,,,,,,,,,,,713.78,,,,ground,M=H;rho=H;debye=H
quartz,mineral,SiO2,60.083,2.5384,,,,,,,,,,,97.8,816.33,,,,ground,M=H;rho=H;debye=H;lamI=R
gas_O2,gas,O2,31.998,,54.8,90.2,54.33,,154.58,50.43,0.44,6.82,,,,,,,,atmosphere,M=C;Tm=R;Tb=R;Tt=R;Tc=R;Pc=R;Hfus=H;Hvap=H
gas_CO2,gas,CO2,44.009,,216.58,194.67,216.58,5.185,304.18,73.8,9.02,16.7,26.1,,88.9,,,,,atmosphere,M=C;Tm=R;Tb=H;Tt=R;Pt=R;Tc=R;Pc=R;Hfus=H;Hvap=R;Hsub=R;lamI=R
gas_CO,gas,CO,28.010,,68.1,81.63,68.09,,134.45,34.99,0.833,6.04,,,,,,,,atmosphere,M=C;Tm=R;Tb=R;Tt=R;Tc=R;Pc=R;Hfus=H;Hvap=H
gas_N2,gas,N2,28.013,,63.3,77.34,63.14,0.1252,126.2,33.98,0.71,5.57,,,89.7,,,,,atmosphere,M=C;Tm=R;Tb=R;Tt=R;Pt=R;Tc=R;Pc=R;Hfus=H;Hvap=H;lamI=R
gas_Ar,gas,Ar,39.948,,83.8,87.28,83.8,0.689,150.86,48.98,1.18,6.43,,,,,,,,atmosphere,M=C;Tm=R;Tb=R;Tt=R;Pt=R;Tc=R;Pc=R;Hfus=H;Hvap=H
gas_H2,gas,H2,2.016,,13.95,20.39,13.95,0.0721,33.18,13.0,0.12,0.90,,,,,,,,atmosphere,M=C;Tm=R;Tb=H;Tt=R;Pt=R;Tc=R;Pc=R;Hfus=H;Hvap=H
gas_H2O,gas,H2O,18.015,,273.15,373.15,273.16,0.0061,647.1,220.64,6.01,40.65,,,83.3,,,,,atmosphere,M=C;Tm=H;Tb=R;Tt=H;Pt=R;Tc=H;Pc=R;Hfus=H;Hvap=H;lamI=R
gas_CH4,gas,CH4,16.043,,90.67,111.66,90.67,0.1169,190.56,,0.94,8.19,9.7,,,,,,,atmosphere,M=C;Tm=R;Tb=H;Tt=R;Pt=R;Tc=H;Hfus=H;Hvap=H;Hsub=R
gas_H2S,gas,H2S,34.08,,190.85,212.87,187.66,0.232,373.3,89.70,2.38,18.67,,,,,,,,atmosphere,M=C;Tm=R;Tb=R;Tt=R;Pt=R;Tc=R;Pc=R;Hfus=H;Hvap=H
gas_SO2,gas,SO2,64.064,,200.75,263.13,197.64,0.0167,430.34,,7.4,24.94,,,,,,,,atmosphere,M=C;Tm=R;Tb=H;Tt=R;Pt=R;Tc=R;Hfus=Q;Hvap=H
```

CSV caveats: the rock rows `Tm` (1400 K, a solidus near 1100 to 1150 C) and `Hfus_kg` (400 kJ/kg) are RECALLED placeholders. Sulfur melting point 388.4 K and iron melting point 1811 K are RECALLED. Grain densities marked Q are typical textbook values and are not checked. Albedo columns hold the mean reflectance values of table 2.6, not geometric albedo.

## 3. How the generator should compute or draw it

### 3.1 Two routes

Route 1 (keep the drawn bulk figures, fill a composition that fits). The planet already has mass, radius and density and, once part B lands, a core mass fraction. Pick the ground type from rules (3.2), take the oxide model for that type from table 2.3, draw K, Th, U from the star (3.3) and compute heat, gamma and albedo from them. The oxide table is a choice among end-members, and the draw lives in the type choice and in a few perturbation parameters (mantle FeO and alkali shift, Th/U offset, K depletion).

Route 2 (build up from the basic abundances). Take the star's Mg, Si, Fe, Al, Ca, Na, K, Th, U (part A), form a mantle mineralogy and run a melting model to get the crust composition ([R1719] describes the thermodynamic approach). Then heat, gamma and albedo follow from the computed crust.

Recommendation: Route 1 for the ground, with a Route 2 hook for the K, Th, U and Al inputs. Reasons: (a) the melt models are reliable only for Earth-like mantles ([R1719], READ); (b) a full melting calculation per body is far heavier than anything else in the pipeline and gives a result that cannot be tested outside the Solar System; (c) the end-member table (2.3) already spans the range between ultramafic (almost no heat) and felsic (heat producers 100 times richer); (d) what the later radiation and heat calculations need from the crust is only four numbers (K, Th, U, mean oxide set for the shielding mass coefficient), which Route 1 gives. Keep the interface so that part B's mantle composition can move the end-member oxides (more FeO gives lower Mg#, more alkali shifts the solidus) without a change in the ground code.

### 3.2 Choosing the ground type (decision rules)

Inputs already stored: mass, radius, density, [Fe/H], age, surface temperature, orbital distance, snow line, `water_mass_fraction`, atmospheric pressure, mobile-lid chance, surface ocean flag. New inputs needed: core mass fraction and mantle FeO (part B), ice inventory (parts C and D), tidal heating (part D), surface oxidation state (atmosphere lane).

1. Body without a solid surface (gas or ice giant): no ground, every field empty.
2. Ice-dominated surface (parts C and D give ice mass fraction above rock in the outer layer, and surface temperature is below the ice stability limit: below 220 K in vacuum per the repo; below 273 K under an atmosphere): ground type `ice`, with species from the ice inventory (water first, then CO2, NH3, CH4, N2 by condensation temperature at the body's surface temperature, from the phase data of 2.10). Albedo from the ice rows with a dust-and-organics darkening draw.
3. Silicate body. Primary crust from magma-ocean solidification: ultramafic if mantle overturn happened, feldspathic or anorthositic if not (Hu 2012). The overturn rule is a parametrisation: draw `overturn` with a probability tied to mantle FeO (more FeO means denser late cumulates and overturn) and not to a fixed threshold; the plagioclase-flotation pressure limit is RECALLED and listed in section 5. Primary crust is kept only on bodies with no later volcanism.
4. Secondary crust, basaltic: present if the body was convecting and melting; with the stored `internal_heat_flux_w_m2` and age, require the flux above a floor `q_min` (parametrisation; log-uniform draw 1e-4 to 1e-3 W/m2 to be tuned, not a published value). Decompression melting works in any convecting mantle [R1719], so the floor is a convenience.
5. Tertiary crust, granitoid: needs tectonic reprocessing with water (Hu 2012). Draw it with the existing `mobile_lid_chance`; the areal fraction `f_felsic` is a uniform draw in 0 to 0.5 (parametrisation) when the lid is mobile, 0 otherwise.
6. Modifications of the dominant silicate crust, each an areal fraction: clay if liquid water ever existed on the surface (ocean flag or past ocean), Fe-oxidised if the atmosphere is oxidising (O2 partial pressure above the repo's ozone floor 0.02 kPa) or if the surface carries water-vapour photolysis products, evaporite only where a dried basin is likely (small areal fraction).
7. Special: sulfur and SO2 frost if tidal heating exceeds the radiogenic flux and the part D composition carries sulfur; carbon-rich crust (graphite or carbide) only if the star's C/O is above a threshold from part A (rule NOT REACHED; threshold in section 5).
8. Grain state: lava and young surfaces (age since resurfacing below about 0.1 Gyr; parametrisation) use slab or crushed rows; older airless surfaces use powder rows with the space weathering factor of 2.7; surfaces under atmosphere use the powder row without weathering darkening beyond what the atmosphere lane assigns.
9. The ground colour is the area-weighted mean of band albedos of the chosen rows.

### 3.3 Radiogenic content

Draws (all seeded, in this order is not required):

1. `x_Th`: normal with sigma 0.22 dex, truncated to +/- 0.5 dex, giving the Th and U offset relative to Mg (2.8d). If part A stores a proper [Th/Mg] or [Eu/Mg] from a real catalogue, use that and skip this draw.
2. `C_Th0`, `C_U0` (concentrations at formation in the silicate part): reference at formation = reference bulk silicate Earth values multiplied by `exp(lambda * 4.567)` (Th 1.253, U-238 2.031, U-235 89.7) and by `10^x_Th`, and by the refractory enrichment of the body (Al2O3 of the chosen composition over the bulk silicate Earth's 4.5, per 2.8b).
3. `C_K0`: from the K depletion factor `d_K` against the Th/U track: `K/U` at the body = 67 407 (CI value) times `d_K`, with `d_K` drawn log-uniform between 0.03 and 1 (parametrisation; the Earth value 0.168 lies inside it; the range itself is a placeholder tied to the large K differences of Solar System bodies stated by [R3293], READ, not a measured spread).
4. `Al0` ratio: log-uniform 0.1 to 10 times 5.23e-5 ([R2678], READ), applied only to bodies whose accretion time is a few Myr (parametrisation).
5. `Ur(t)`: draw or compute from `q_total` and `q_rad`: since total flux is `q_rad + q_cool`, store `q_cool` from the planet's thermal history if part B supplies it; else `Ur = 0.47` at 4.567 Gyr scaled by `H(t)/H(4.567)` is a rough placeholder. Do not use 0.67.

Heat: `H(t) = M_sil * sum_i C_i * h_i * f_i(t)` (2.8e) with `t` equal to the body's age. If the body is a rogue planet use the rogue age. Surface radiogenic flux `q_rad = H / (4 pi R^2)`. The crust holds a fraction of the heat producers; for the ground dose only the crust matters (3.4), and for the mantle heat subtract the crust inventory: `H_mantle = H - H_crust` with `H_crust` from the crust mass and its enrichment (continental crust of the bulk silicate Earth holds about 40 percent of Th and U; this share is RECALLED, check section 5).

### 3.4 Surface radiation from composition

1. Ground gamma: `D_gamma = 12.35 * 0.463 * C_U(t) + 4.06 * 0.604 * C_Th(t) + 316 * 0.042 * C_K(t)%` in nGy/h with `C_U`, `C_Th` in ppm and `C_K` in weight percent of the crust at the surface type chosen in 3.2, each at time `t`; for mixed surfaces, area-weight. Convert to mSv/yr with `8766e-6 * q` where `q` is a documented dose-equivalent factor (RECALLED 0.7 Sv/Gy; part of the placeholder).
2. Radon: only when there is an atmosphere above 1000 Pa and soil is open (repo rule): `D_radon = D_radon_ref * C_U(t) / C_U,ref`, with `D_radon_ref` and `C_U,ref` a documented pair to be taken from a primary source (NOT REACHED). Do not draw `a_rad` any more.
3. Column: replace the atmosphere-only column in `gcr_dose_msv_yr` by `X_eff = X_atm * (90.1/lambda_gas) + sum_ground X_k * (90.1/lambda_k)` for a point under or inside the ground. For the surface the ground adds nothing; for a buried habitat at depth `z` use `lambda_rock` of the chosen crust (about 101 g/cm2 for most silicates, table 2.9d). The gas `lambda` comes from the table (N2 89.7, CO2 88.9, or by the gas mix with `1/lambda = sum w/lambda`).
4. Shielding of gamma and X-rays by the ground and atmosphere: `exp(-mu_rho * X)` with the coefficients of table 2.9b.
5. Proton shielding for the stellar-particle law: use the PSTAR table 2.9c for the actual medium, or the repo's range law with a medium multiplier (rock needs 1.22 times water's mass for the same range).
6. The existing magnetic and stellar-activity rules stay as they are (2.1).

### 3.5 What to store per body

| Quantity | Name (proposed) | Unit | Range | Per | Label |
|---|---|---|---|---|---|
| Ground type of the dominant surface | `ground_type` | enum | ultramafic, anorthositic, basaltic, felsic, ice, sulfur, carbon, evaporite, clay, oxidised, metal | body | new |
| Areal mix of surface types | `ground_mix` | fractions | sum 1 | body | new |
| Grain state | `ground_grain` | enum | slab, crushed, powder | body | new |
| Mean band albedo | `ground_albedo_040_050`, `_050_060`, `_060_070` | 1 | 0 to 1 | body | new |
| Radiogenic concentrations at formation | `c_k0_ppm`, `c_th0_ppm`, `c_u0_ppm` | ppm of the silicate | K 20 to 2000, Th 0.01 to 0.5, U 0.003 to 0.15 | body | new |
| Early Al ratio | `al26_ratio` | 1 | 5e-6 to 5e-4 | body | new |
| Radiogenic heat at the body's age | `radiogenic_heat_tw`, `radiogenic_flux_w_m2` | TW, W/m2 | derived | body | derived (could be stored for speed) |
| Urey ratio | `urey_ratio` | 1 | 0.1 to 0.9 | body | new |
| Ground gamma | `ground_gamma_ngy_h` | nGy/h | 0.1 to 100 | body | derived |
| Regolith grain radius and filling factor | `regolith_grain_um`, `regolith_phi` | um, 1 | 10 to 10000, 0.3 to 0.6 | body | new (only if thermal inertia is modelled) |

Already stored and reused: mass, radius, density, [Fe/H], age, surface temperature, orbital distance, atmospheric pressure, gravity, magnetic moment, water mass fraction. `a_rad` is removed in favour of the stored concentrations.

### 3.6 Test invariants

1. Heat consistency: `radiogenic_flux_w_m2 * 4 pi R^2` equals `H(t)` from the stored concentrations to 1 percent, and a doubling of `c_th0_ppm` and `c_u0_ppm` raises both the heat flux and the ground gamma (so there is no independent dose draw).
2. Scaling: `D_gamma(2C) = 2 D_gamma(C)`.
3. Mass balance: oxide sums of each ground type are 100 +/- 1.5 percent, and area fractions sum to 1.
4. Column additivity: with zero ground depth, the dose equals the old atmosphere-only dose to within 1 percent for a pure N2 atmosphere; for a CO2 atmosphere it changes by less than 2 percent (lambda 88.9 against 90.1).
5. Albedo bounds: every band mean lies in 0 to 1, and the area-weighted mean of two types lies between them.
6. Ordering: ground gamma of ultramafic below basaltic below felsic for the same age.

## 4. Calibration and test targets

Measured bodies, used only to check the laws, never as presets.

| Target | Value | Source | Used for | Label |
|---|---|---|---|---|
| Bulk silicate Earth heat at 4.567 Gyr | 21.7 TW (Palme and O'Neill), 19.8 TW (McDonough and Sun); tolerance +/- 5 percent of the chosen model | table 2.8b | heat formula check | COMPUTED |
| Earth radiogenic heat | about 22 TW total, about 15 TW in the mantle ([R3293] via [R3363]) | [R3293] | order of magnitude | READ |
| Earth total heat flow | 42 to 47 TW | [R3293] | Urey ratio check 0.47 to 0.52 | READ |
| Geoneutrino total | 38.2 (+13.6, -12.7) TW | [R3293] | broad check | READ |
| Heat production 4 Gyr ago against present | 3.5 times | [R3293] | time evolution; mine 3.67 | READ |
| Lunar interior heat flow | 0.018 W/m2 | [R1861] | low-Th smaller body | READ |
| Lunar regolith density, conductivity, inertia | 1100 and 1800 kg/m3; 7.4e-4 and 3.4e-3 W/m/K; I = 36.6 to 100 (at 273 K) | [R1861] and table 2.4 | regolith law | READ |
| Mars dose under 16 g/cm2 | 6.2 cGy/yr absorbed, peak 6.6 at 40 g/cm2 | [R0996] | dose-depth shape | READ |
| Lunar orbit absorbed dose | 0.22 to 0.27 mGy/day | [R4718] | free-space GCR | READ |
| Arrokoth geometric albedo | 0.21 (550 nm), 0.24 (610 nm) | [R1970] | dark red ice | READ |
| World soil gamma | 60 nGy/h (from 33 / 45 / 420 Bq/kg; UNSCEAR 2008 gives 33 / 45 / 412 and 58 nGy/h outdoors) | activities now READ (part H) | gamma conversion | COMPUTED |
| Proton range | PSTAR table 2.9c | NIST PSTAR | stopping | READ |
| Elemental ratios of chondrites | K/U 67 407 (CI) and 11 354 (BSE); Th/U 3.7 | table 2.8b | refractory/volatile rule | COMPUTED |
| Asteroid heat sources | Al-26 plus long-lived heat of chondrites | part D section 2.5 | cross-check with D | READ |
| Space weathering timescale | 960 to 2050 Myr | [R4713] | `tau_w` draw range | READ |

## 5. Gaps

### Boss must research (only what is still unreachable)

Ranked by effect on the design.

1. A K depletion survey for rocky bodies: [R2726] (the IOP page opened but only the abstract; the K condensation temperature of 1006 K is RECALLED). Gives the distribution of `d_K`.
2. Original geochemical references behind publisher walls, to turn the pyrolite-package numbers (2ND-HAND) into READ: [R3417] (Treatise on Geochemistry 3.1, Elsevier), [R3842] (Treatise 4.1), McDonough and Sun 1995 (Chem. Geol. 120, 223), [R1494] (Geochem. Geophys. Geosyst. 14, 489, AGU), [R4775]. Elsevier and AGU pages were refused.
3. Komatiite, Mars, Venus and lunar crust oxide tables: Taylor and McLennan), Wanke and Dreibus on Martian soil, Arndt on komatiite; each would give a calibration row.
4. Radon: UNSCEAR 2008 Annex B (Ra-226 to radon emanation, indoor and outdoor dose, soil Ra-226 activities) to replace the flat 0.5 mSv/yr and the 33 / 45 / 420 Bq/kg world soil means. Done in part of GEN.465 (part H section 3.3): the soil means are READ, indoor radon is log-normal with geometric mean 45 Bq/m3 and GSD 2.1; the radon emanation fraction from soil Ra-226 is still not found.
5. Ground-depth dose profile for non-Mars cases: Geant4 or published dose-depth curves for the Moon and thick atmospheres (Dartnell's Fig. 6 data digitised; Schwadron, Wilson or Dartnell follow-ups). Gives the build-up amplitude and Lambda in the ground.
6. Space weathering darkening amplitude (the `d` in 2.7) per mineral type: Hapke-type models and lunar sample reflectance (Pieters et al.; Lucey), and the fresh-to-weathered ratio for anorthositic, basaltic and ultramafic regolith.
7. Visible reflectance of CO2, CH4, NH3, N2, SO2 and tholin frosts and of Mars dust, with grain-size dependence (Grundy and Fink; Cruikshank; Mars Pancam dust papers behind publishers).
8. Debye temperatures and specific heats of water ice, CO2, CH4, NH3, N2, SO2, S and iron (CRC Handbook or Giauque tables), and densities of the solid ices; liquidus and solidus of basalt, peridotite and iron-sulfur at 1 atm, with latent heats (NIST has only pure compounds).
9. The Th/U heat-production time curve discrepancy against [R4469] (their 51.6 TW at formation against my 62 TW without K-40, and 108 TW with it); their Methods (not just the text read) would resolve whether they exclude K-40.
10. Conversion from lab reflectance to geometric and Bond albedo for particulate surfaces ([R1806] book or an arXiv equivalent).
11. Origin[R4232] paper (Stixrude and Lithgow-Bertelloni, Geophys. J. Int. 184, 1180; Oxford Academic refused); the parameters used here were read from open-source code that cites it.
12. Plagioclase-flotation depth/pressure limit and mantle overturn criterion for magma-ocean crusts ([R1247] and later).

### To verify (RECALLED items), ranked

1. Solid grain densities in the CSV for rocks, evaporites, clays, hematite, goethite, ilmenite (Q).
2. Basalt and peridotite solidus 1400 K and latent heat 400 kJ/kg (CSV).
3. Sulfur melting point 388.4 K, iron melting point 1811 K, iron enthalpy of vaporisation 340 kJ/mol, SO2 enthalpy of fusion 7.4 kJ/mol.
4. Dose-equivalent factor 0.7 Sv/Gy and the UNSCEAR 0.48 mSv/yr total terrestrial dose (it explains the repo's `GROUND_MSV_YR`).
5. World soil means 33 / 45 / 420 Bq/kg.
6. Condensation temperature of K, 1006 K (Th, U, Na are READ in 2.8b). Now READ in [R2726] (part H section 3.4).
7. Calorimetric Cp of forsterite near 300 K (about 840 J/kg/K) used to state the Debye overestimate.
8. Continental crust holds about 40 percent of Th and U of the bulk silicate Earth.
9. Gas phase data values in the CSV that are marked H come from CRC files read in the `chemicals` package; a few (CO Hvap, Ar Hvap, H2 Hvap) were not cross-checked with NIST.

## 6. Proposed work items

| # | Title | One-line scope | Depends on | Size |
|---|---|---|---|---|
| E1 | Research: primary crust and radiogenic sources | Obtain the items 1 to 4, 6, 7 of "Boss must research" and relabel table values READ | none | medium |
| E2 | Design: ground and radiogenic section of the composition note | Turn sections 2 and 3 of this part into the design doc section with decision rules | E1 partly; B's core mass fraction | medium |
| E3 | Table: `materials.csv` | Add the CSV of 2.10 to the repo with a loader, label column kept; test that every `labels` entry names an existing column | E2 | small |
| E4 | Table: crust oxide and heat-producer reservoir table | Load table 2.3 and 2.8b as data (K, Th, U, oxides per ground type) | E2 | small |
| E5 | Table: band albedo table | Add the albedo rows (2.6) to the materials table and a grain-state flag | E3 | small |
| E6 | Drawing rule: K, Th, U and Al-26 | Seeded draws of 3.3 from stored [Fe/H] and (when A provides it) [Th/Mg]; stores the new fields | A's abundance fields; E4 | medium |
| E7 | Drawing rule: ground type | Rules of 3.2 using mobile-lid chance, ocean, ice inventory, heat flux | B, C, D interface; E4 | medium |
| E8 | Derivation: radiogenic heat module | `H(t)` isotope-wise with the isotope constants of 2.8a; replaces `ROGUE_RADIOGENIC_ISOTOPES`; Urey ratio draw | E6 | medium |
| E9 | Derivation: ground gamma and radon | Per-ppm dose from 3.4; remove `a_rad`; keep radon flat until a source is found | E6, E7 | small |
| E10 | Derivation: air-equivalent column | `X_eff` rule of 2.9d in `gcr_dose_msv_yr` and `sep_dose_msv_yr` (ground depth argument) | E3 | small |
| E11 | Derivation: band albedo and space weathering | Area-weighted colour, weathering factor of 2.7 | E5, E7 | medium |
| E12 | Derivation: regolith thermal properties | Hayne laws and Gundlach-Blum conduction from grain radius and filling factor; ice laws | E3 | medium |
| E13 | Storage: schema columns | Add the fields of 3.5 and a migration, keeping existing columns until the old path is removed | E6, E7 | medium |
| E14 | Tests | The six invariants of 3.6, plus the calibration targets of section 4 with stated tolerances | E8, E9, E10 | medium |
| E15 | Docs | Update the radiation and magnetism design notes and the composition note to refer to the new quantities | E2, E13 | small |
| E16 | Cross-check with D | Make D's chondrite heating and E's formula use one function; add the Al-26 term of 2.8f | E8, part D | small |

## 7. Verification summary

Counts of labelled entries in this file (tables 2.1 to 2.9 and section 4; the CSV is counted separately by cell group because each CSV row carries a `labels` string).

| Where counted | READ | 2ND-HAND | COMPUTED | SNIPPET | RECALLED | NOT REACHED |
|---|---|---|---|---|---|---|
| Table rows in section 2 (192 rows) | 72 | 9 | 106 | 4 | 0 | 1 |
| Table rows in section 4 (15 rows) | 12 | 0 | 3 | 0 | 0 | 0 |
| CSV cell groups in 2.10 (282 groups on 44 rows; each group is one column or column family named in a `labels` string) | 94 | 74 | 62 | 0 | 26 | 26 |
| Inline RECALLED items in the text (komatiite MgO range, forsterite Cp, K condensation temperature, world soil activities, UNSCEAR 0.48 mSv/yr, 0.7 Sv/Gy, crustal share of Th and U, plagioclase-flotation limit) | 0 | 0 | 0 | 0 | 8 | 0 |
| Total | 178 | 83 | 171 | 4 | 34 | 27 |

Summary: READ plus 2ND-HAND plus COMPUTED is 432 entries; RECALLED plus NOT REACHED is 61; SNIPPET is 4. READ entries come from the arXiv papers listed in the findings (Dye, Unterborn, Nimmo, Lichtenberg, Lugaro, Hayne, Hu, Guimond, Malins, Dartnell, Wimmer-Schweingruber, Willman, Hofgartner, Lodders and others), NIST XCOM, PSTAR and WebBook pages, PDG nuclear-property pages, USGS splib07a and the Paragas Zenodo data, and the planetGen repo files named in 2.1. 2ND-HAND entries are numbers read in the `pyrolite`, `chemicals` and `burnman` package data or in a paper that cites the original. COMPUTED entries were produced with Python from those inputs (scripts kept outside the repo). A CSV label R on a NIST WebBook value means the WebBook page, a compilation, was read, not the original measurement paper cited inside it.
