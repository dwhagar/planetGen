# Part G1: Boss's second note, section 6 (comet volatiles and photochemical lifetimes)

Task G1 of brief G. Each claim in the note is treated as a pointer, not a source. Labels as in the common brief: READ (text or table read in an arXiv PDF or a reachable open file), 2ND-HAND (read in a reachable paper or open code that cites the original), COMPUTED (re-derived here with Python from READ inputs), SNIPPET, RECALLED, NOT REACHED.

Scripts and outputs are in `parts/G-files/` with the prefix `G1_`:
`G1_robinson_stats.py` and `.out.txt` (per-species statistics), `G1_class_test.py` and `.out.txt` (dynamical type against composition), `G1_co_vs_size_and_hypervolatile_corr.py` and `.out.txt` (CO against nucleus size; rank correlations), `G1_heays_solar_lifetimes.py` and `.out.txt` (quiet-Sun photolysis lifetimes).

Sources actually read (all from reachable hosts):

| Short name | Reference | How reached | Used for |
|---|---|---|---|
| HP22 | Harrington [R3535] (CO, CO2, H2O in comets) | arXiv PDF | CO and CO2 ranges, dynamical-class split, bands |
| Rub19 | [R3835] (67P ROSINA, Table 2) | arXiv PDF | 67P bulk ratios |
| Rob24 + DB | Robinson, Malamud, Opitom, Perets, [R3782], and its database (Edinburgh DataShare, [R4972], `df_comp_all.csv`, 909 rows, 69 comets, CC BY 4.0) | arXiv PDF and CSV | per-species statistics; most CO values are Dello [R1063], others [R3835], Biver and Bockelee-Morvan papers |
| Coch15 | [R0835], "The Composition of Comets", Space Sci. Rev. 197, 9, [R0835] | arXiv PDF | taxonomies (optical and IR), CN/OH bands |
| Eis19 | [R1234] | arXiv PDF | chemical-evolution taxonomy |
| Cro16 | [R0937], "Comets at radio wavelengths", [R0937] | arXiv PDF | OH 18 cm, HCN 88.6 GHz, NH3 radio |
| dVB13 | [R1088] (C/2001 Q4 and C/2002 T7 at mm wavelengths), Table 2 | arXiv PDF | photodissociation rates at 1 AU (quiet Sun), CO 230.538 GHz, r^-2 scaling |
| Odin07 | [R0407], Planet. Space Sci., [R0407] | arXiv PDF | NH3 572.498 GHz, NH3 lifetime 6700 s |
| Hat05 | Hatchell, Bird, van [R1844] | arXiv PDF | NH3 inversion lines, lifetime 5400 to 5800 s |
| Hea17 | [R1873] | arXiv PDF | Tables 18, 19, 22: rates, ionisation, branching, spectral dependence |
| sbpy | sbpy 0.6.0 source (PyPI), `activity/gas/data/__init__.py` and `core.py` | files.pythonhosted.org | [R2040] quiet/active values (2ND-HAND), [R2041] |
| Rag20 | Raghuram, [R3647: Table 6] | arXiv PDF | photodissociation frequencies into excited O and C |
| Tegler | Long-slit NH2 paper, NTRS 19930010041 | NTRS PDF | NH3 7.0e3 to 7.7e3 s |
| Cor25 | [R0890] (JWST, 3I/ATLAS) | arXiv PDF | interstellar object |
| Bel26 | [R4752, R4754], [R0336] (JWST MIRI, 3I/ATLAS) | arXiv PDF | CH4 nu3 band |

Not reached: [R0435] (Phil. Trans. R. Soc. A 375, 20160252; royalsocietypublishing.org returns 403, the HAL record has no PDF link, not on PMC), Dello [R1063, R2040, R3206] (PSS 106, 11; ScienceDirect refused), [R0410, R0934, R2041] (no single paper identified). Where a number from them is used below it came through one of the papers in the table above and is labelled 2ND-HAND.

## Verdicts

Rule used for G1.1: CONFIRMED when the observed range at heliocentric distance r_h at or below 2 AU lies inside the claimed range or within a factor 1.5 of its ends; PARTLY otherwise. "Observed" means comet-to-comet spread of the coma production-rate ratio to water (the quantity every table below reports), not a measured nucleus bulk ratio, except 67P.

### G1.1 Abundances relative to H2O (percent, per species)

| id | Claim as written | Verdict | The right figure (percent of H2O) | Source | Label |
|---|---|---|---|---|---|
| G1.1a | CO 0.5 to 30 | PARTLY | 0.3 to 26 within 2.5 AU (32 comets), median 3 +- 1 (3.3 +- 1.3 in the figure caption); mean 6.4 +- 1.3 (Dello [R1063] mean for 24 comets 5.2 +- 1.3). Lowest: 103P/Hartley 2 0.3 (Dello [R1063]), 0.15 to 0.45; 8P/Tuttle 0.4, 2P/Encke 0.43. Highest at r_h at most 2 AU: Hale-Bopp 26. Values above 30 exist only far out (C/2010 G2 Hill 91 at 2.5 AU; C/2006 W3 230 at 3.4 AU; 29P 460 at 6.2 AU) and mostly reflect the water cut-off, not the nucleus. The lower bound 0.5 is too high; 30 is fine for r_h at most 2.5 AU | HP22 abstract, Sect. 3.1, Table 3; Rob24 DB (COMPUTED); Coch15 Sect. 2.3 | READ / COMPUTED |
| G1.1b | CO2 2 to 20 | PARTLY | 2 to 30 within 2.5 AU (23 comets), median 12 +- 2 (HP22 Table 3). In the database at r_h at most 2 AU: 4 to 37 (N = 14, median 11). Hartley 2 about 20 ([R0047] as quoted in Coch15). 67P bulk 4.7 +- 1.4. Upper bound 20 is too low; values near 87 at 6 AU are outside the usable regime | HP22; Rub19 Table 2; Rob24 DB | READ / COMPUTED |
| G1.1c | CH4 0.2 to 1.5 | CONFIRMED | 0.15 to 2.2 at r_h at most 2 AU (N = 12): C/2001 A2 LINEAR 0.15, 67P 0.34 (in situ), Hale-Bopp 1.22, 73P/Schwassmann-Wachmann 3 2.17; inter-quartile 0.36 to 1.08, median 0.95, geometric mean 0.68. (C/2010 G2 Hill 9 at 2.5 AU is an outlier.) The Les Houches table gives comets 0.2 to 1.2 | Rob24 DB (COMPUTED); Rub19; BM-slides | COMPUTED / READ |
| G1.1d | C2H6 0.1 to 1.5 | CONFIRMED | 0.09 to 0.82 at r_h at most 2 AU (N = 17), median 0.39, geometric mean 0.37; 67P 0.29. Holmes 1.85 (2.5 AU) and Hill 3.8 (2.5 AU) are outliers. Slides: 0.4 to 1.2 | Rob24 DB (COMPUTED); Rub19; BM-slides | COMPUTED / READ |
| G1.1e | NH3 0.2 to 1.5 | PARTLY | 0.125 to 0.9 (N = 14, r_h at most 2 AU), median 0.66; the lowest comet is 73P/Schwassmann-Wachmann 3 at 0.125, then Lulin 0.24. 67P bulk 0.67 +- 0.20 (the hemisphere values quoted from Le Roy are 0.06 to 0.15). Odin tentative detections 0.5 and 0.3. Upper bound 1.5 not reached in the database (slides: 0.6 to 1.8). Radio upper limits 0.13 to 0.74 (Hat05) | Rob24 DB (COMPUTED); Rub19; Odin07; Hat05; BM-slides | COMPUTED / READ |
| G1.1f | HCN 0.05 to 0.5 | PARTLY | 0.03 to 0.36 at r_h at most 2 AU (N = 21), median 0.14, geometric mean 0.13; lowest 6P/d'Arrest 0.03, 19P/Borrelly 0.06; highest Hale-Bopp 0.36 (Holmes 0.54 after outburst at 2.5 AU). 67P 0.14 +- 0.04. Range is narrow ([R0936], via Coch15). Lower bound 0.05 too high by a factor 1.7 | Rob24 DB (COMPUTED); Rub19; Coch15 | COMPUTED / READ |

Added species for the weighted table T16 (not claims; same basis, r_h at most 2 AU, from the same database, geometric mean and 1-sigma scatter in dex):

| Species | N | Observed range | Median | Geo-mean | sigma (dex) | 67P bulk (Rub19) | Label |
|---|---|---|---|---|---|---|---|
| H2S | 9 | 0.09 to 1.5 | 1.1 | 0.71 | 0.37 | 1.10 +- 0.46 | COMPUTED / READ |
| CH3OH | 20 | 0.21 to 4.5 | 2.15 | 1.89 | 0.32 | 0.21 +- 0.06 | COMPUTED / READ |
| H2CO | 15 | 0.09 to 1.1 | 0.22 | 0.245 | 0.34 | 0.32 +- 0.10 | COMPUTED / READ |
| C2H2 | 11 | 0.03 to 0.28 | 0.074 | 0.09 | 0.28 | not tabulated here | COMPUTED |
| N2 | 1 comet | no population data in any file read | n/a | n/a | n/a | 0.089 +- 0.024 (Rub19 Table 2) | READ |
| O2 (for the record) | 1 comet | 67P only | n/a | n/a | n/a | 3.1 +- 1.1 | READ |

Note the 67P values are not typical for CH3OH (0.21 against a population median about 2) and are lower than the median for CH4 (0.34 against 0.95). See "Changes to parts A to F".

### G1.2 Photochemical lifetimes at 1 AU (seconds; photodissociation plus photoionisation)

All source values are for an optically thin coma and a stated solar activity; the lifetime scales as r_h squared (rate as r_h^-2) (dVB13 Sect. 3: "photodissociation rates were scaled by the heliocentric distance as r_h^-2", READ).

| id | Claim | Verdict | Quiet Sun | Active Sun | Source | Label |
|---|---|---|---|---|---|---|
| G1.2a | H2O 8.0e4 | CONFIRMED (this is the quiet-Sun value) | 8.29e4 | 4.94e4 | [R2040] as coded in sbpy (`[8.294e4, 4.939e4] s`, ordered quiet, active); also [R0837] 5.2e4 "intermediate activity" in sbpy; Hea17 computed 1.13e5 (unshielded, factor 1.4 higher) | 2ND-HAND (sbpy), COMPUTED (Hea17) |
| G1.2b | CO 1.5e6 | CONFIRMED (quiet Sun) | 1.5e6; dVB13 Table 2 beta0 = 6.5e-7 1/s from [R0934] gives 1.54e6 | not given | sbpy; dVB13 Table 2 | 2ND-HAND / READ |
| G1.2c | CO2 5.0e5 | CONFIRMED (quiet Sun) | 4.95e5 (HM15), 5.0e5 (CE83) | 2.10e5 | sbpy. Hea17 computed 1.2e6 is 2.4 times longer; the CO2 to O(1D) channel alone is 1.6e-6 1/s (Rag20 Table 6), so 5e5 is the safer value | 2ND-HAND; COMPUTED |
| G1.2d | CH4 1.2e5 | PARTLY | 1.2e5 appears in the 1976 table (Jac76 Table 1, from Samson cross sections) but no modern comet table was reached; Hea17 computed quiet-Sun total 1.9e5 (photodissociation 2.0e5, ionisation adds 7 percent) | active Sun about 1.7 to 2.4 times shorter (see below) | Jac76 (READ, OCR of exponent partly lost); Hea17 (COMPUTED) | READ / COMPUTED |
| G1.2e | C2H6 1.0e5 | UNVERIFIABLE as stated | no 1 AU comet table found. Hea17 computed 1.45e5 (photodissociation 1.66e5, ionisation adds 13 percent). Tried: sbpy data (no C2H6), dVB13, Cro16, Odin07, Jac76 Table 1 (not listed), NTRS search; missing: [R2041] [R0934] | | | COMPUTED |
| G1.2f | NH3 2.5e4 | CONTRADICTED | 6.7e3 (Odin07, "6700 s at 1 AU"); 5.4e3 to 5.8e3 near solar minimum (Hat05); 7.0e3 to 7.7e3 (used in the Tegler paper); Hea17 computed 6.5e3. The 1976 table gives a value of order 2e3 (exponent garbled by OCR). Use 6.7e3 (range 5.4e3 to 7.7e3). Boss's value is about 3.7 times too long. Reconciliation with the atmosphere lane's 5.6e3 s (GEN.496, 6.5e3 to 5.6e3 range): their 5.6e3 comes from k = 1.8e-4 per s, the mean of quiet and active Sun ([R2041] and HM15 as quoted by [R1163]); mine is the quiet Sun from [R1873] (rate 1.55e-4 per s from cross sections times the quiet solar field, 91.2 to 200 nm only, scale 37 700) and [R0407]. Ratio 1.2, no disagreement in data; the HM15 quiet and active NH3 rates themselves were not read by either lane | | Odin07, Hat05, NTRS 19930010041, Hea17 | READ / COMPUTED |
| G1.2g | HCN 7.0e4 | CONFIRMED (quiet Sun) | 6.7e4 ([R0934] via dVB13: beta0 = 1.5e-5 1/s). A search snippet puts the range at 3.2e4 to 7.7e4 over solar activity ([R2041], not read). Hea17 computed 1.2e5 (lower rate) | likely 3.2e4 | dVB13 Table 2 | READ; SNIPPET for range |

Other species from the same sources (dVB13 Table 2, quiet Sun, 1 AU, from [R0934]): CH3OH 1.3e-5 1/s (7.7e4 s), H2CO 2.0e-4 1/s (5.0e3 s; as a parent), HNC 1.9e-5, CS 2.5e-5 (parent), CN ([R2041] via sbpy) 3.15e5 s quiet and 1.35e5 s active. Hea17 computed: H2CO 5.2e3 s, H2S 4.3e3 s, N2 2.2e6 s, CH3OH 1.2e5 s.

Solar activity: sbpy's two-element arrays are (quiet Sun, active Sun); the active/quiet shortening is 1.68 for H2O, 2.35 for CO2, 2.33 for CN. Rag20 (quoting [R2041]) says rates vary by about a factor 2 between solar minimum and maximum. dVB13 and [R0408] note that rates for species absorbing at 130 to 200 nm vary little with activity; Lyman-alpha and the extreme UV drive the variation. Therefore: Boss's H2O, CO, CO2, HCN are all quiet-Sun values and must be labelled so.

### G1.3 Photolysis products

| id | Claim | Verdict | The right figure | Source | Label |
|---|---|---|---|---|---|
| G1.3a | H2O to OH + H and O(1D) + H2 | CONFIRMED (add branching) | Solar field: OH + H 0.88, O + 2H or O + H2 0.12 (not state-resolved); Lyman-alpha field 0.74 / 0.26. Emission-line model: H2O to O(1D) 1.0e-6 1/s at 1 AU, about 8 percent of the total rate | Hea17 Table 22; Rag20 Table 6 | READ |
| G1.3b | CO to C + O and CO+ + e | CONFIRMED | Photodissociation C + O, rate 3.0e-7 1/s; photoionisation CO+ + e, 2.45e-7 1/s (so about 45 percent of the destruction is ionisation); C(1D) from CO 7.6e-8 1/s | Hea17 Tables 1, 18, 19 (COMPUTED x 37 700); Rag20 | READ / COMPUTED |
| G1.3c | CO2 to CO + O(1D) | PARTLY | CO + O is the dominant channel; O(1D) is one excited state of it (1.6e-6 1/s of a total near 2.0e-6 at the quiet Sun, Rag20 Table 6, baseline-model values; do not sum with the O(1S) channel in the same table). Photoionisation to CO2+ + e is about 2e-7 1/s (Hea17 COMPUTED) | Hea17 Table 1; Rag20 Table 6 | READ |
| G1.3d | CH4 to CH3 + H and CH2 + H2 | CONFIRMED | Both channels are real. Hea17 lists CH3 + H as the main channel; Jac76 gives CH2 + H2 0.42 to 0.52 and CH3 + H 0.48 to 0.58 at 121.6 nm (OCR of the inequalities is partly garbled). A third channel CH + H2 + H exists (RECALLED, not checked) | Hea17 Table 1 and Sect. 4.3.19; Jac76 | READ |
| G1.3e | C2H6 to C2H4 + H2 and CH3 + CH3 | UNVERIFIABLE | Not in any file read: Hea17 gives only the cross section and ionisation fraction for C2H6. Tried Hea17 Sect. 4.3.23, Jac76, Rag20. Missing: Lee et al. or Gans et al. branching at Lyman-alpha. RECALLED: C2H4 + H2 and C2H4 + 2H dominate, CH3 + CH3 and CH4 + CH2 are minor | | RECALLED |
| G1.3f | NH3 to NH2 + H | CONFIRMED | Solar field: NH2 + H 0.99, NH + H2 or NH + 2H 0.01. At Lyman-alpha the split reverses (NH2 + H 0.05, NH + H2 0.95). The NH2 daughter then lives 3.4e4 s | Hea17 Table 22; NTRS 19930010041 | READ |
| G1.3g | HCN to CN + H | CONFIRMED | H + CN is the listed dissociation product (grade A); photoionisation rate 2.8e-7 1/s (COMPUTED) | Hea17 Table 1 | READ |

### G1.4 Diagnostic bands

| id | Claim | Verdict | The right figure | Source | Label |
|---|---|---|---|---|---|
| G1.4a1 | OH 18 cm | CONFIRMED | OH radical lines at 18 cm (Nancay); the 1665 and 1667 MHz frequencies are not printed there | Cro16 Sect. 2 | READ (frequencies RECALLED) |
| G1.4a2 | OH 309 nm | PARTLY | The OH ultraviolet band is quoted as 3080 Angstrom (308 nm); 309 nm is the upper edge | Coch15; HP22 | READ |
| G1.4b | CO 4.67 um and 230 GHz (J = 2-1) | CONFIRMED | 4.67 um fundamental; CO 2-1 at 230.538 GHz | HP22; dVB13 Table 1 | READ |
| G1.4c | CO2 4.26 um | CONFIRMED | nu3 antisymmetric stretch at 4.26 um | HP22 | READ |
| G1.4d | CH4 3.30 um (nu3) | CONFIRMED | nu3 asymmetric stretch "near 3.3 um" | Bel26 | READ |
| G1.4e | C2H6 3.35 um (nu7) | PARTLY | The nu7 band lies in the 3.3 um region (NTRS 20110023457 title "nu7 Band at 3.3 Micrometers"); the 3.35 um head is not printed in any file read | NTRS 20110023457 | READ (3.3 um); 3.35 um SNIPPET |
| G1.4f1 | NH3 1.15 mm sub-mm rotational | CONTRADICTED | The ground-state rotational line JK = 1_0 to 0_0 is at 572.498 GHz (0.524 mm), searched by Odin. The inversion lines (J, K = J, J = 1 to 4) are in the 18 to 26 GHz range (1.2 to 1.7 cm; 23.694 GHz for (1,1) is RECALLED); NH3 was detected at centimetre wavelengths with the Effelsberg 100-m telescope (Hale-Bopp). No NH3 line near 1.15 mm (261 GHz) exists in these sources | Odin07 (Sect. 2, 6); Hat05 abstract; Cro16 Sect. 3 | READ |
| G1.4f2 | NH3 3.0 um | UNVERIFIABLE | iSHELL covers 2.8 to 5.3 um and detects NH3 ([R1312] abstract, SNIPPET); the band position was not read | | SNIPPET |
| G1.4g | HCN 88.6 GHz (J = 1-0) | CONFIRMED | three hyperfine components of J = 1-0 at 88.6 GHz (IRAM) | Cro16 Fig. 3 | READ |
| G1.4h | CN violet 388 nm | CONFIRMED | CN band at 3880 Angstrom | Coch15 Fig. 2 caption | READ |

### G1.5 Classify comets by chemistry (R9 suggestion)

| id | Claim | Verdict | What the literature gives | Source | Label |
|---|---|---|---|---|---|
| G1.5 | Classify by hyper-volatiles (CO/H2O, CH4/H2O) and D/H instead of dynamical class | PARTLY | Sound in spirit, not an established taxonomy. See Findings 2.5 | Coch15; Eis19; HP22; Rob24 DB; Alt15 | READ / COMPUTED |

## Figures for the generator

### 2.1 Species table (constants, not per comet)

Keep these on the species record (all PARAMETRISATION or table constants, not draws):

| Field | Source of the numbers | Note |
|---|---|---|
| `beta0_diss`, `beta0_ion` (1/s at 1 AU, quiet Sun, unshielded) | Hea17 Tables 18 and 19 times 37 700 (COMPUTED), cross-checked with HM15 via sbpy and with [R0934] via dVB13 | Values in the table in 2.2 |
| `lifetime_1au_s` (quiet) and `activity_factor` (about 1.7 to 2.4 for the active Sun; PARAMETRISATION from HM15 and Rag20) | sbpy; Rag20 | H2O 1.68, CO2 2.35 |
| `spectral_ratio_cool` and `spectral_ratio_hot` (rate for a 4000 K and a 10 000 K black body divided by the solar rate at equal 91.2 to 200 nm flux) | Hea17 Table 18 (COMPUTED) | See 2.2: the correction is species-specific and spans five orders of magnitude |
| photolysis `branching` table | Hea17 Table 22 for H2O and NH3; Table 1 for the others | For C2H6 no source read |
| `band` list (wavelength or frequency, transition) | section G1.4 | Display and detection only |

Scaling law for the generator (PARAMETRISATION assembled from READ parts): beta_i(r, star) = beta0_i x S_i(T_eff) x A_activity / r^2, with r in AU, S_i interpolated in log between the three template spectra of the species record (4000 K, Sun, 10 000 K) by the star's effective temperature, A_activity = 1 for a quiet star and the species factor above for an active one. Comet-coma length scale L = v_exp / beta. The r^-2 law is READ (dVB13); S_i and A_activity are not a published law, they are a reduction of Hea17's templates.

### 2.2 Computed quiet-Sun photolysis table (unshielded, 1 AU)

Hea17 text: the solar field is a "quiet period" spectrum (SOHO-SUMER plus UARS SOLSTICE); the tabulated rates are for the Draine-normalised field and the solar rates "should be scaled by 37 700" (flux 0.098 W per square metre between 91.2 and 200 nm). Output of `G1_heays_solar_lifetimes.py`:

| Species | beta_diss 1/s | beta_ion 1/s | tau total s | Rate ratio 4000 K / Sun | Rate ratio 10 000 K / Sun | Compare |
|---|---|---|---|---|---|---|
| H2O | 8.67e-6 | 1.58e-7 | 1.13e5 | 0.70 | 2.0 | HM15 8.3e4 |
| CO | 3.02e-7 | 2.45e-7 | 1.83e6 | 3.5e-4 | 2.2 | 1.54e6 (dVB13) |
| CO2 | 6.41e-7 | 1.96e-7 | 1.20e6 | 0.14 | 6.5 | HM15 4.95e5 |
| CH4 | 4.90e-6 | 3.51e-7 | 1.90e5 | 0.013 | 2.1 | Jac76 1.2e5 |
| C2H6 | 6.03e-6 | 8.67e-7 | 1.45e5 | 0.037 | 3.0 | none read |
| NH3 | 1.55e-4 | 3.62e-7 | 6.46e3 | 0.88 | 0.41 | Odin07 6.7e3 |
| HCN | 7.92e-6 | 2.83e-7 | 1.22e5 | 0.027 | 1.3 | dVB13 6.7e4 |
| CH3OH | 7.92e-6 | 6.03e-7 | 1.17e5 | 0.57 | 2.5 | dVB13 7.7e4 |
| H2CO | 1.92e-4 | 4.90e-7 | 5.19e3 | 0.78 | 0.19 | dVB13 5.0e3 |
| H2S | 2.30e-4 | 9.80e-7 | 4.33e3 | 0.69 | 0.36 | none read |
| N2 | 4.15e-7 | 4.90e-8 | 2.16e6 | 2.9e-5 | 1.0 | none read |

Reading it: Hea17-derived lifetimes agree with comet-community values within a factor 1.1 for NH3, H2CO, within 1.2 for CO, and are 1.4 to 2.4 times longer for H2O, HCN, CH3OH, CO2. Use the HM15/Crovisier values where they exist (G1.2 table) and the Hea17 values only for C2H6, H2S, N2, with an uncertainty factor 2. Shielding (self-shielding of CO and N2 in a dense inner coma) is not included in any number here. Spectral dependence is the larger problem: for a 4000 K star the CO and N2 photodissociation rates fall by 3.5e-4 and 2.9e-5 relative to the Sun at the same 91.2 to 200 nm flux, CH4 by 80, HCN by 40, while NH3 and H2CO barely change. A single "stellar UV" scalar is wrong by orders of magnitude for CO, N2, CH4, HCN, CO2 around cool stars.

### 2.3 Composition distributions for comets (PARAMETRISATION, from the database, r_h at most 2 AU, per-comet median)

Each species: log-normal in the ratio to water with the geometric mean and sigma below. Hard bounds are the observed extremes. These are coma ratios; 67P is the only bulk estimate.

| Species | Geo-mean % | sigma dex | Hard bounds % | P-type median | C-type median | Label |
|---|---|---|---|---|---|---|
| CO | 2.6 | 0.64 | 0.15 to 30 | 1.2 | 8.4 | COMPUTED |
| CO2 | 11 | 0.27 | 2 to 37 | 14.5 | 8 (N = 3 only) | COMPUTED |
| CH4 | 0.68 | 0.31 | 0.15 to 2.2 | 0.54 | 0.95 | COMPUTED |
| C2H6 | 0.37 | 0.29 | 0.09 to 0.9 | 0.30 | 0.62 | COMPUTED |
| NH3 | 0.55 | 0.22 | 0.12 to 0.9 | 0.66 | 0.48 | COMPUTED |
| HCN | 0.13 | 0.25 | 0.03 to 0.4 | 0.13 | 0.165 | COMPUTED |
| CH3OH | 1.9 | 0.32 | 0.2 to 4.5 | 1.95 | 2.4 | COMPUTED |
| H2CO | 0.25 | 0.34 | 0.09 to 1.1 | 0.19 | 0.30 | COMPUTED |
| H2S | 0.71 | 0.37 | 0.09 to 1.5 | 0.71 | 1.2 | COMPUTED |

Limits of the database: it contains the comets that also have a measured nucleus size, mixes several instruments, and 67P enters with in situ data. The type split is Periodic (P: JFC and Halley type) against Long-period or Hyperbolic (C); it does not separate JFC from Halley type or dynamically new from old. Dello [R1063] itself was not read.

### 2.4 Relations between species

| Relation | Value | Source | Label |
|---|---|---|---|
| (CO + CO2)/H2O | median 18 +- 4 percent for most comets (18 comets) | HP22 abstract, Table 3 | READ |
| Rank correlation CO against CO2 | rho = -0.86 (N = 8): CO and CO2 trade off | G1 script | COMPUTED |
| Rank correlation CO against CH4 | rho = 0.34 (N = 12); CO against C2H6 0.11; CH4 against C2H6 0.28; CO against HCN 0.48 (N = 14); CO against CH3OH -0.12 (N = 13). The hyper-volatiles CO and CH4 are weakly or not correlated, consistent with a search snippet that they vary most and are not correlated | G1 script; SNIPPET | COMPUTED |
| CO/H2O against nucleus radius, r_h at most 2 AU | log10(CO percent) = 0.08 + 0.91 log10(R in km), residual 0.52 dex, Pearson 0.66, N = 16 | G1 script on Rob24 DB. Rob24 reports the CO correlation as the only strongest one (three-sigma level) and stronger for ecliptic comets than nearly isotropic ones | COMPUTED / READ |

The size relation is confounded: Oort-type comets are larger and have higher CO. Treat as a candidate driver (Rob24's model: radiogenic heating redistributes hyper-volatiles in larger nuclei), not a law.

### 2.5 Taxonomy facts (answers to G1.5)

| Statement | Number | Source | Label |
|---|---|---|---|
| Optical (daughter-species photometry and spectroscopy: [R0049] 85 comets, [R2559] 130 comets, about 250 in total) | the majority are alike; 25 to 30 percent are different, mostly depleted in C2 and C3. Typical share quoted 75 to 91 percent and about 70 percent (Schleicher and Bair) | Coch15 Sect. 2.2; Eis19 | READ |
| Number of optical classes | 7 (typical about 70 percent; strongly C2+C3 depleted; weakly C2+C3 depleted; C2 only; C3 only; NH only; CN-low) | Coch15 | READ |
| Dynamical type of depleted comets | any type, "more common for Jupiter Family"; two thirds of depleted comets are JFC and one third long period; depleted fraction 37 percent of JFC against 18.5 percent of long-period | Coch15 | READ |
| Infrared (primary parent volatiles, [R3206], Dello [R1063]) | three groups: typical, organic-enriched, organic-severely-depleted; 8P/Tuttle is typical in the optical but peculiar in the IR | Coch15 Sect. 2.3 (Mumma and Charnley read second-hand) | 2ND-HAND |
| Radio ([R0936], about 30 comets) | no evidence for the three IR categories; no JFC/other difference; wide spread for most species; HCN narrow | Coch15 | 2ND-HAND |
| Dello [R1063] (30 comets, 8 volatiles) | JFCs depleted in volatiles relative to long-period (press summary only) | SNIPPET | SNIPPET |
| Dynamical class against composition (G1 test, permutation on per-comet log ratios, r_h at most 2 AU) | share of variance explained by P against C: CO 0.28 (p = 0.035; N = 9 against 7), H2S 0.25 (p = 0.16), CO2 0.14 (p = 0.17), HCN 0.08, NH3 0.06, CH3OH 0.05, H2CO 0.05, C2H6 0.04, CH4 0.00 (all p above 0.2) | G1 script | COMPUTED |
| CO against CO2 by family (HP22 Table 2) | JFC: 8 of 9 CO2-dominated; OCC: 5 CO2-dominated, 6 CO-dominated, 3 comparable; dynamically new Oort comets have more CO2 relative to CO than older ones | HP22 | READ |
| A chemical taxonomy from formation location | 15 comets fall on a "trail" tracing the CO iceline at about 30 AU early to 12 AU later; proposed as an alternative taxonomy and a formation-time classification | Eis19 abstract | READ |
| D/H in water | 67P (5.3 +- 0.7) e-4; 103P (1.61 +- 0.24) e-4; 45P below 2.0e-4 (upper limit); Halley about 3e-4; protosolar nebula (2.1 +- 0.5) e-5 | Alt15 text; Coch15 | READ |
| JWST era | no Solar System comet abundance table found. 3I/ATLAS (interstellar): CO2/H2O = 7.6 +- 0.3 and CO/H2O = 1.65 +- 0.09 at r_h = 3.32 AU (ratios, i.e. 760 and 165 percent), CH4 enriched | Cor25; Bel26 | READ |

So: chemistry-based classes exist, but the best-documented ones are on optical daughter species or on IR primary species (3 groups), not on CO and CH4; no source read defines hyper-volatile classes with fractions; D/H has too few values (3 to 4 JFC or Halley comets) to classify and does not follow dynamical class (the two JFCs differ by a factor 3).

## How the generator should compute or draw it

### 3.1 Per comet: store or draw

Already stored in the comet record: nucleus radius, dynamical `period_class`, interstellar flag. New stored quantities (PARAMETRISATION, seeded draws):

| Field | Unit | Draw | Per body or constant |
|---|---|---|---|
| `x_CO`, `x_CO2`, `x_CH4`, `x_C2H6`, `x_NH3`, `x_HCN`, `x_CH3OH`, `x_H2CO`, `x_H2S`, `x_N2` | percent of H2O (number ratio) | log-normal, geometric mean and sigma from 2.3, truncated at the hard bounds | per body |
| `carbon_chain_depleted` | flag | Bernoulli: p = 0.37 for JFC, 0.185 for long-period, 0.25 otherwise (from Coch15) | per body |
| `dh_water` | ratio | PARAMETRISATION only: uniform in log between 1.5e-4 and 5.5e-4 (the span of three measured JFC or Halley values); not tied to class | per body (optional) |
| coma `beta_i(r)` | 1/s | computed from the species constants and the star, never stored per comet | derived |

Constants of the species table (not per comet): `beta0_diss`, `beta0_ion`, quiet lifetime, `activity_factor`, spectral templates, branching, bands. A photolysis lifetime is a species constant at 1 AU scaled by r squared and by the star's spectrum, not a body property.

### 3.2 Route recommendation for comets (route 2)

Draw composition first, then derive bulk density, as part D already recommends. For the hyper-volatile block I recommend a three-step draw that keeps the observed anticorrelation and the 18 percent carbon budget:

1. S = (CO + CO2)/H2O: log-normal with median 18 percent and sigma 0.2 dex (GEN.629: the scatter among single comets from the [R3782] database, 0.197 dex for 9 comets within 2.5 AU and 0.138 within 2.0 AU, the widest kept; the first draft took the +-4 of 18 in HP22, 0.1 dex, which is the median's uncertainty and not the comets' spread).
2. f = CO/(CO + CO2): logit-normal with class median 0.08 for JFC and 0.5 for Oort type (from P and C medians 1.2/(1.2 + 14.5) and 8.4/(8.4 + 8), computed from different comet subsets, so rough), logit sigma about 1.5 so that 0.01 to 0.9 is covered (PARAMETRISATION).
3. All other species independent log-normals from 2.3.

Alternative (single-variable): draw log10(CO percent) = 0.08 + 0.91 log10(R) + N(0, 0.52), clipped to 0.15 to 30, and CO2 from S minus CO. Use either, not both. The dynamical-class split of CO (28 percent of variance) is the largest class effect found; all others are below 10 percent and not significant, so apart from CO (or the f step) and the depleted flag, composition draws should not depend on class.

A classification by CO/H2O and CH4/H2O alone is not supported as a taxonomy by what was read. Proposed instead: store the continuous fractions, derive a descriptive label at read time (for example CO/H2O below 1 percent "CO-poor", above 10 percent "CO-rich"; these thresholds are a PARAMETRISATION from the observed range) and keep dynamical class as the source of the `carbon_chain_depleted` probability.

### 3.3 Consistency tests a test could assert

1. For every comet: sum of x_i for the drawn volatiles stays below a cap implied by the ice mass fraction (part D, 40 to 50 percent ices by mass).
2. beta_i(r) = beta0_i x S_i x A / r^2: doubling r divides beta by 4 (exact).
3. Quiet-Sun lifetimes reproduce the table in G1.2 within the stated tolerance (2 percent for H2O, CO, CO2, HCN, CH3OH, H2CO; 10 percent for NH3).
4. The star at 4000 K gives CO beta below 1e-3 of the solar value; at 10 000 K above the solar value (Hea17 ratios).

## Calibration and test targets (Solar System only, never presets)

| Target | Value | Use |
|---|---|---|
| 67P bulk ratios (CO2 4.7, CO 3.1, CH4 0.34, C2H6 0.29, NH3 0.67, HCN 0.14, CH3OH 0.21, H2CO 0.32, H2S 1.10, N2 0.089 percent of H2O) | Rub19 Table 2 (READ) | the draw should be able to reach each value inside its 2-sigma |
| 103P/Hartley 2: CO 0.3 (0.15 to 0.45), CO2 about 20 | Dello [R0047, R1063] | CO-poor, CO2-rich corner |
| Hale-Bopp: CO 26, CH4 1.22, HCN 0.36, H2S 1.5 | Rob24 DB | CO-rich corner |
| Lifetimes | table G1.2 | species constants |

## Changes to parts A to F

1. Part D, 3.4.2 step 3: "other species as lognormal around the 67P table with sigma 0.5 dex" should change. The population geometric means differ from 67P by 10 times for CH3OH (1.9 against 0.21) and 2 times for CH4 (0.68 against 0.34), and the scatter is 0.22 to 0.37 dex (CO 0.64), not 0.5 dex. Use the table in 2.3 and keep 67P as a test target.
2. Part D, 3.4.2 step 6 (interstellar CO/H2O log-uniform on 0.01 to 1.6): the 3I/ATLAS value 1.65 +- 0.09 (at 3.3 AU, Cor25) sits at the top edge and the CO2/H2O value 7.6 is far beyond the 2 to 30 percent range of Solar System comets; widen the interstellar range or draw CO2 separately (one measurement; PARAMETRISATION).
3. Part F row 8: the lower CO bound 0.15 percent is real (Hartley 2, Coch15), so "not in what was read" is now closed; the CO figure in Part D (median 3.3, 0.3 to 26) is confirmed against HP22 and against Rob24.
4. Part D, comet record (`period_class` weights jupiter_family 50, halley_type 30, long_period 20): not checked here; unchanged.

## Still unreachable (Boss must research)

1. Dello [R1063]: the tier definitions, number of classes and fractions of the IR taxonomy. All of its per-comet values are in the Rob24 database but not the classes.
2. [R3206]: the three-group IR taxonomy and abundance ranges for about 20 species (only second-hand here).
3. [R0435], Phil. Trans. R. Soc. A 375, 20160252: the range table for 28 molecules; HAL record obspm-02268099 shows a link only to the refused publisher PDF.
4. [R2040] (PSS 106, 11) and [R2041] (Ap&SS 195, 1): quiet and active photo rates for CH4, C2H6, NH3, H2S, N2 and the full branching tables; [R0934]. The PHIDRATES site (phidrates.space.swri.edu) resets the connection.
5. C2H6 photolysis branching at Lyman-alpha and the NH3 3.0 um band position in a comet paper [R1312].
6. [R0410] (the note's reference): no single paper was identified; [R0408] is read for HNCO, HCOOH, CO below 2 percent and H2S 0.15 percent in C/2021 A1 and C/2022 E3 only.
7. Solar System JWST comet papers (for example abundance tables for 12P/Pons-Brooks or C/2017 K2): searches returned only 3I/ATLAS papers and ground-based C/2017 K2.

## Counts

31 verdict rows (G1.1 six, G1.2 seven, G1.3 seven, G1.4 ten, G1.5 one). The five added species rows are not claims and are not counted.

| Verdict | Count | Rows |
|---|---|---|
| CONFIRMED | 17 | G1.1c, G1.1d; G1.2a, b, c, g; G1.3a, b, d, f, g; G1.4a1, b, c, d, g, h |
| PARTLY | 9 | G1.1a, b, e, f; G1.2d; G1.3c; G1.4a2, e; G1.5 |
| CONTRADICTED | 2 | G1.2f (NH3 lifetime); G1.4f1 (NH3 1.15 mm) |
| UNVERIFIABLE | 3 | G1.2e (C2H6 lifetime); G1.3e (C2H6 products); G1.4f2 (NH3 3.0 um) |

Approximate label counts for the figures in this file: READ about 55 items, 2ND-HAND about 8, COMPUTED about 60 (table cells), SNIPPET 4, RECALLED 4 (OH 18 cm frequencies, NH3 23.694 GHz, CH + H2 + H channel, C2H6 products), NOT REACHED 7 sources (list above).
