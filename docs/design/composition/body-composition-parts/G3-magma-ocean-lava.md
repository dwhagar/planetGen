# G3: magma-ocean redox, carbon solubility, outgassed split and lava-world vapour (check of Boss's second note, section 2)

Scope: claims G3.1 to G3.7 of Boss's note (messages cmsg_01RuNTFUNB92A8cm5CbojQAiBH5ozwSDue3azzQUisow2D), plus what the generator must store and a parametrisation of the mantle delta-IW draw. Read-only research. Scripts, extracted tables and outputs are in `/mnt/project-files/design-drafts/parts/G-files/` with prefix `G3_`. Existing work reused: F items 5 and 6 (Sossi et al. 2020, PMC7688334, re-read here, all figures below are from that re-read).

Labels as in `_brief.md`: READ, 2ND-HAND, COMPUTED, SNIPPET, RECALLED, NOT REACHED. "Code" below means open-source code that implements a published equation (original paper not reached), so the constants are 2ND-HAND.

Notation: dIW = log10 fO2 minus log10 fO2 of the iron-wustite buffer at the same T and P. "H21" = Hirschmann 2021 IW calibration (see Figures, F1). QFM = FMQ = quartz-fayalite-magnetite.

## Verdicts

| id | claim as written (short) | verdict | the right figure (with unit) | source | label |
|---|---|---|---|---|---|
| G3.1a | Oxidised mantle (dIW >= +2): carbon dissolves as carbonate CO3^2- | CONFIRMED | Carbonate is the dissolved carbon species above about dIW +1 (Gaillard 2022: "C dissolves as CO (fO2 < IW) and as CO3^2- (fO2 > IW+1)"); Sossi 2020: at dIW +0.5 H, C, N dissolve as CO3^2-, H2O/OH, N2. Both include dIW >= +2 | Gaillard et al. 2022 EPSL 577:117255 (HAL insu-03383209, results section 5.3); Sossi et al. 2020 Sci. Adv. (PMC7688334), text after Eq. 4 | READ |
| G3.1b | dIW >= +2 "corresponding to or exceeding the QFM buffer" | CONTRADICTED | QFM is above dIW +2 by a wide margin: QFM minus IW = 3.63 / 3.36 / 3.23 / 3.12 log units at 1500 / 1800 / 2000 / 2200 K and 1 bar. dIW +2 is about QFM -1.1 to -1.6. Sossi 2020 and Deng 2020 themselves put Earth's present upper mantle at about IW +3.5 | Table in F1 below (COMPUTED from Frost 1991 constants, 2ND-HAND via PetroBuffer 0.0.1, and H21 IW); Sossi 2020 text "spinel-bearing peridotites ... ~dIW +3.5"; Deng et al. 2020 Nat. Commun. 11:2007 | COMPUTED, READ |
| G3.1c | QFM versus IW offset "about 3.5 to 4" at magma-ocean temperature (task note) | PARTLY | 3.1 to 3.6 at 1 bar over 2200 to 1500 K. It reaches 3.4 to 4.0 only at 1 GPa and 3.9 to 4.7 at 3 GPa, because the solid-buffer volume term of Frost (0.110 vs 0.055 per bar/T) adds about 0.5 to 1.5 log units per GPa. At magma-ocean T the QFM is metastable (fayalite melts at 1477 K, H21); using it above that is standard practice but the high-P number is an extrapolation. Doyle 2019 rounds the conversion to "dIW = dQFM - 4" | F1 table; Hirschmann 2021 text on metastable buffers; Doyle et al. 2019 (arXiv 1910.12989) | COMPUTED, READ |
| G3.1d | Buffer equations: Frost 1991 / Hirschmann 2021 | PARTLY (source detail) | A polynomial fits QFM free energy only for 900 to 1420 K (polynomial given in the paper, scan OCR garbled; not used). Frost 1991 form log10 fO2 = A/T + B + C (P-1)/T, P in bar, constants in F1. H21 is the IW source for 1000 to 3000 K, 1e-4 to 100 GPa; there is no matching QFM calibration at magma-ocean T | Frost constants 2ND-HAND (PetroBuffer code); H21 Table 1 READ | READ, 2ND-HAND |
| G3.2a | Reduced (dIW <= -1.5), graphite-saturated high-P experiments: carbon is molecular CO, carbonyl, dissolved CH4 when water is present | PARTLY | Observed (FTIR, Raman, SIMS): in anhydrous reduced glass CO is the main C species (Armstrong 2015, <120 ppm CO); carbonate "expected only in very small quantities"; CH4 appears when the glass has about 400 ppm H2O or more and stays stable from dIW -0.4 to -3 (Dalou 2019 summary of Li 2015, Armstrong 2015, Ardia 2013, Kadik). Below dIW -3.0 C is mostly C-H species; between -3.0 and -1.5 C is C=O (molecular CO or Fe-carbonyl) or C-H (Grewal 2020, SIMS+FTIR+Raman, IW-4.3 to IW-0.8, 1-7 GPa, 1500-2200 C). Not all experiments agree: Dasgupta-group abstracts report carbonate dominant down to about IW-0.4 to IW-2 in graphite-saturated basalt (SNIPPET) | Dalou et al. 2019 (Raman, IW-2.1 to IW-0.4; par.nsf.gov purl 10137270), sections 4.1 and 4.2; Grewal et al. 2020 abstract (cosmochemistry-papers.com copy); Duncan and Dasgupta AGU 2015 | READ, SNIPPET |
| G3.2b | "as fO2 drops from dIW -1.4 down to -4.7 carbonate solubility becomes negligible while CH4 and non-carbonate dissolution increase with melt hydration" | UNVERIFIABLE (range) / direction supported | No carbon study found with the range -1.4 to -4.7. The only -4.7 found is the lower end of the N2-solubility experiments (fO2 IW+4.9 to IW-4.7), which looks like the source of the figure. Carbon experiments: Armstrong 2015 IW-3.65 to +1.46 (1.2 GPa, 1400 C); Dalou 2019 IW-2.1 to -0.4; Grewal 2020 IW-4.3 to -0.8; Hirschmann 2016 CO/CH4 law calibrated IW-3.6 to +2.3 (350-32000 bar, 1400-1800 C); Ardia 2013 CH4 at 0.7-3 GPa, 1400-1450 C. Direction: Ardia 2013 model has CH4 solubility rising with H content and P and falling fO2 (Dalou summary); Gaillard 2022 has C-solubility in the melt lowest at intermediate fO2 (2 to 15 ppm C between IW-4 and IW) | atmodeller 1.0.2 code docstrings (Armstrong 2015 range); Maurice et al. 2024 (arXiv 2405.09284) Table A.1; Gaillard 2022 text | 2ND-HAND, READ |
| G3.2c | The experiments are "Hirschmann 2012; Gaillard 2014, 2021" | PARTLY | These are models, not experiments. Gaillard and Scaillet 2014 (EPSL 403:307) is a thermodynamic degassing framework for volcanic gases (fO2 FMQ-4 to FMQ). Gaillard 2022 (the "2021" paper, EPSL 577:117255) extrapolates solubility laws to magma-ocean conditions. Hirschmann 2012 (EPSL 341-344:48) is a magma-ocean/atmosphere model (original not reached). The experiments are Ardia 2013, Armstrong 2015, Dalou 2019, Grewal 2020, Yoshioka 2019, Duncan 2017 (IW+1.6, 1 GPa, 1600 C, 0.4 to 0.9 wt% CO2 across basalt-peridotite join, SNIPPET) | HAL texts of Gaillard and Scaillet 2014 and Gaillard 2022; snippets | READ, SNIPPET |
| G3.2d | Reaction C + 2 H2O(melt) = CH4(melt) + O2 | PARTLY | Balanced as a net reaction. The melt-species form used with Raman data is C(graphite) + H2(melt) + 2 OH-(melt) = CH4(melt) + O2-(melt) (Dalou 2019 Eq. 3). Mass action then gives X_CH4 rising with H content and falling fO2, as Ardia 2013 found. Gas-phase check: CO + 3 H2 = CH4 + H2O has log10 K = -5.35 (1500 K), -7.29 (2000 K), -7.75 (2173 K) in bar^-2, so CH4/CO < 0.1 at 2000 K and <= 100 bar even at dIW -4.7 (ideal gas) | Dalou 2019 section 4.1; G3_out_ch4_equilibrium.txt | READ, COMPUTED |
| G3.3a | Oxidised: C_melt proportional to fCO2 times NBO/T; valid to about 0.7 GPa for NBO/T > 2 ultramafic melts | PARTLY | Physical part: mass action, dissolved CO3^2- proportional to fCO2 with a volume term: X_CO3 = X0 fCO2 exp(-DeltaV (P-1)/(R T)), DeltaV = 23 cm3/mol, T0 = 1473 K (Dixon 1997, via code). Linear in P is observed to about 0.7 GPa for basalt and rhyolite with a common slope of about 0.57 ppm CO2 per bar (Ni and Keppler Fig. 1a), and "may apply to 4 GPa". The composition term is exponential in NBO/O, not linear in NBO/T: Iacono-Marziano 2012 ln C = sum(d_i x_i) + 1.0 ln fCO2 + 15.8 NBO/O - 5.3 + 0.14 P/T (anhydrous set). No source found for "valid to 0.7 GPa for NBO/T > 2". NBO/T alone is not enough (Ni and Keppler: Mg and Ca roles differ, Duncan and Dasgupta AGU abstract) | VESIcal 1.2.12 dixon.py and iaconomarziano.py (2ND-HAND) | READ, 2ND-HAND |
| G3.3b | Moderately reduced: C_melt proportional to P times exp(-dIW) | PARTLY (schematic) | Not a measured law. Mass action gives, for dissolved CO: C_CO proportional to fCO^0.876 (Armstrong, 1.2 GPa) to fCO^0.8 (Yoshioka, 3 GPa); and if the gas is CO2-dominated, fCO = fCO2 / (K 10^(0.5 dIW)), so C_CO falls as 10^(-0.5 dIW) = exp(-1.15 dIW), not exp(-dIW). For carbonate the sign is the opposite (proportional to fCO2, which rises as 10^(+0.5 dIW) at fixed fCO). The sign of the dIW term depends on which species and which fixed reference | derived here from ratio table F3; Armstrong and Yoshioka forms from atmodeller code | COMPUTED, 2ND-HAND |
| G3.3c | Highly reduced: C_melt proportional to fH2^2 times X_H2O | PARTLY (schematic) | fH2^2 follows from C + 2 H2 = CH4 with Henry's law for CH4 (CH4 solubility in the Henry form is linear in fCH4: Lichtenberg 2021 beta = 1.0). The X_H2O factor reflects the observed rise of CH4 with melt H (Dalou, Ardia) and the melt reaction in G3.2d (which would give a product fH2 times fH2O). No published closed form with this exact product found | Lichtenberg et al. 2021 Table 1; Dalou 2019 Eq. 3 | READ |
| G3.4a | Oxidised: outgassed gas is CO2, H2O, N2 | PARTLY | Deng 2020 at dIW +2 (1 bar, 1800 K, H/C mass 0.5, no solubility): H2O about 70 mol%, CO2 about 15 mol%. Gaillard 2022: > IW+2 gives a dry, C-N-S-rich atmosphere (H2O stays dissolved; SO2 appears), and Sossi 2020: about 99 % of the BSE water is dissolved at 2173 K. So H2O in the gas depends on solubility and inventory, not only on dIW | Deng 2020 Nat. Commun. section "redox states of the MOs"; Gaillard 2022 abstract; Sossi 2020 | READ |
| G3.4b | Moderately reduced: CO, H2O, H2, CO2 | CONFIRMED | Sossi 2020 (dIW +0.5, 2173 K, H/C = 0.22, H/N = 5.8): about 100 bar CO, 16 bar CO2, 7 bar H2, 7 bar H2O, 2.5 bar N2, total about 140 bar. Deng 2020 Mars (dIW -0.3): H2O 40, H2 40, CO 15, CO2 5 mol%. Maurice et al. 2024: lava-world magma oceans vary "from CO to CO2-N2-SO2 with increasing fO2" | Sossi text after "Given estimates of H, C, and N"; Deng 2020; Maurice 2024 abstract | READ |
| G3.4c | Highly reduced (<= -2): H2, CO, CH4, H2S | PARTLY | H2 and CO yes: Deng 2020 Moon (dIW -2): H2 > 70, CO about 20, H2O 10 mol%. Gaillard 2022 (< IW-2): H-C(+N)-rich atmosphere, graphite stable below about IW-3. CH4 is minor at magma-ocean temperature (Sossi: X_CH4 about 1e-3 at dIW +0.5; computed CH4/CO < 0.1 for T >= 2000 K, P <= 100 bar even at dIW -4.7) and appears on cooling or at high H/C (Sossi Fig. 4: CH4 dominates for log(H/C) > 0.7, frozen in below about 900 to 1000 K). H2S is not supported: Gaillard 2022 says S stays in the melt at low fO2 and "only becomes volatile at the uppermost fO2" (as SO2) | Deng 2020; Gaillard 2022 section 5.1; Sossi 2020; G3_out_ch4_equilibrium.txt | READ, COMPUTED |
| G3.4d | Sossi's CO-rich gas at dIW +0.5 lies in the "moderately reduced" row, not the oxidised one | CONFIRMED | Yes. Computed crossover (ideal gas, NASA Glenn fits, H21 IW): CO2 = CO at dIW +1.00 (1500 K), +1.20 (1800 K), +1.37 (2173 K); H2O = H2 at +0.17, +0.03, -0.06. So Boss's "oxidised: dIW >= +2" threshold is right for CO2 > CO at 2173 K (CO2/CO = 2.1 at +2). On Sossi's own scale the crossover is +1.9 | F3 below | COMPUTED |
| G3.4e | Split is set by redox | PARTLY | For magma-ocean degassing redox is first order (Gaillard 2022, Sossi 2020, Maurice 2024). For volcanic degassing Gaillard and Scaillet 2014 find surface pressure first order: high P gives N2-CO2-rich dry gas, about 1 bar gives H2O-dominated gas, low P gives S-rich gas, and different source redox states "converge toward a unique planetary fO2" (about FMQ-2 +/- 0.5 at 10 to 0.01 bar) | Gaillard and Scaillet 2014 abstract and Fig. 2 text | READ |
| G3.4f | Compare with Schaefer and Fegley 2010, Bower 2019, Lichtenberg 2021, Nicklas 2019, Gaillard 2022 (arXiv versions) | UNVERIFIABLE (partly) | Bower 2019 (arXiv 1904.08300) and Lichtenberg 2021 (arXiv 2101.10991) do not speciate by fO2: Bower tracks CO2 and H2O only; Lichtenberg runs end-member single-gas atmospheres (H2, H2O, CO2, CH4, CO, O2, N2) and a Henry table (F5). Schaefer and Fegley 2010 (Icarus 208:438): NOT REACHED (no arXiv version found). Nicklas 2019: not identified; Gaillard 2022 has no arXiv version; the HAL open copy was read. Schaefer and Fegley 2017 (ApJ 843:120, par.nsf.gov purl 10038075) is the readable redox-of-outgassing source: CI chondritic material about QFM -1, EH chondritic about IW -1 | as listed | READ, NOT REACHED |
| G3.5a | LavAtmos is arXiv:2210.10463 (Zilinskas et al.) | PARTLY | id and title correct: "LavAtmos: An open source chemical equilibrium vaporisation code for lava worlds". Authors are van Buchem, Miguel, Zilinskas and van Westrenen (Meteoritics and Planetary Science 58:1149, 2023, doi 10.1111/maps.13994); Zilinskas is third author. Zilinskas et al. 2022 is a different paper ("Observability of evaporating lava worlds", arXiv 2202.04759, A&A 661:A126). Sequel: LavAtmos 2.0 with volatiles, arXiv 2408.10863. Oxides: 9 (SiO2, MgO, Al2O3, TiO2, Fe2O3, FeO, CaO, Na2O, K2O; no Cr2O3), 31 reactions, MELTS activities via ThermoEngine, GPL-3.0 | LavAtmos paper text, Table 1 | READ |
| G3.5b | VapoRock oxide set SiO2-MgO-FeO-Fe2O3-Al2O3-CaO-Na2O-K2O-TiO2-Cr2O3 | CONFIRMED | Exactly this set (text Section 2.2). Authors Wolf, Jaggi, Sossi, Bower; arXiv 2208.09582;. 34 gas species. Liquid: MELTS through ThermoEngine, liquid components Table 1. Gas thermodata: NIST-JANAF or Lamoreaux and Hildenbrand 1984. species). Licence of the code: GNU AGPL v3 (LICENSE and README on gitlab.com/ENKI-portal/vaporock, branch master); Zenodo record 4594226 is CC-BY-4.0 for the record only | VapoRock paper Table 1; gitlab raw LICENSE; Zenodo API | READ |
| G3.5c | VapoRock "is a Gibbs minimisation tool" | PARTLY | It is not a global Gibbs minimiser. Each vapour species gets a mass-action partial pressure from equality of chemical potentials with the melt oxide components (paper Eqs. 1 to 3): log P_i = [sum nu_ij mu_j + nu_O2 mu0_O2 - mu0_i]/(RT ln10) + nu_O2 log fO2, in the "vapour-limited regime" (melt composition unchanged), ideal gas, fO2 an input. Melt activities from MELTS. Runs without H, C, N, S | VapoRock paper Section 2.3 | READ |
| G3.6a | 1500-2000 K "sub-liquidus": vapour dominated by Na and K | PARTLY | Na is the dominant gas for BSE across 1500 to about 2700 K; second and third are O2 and O (Schaefer and Fegley 2009; Ito 2015; Miguel 2011 Type I Na, O2, O, Fe). K is minor: Ito (BSE at 2000 K, 1e-5 bar) molar fractions Na 0.721, O2 0.189, K 0.0407, O 0.0233, SiO 0.00731, Fe 0.001; VapoRock 1900 K K/Na = 0.11 to 0.15 for tholeiite, komatiite (up to 1.5 for dunite). Older MAGMA runs gave K/Na 100 to 600 times lower (known K error, Wolf 2023 Table 3). Wording: all of these treat a melt; BSE-like peridotite liquidus exceeds 2000 K at 1 bar (Sossi), so 1500-2000 K is below liquidus and the melt is metastable or partial | sf09 text; Ito Fig. 3 caption; VapoRock Table 3; Miguel 2011 section 3.1.1 | READ |
| G3.6b | T > 2200 K: SiO2 -> SiO + 1/2 O2, MgO -> Mg + 1/2 O2, FeO -> Fe + 1/2 O2, TiO2 -> TiO + 1/2 O2 | PARTLY | First three are the right reaction types (LavAtmos Table 1 reactions 5, 20, 18; Sossi/Wolf: FeO(l) = Fe(g) + 1/2 O2; "alkalis evaporate via n = 1 reactions, A+ O0.5 = A0 + 1/4 O2"). The main Mg gas is Mg(g), Fe(g), Ca(g), K(g), AlO(g), CrO2(g). Ti is different: the main Ti gas for BSE is TiO2(g), vaporising without O2 exchange; TiO(g) only becomes significant at lower fO2 (proportional to fO2^-0.5). Ito 2015: TiO mole fraction about 1e-12 at 1500 K and 1e-6 at 3000 K. Reactions begin well below 2200 K; there is no threshold | VapoRock Section 4.3; Ito section 3 | READ |
| G3.6c | Vapour shifts from alkali gas to heavy oxide gas above 2200 K (Na stops dominating) | CONTRADICTED | Na stops dominating at about 2700 K for BSE (VapoRock; SiO takes over), about 2800 K (Ito, BSE and MORB), 2600 K (bulk crust), 2500 K (upper crust); VapoRock crossover 2500 K for Bishop Tuff and above 3000 K for dunite; Miguel komatiite: SiO exceeds Na only above 3170 K (MAGMA). The crossover rises with fO2 and falls with melt SiO2 and alkali content | VapoRock Section 4.2, 4.3; Ito section 2.1; Miguel section 3.1.1 | READ |
| G3.6d | Mean molecular weight of heavy vapour "about 40 to 60 g/mol" | CONTRADICTED | Mean molar mass of the vapour over BSE: 25.4 g/mol at 2000 K, 1e-5 bar and 33.2 g/mol at 3000 K, 1e-2 bar (from Ito's seven listed species). Elemental atoms only at 1900 K (VapoRock Table 3): 25 to 28 g/mol for tholeiite, alkali basalt, komatiite, Bishop Tuff and 33 for dunite. Kite 2016 uses 34.15 g/mol (mean of SiO and Mg) for fractionated vapour. 40 to 60 would need SiO2 (60) and FeO (72) to dominate; no source reports that | G3_out_mean_molar_mass.txt; Kite 2016 Table 1 | COMPUTED, READ |
| G3.6e | Alkali gas M about 23 to 39 g/mol | CONFIRMED | Na 22.99, K 39.10 g/mol (trivial) | standard atomic weights | COMPUTED |
| G3.7 | Line lists TiO, VO, FeH, CrH, Na, K govern the radiative structure | PARTLY | Line-list sources are standard (Kurucz, VALD3, ExoMol, HITRAN). Zilinskas 2022 uses Kurucz for atoms (Al, Ca, Fe, K, Mg, O, Si, Ti), VALD3 for Na, EBJT for SiO, ExoMol Toto for TiO, LiTY for MgO, OYT3 for SiO2, ATP for AlO, VBATHY for CaO, HITRAN for O2. VO, FeH, CrH are not in that set (they are hot-Jupiter and M-dwarf opacities, RECALLED). For BSE the inversion is made by SiO and Na/K; TiO matters for Ti-rich crustal melts and after fractional vaporisation (Zilinskas 2022; Ito: TiO opacity always below Na for BSE) | Zilinskas 2022 Table A.1; Ito section 3 | READ |

## Figures for the generator

### F1 Buffers (stored fO2 reference)

IW, recommended: Hirschmann 2021 (GCA 313:74), valid 1000 to 3000 K, 1e-4 to 100 GPa, r.m.s. 0.0065 log units. log10 fO2 = a + b T + c T ln T + d/T (fcc/bcc iron), with each parameter m = m0 + m1 P + m2 P^2 + m3 P^3 + m4 P^0.5, T in K, P in GPa; hcp-iron set (e,f,g,h) when P > -18.64 + 0.04359 T - 5.069e-6 T^2 (GPa). Magnitudes READ on par.nsf.gov purl 10335619 (Table 1; the text extraction lost the minus signs); signs as coded in atmodeller 1.0.2 `thermodata/_redox_buffers.py` (2ND-HAND, GPL-3.0). Coefficients (m0, m1, m2, m3, m4):

| set | m0 | m1 | m2 | m3 | m4 |
|---|---|---|---|---|---|
| a | 6.844864 | 1.175691e-1 | 1.143873e-3 | 0 | 0 |
| b | 5.791364e-4 | -2.891434e-4 | -2.737171e-7 | 0 | 0 |
| c | -7.971469e-5 | 3.198005e-5 | 0 | 1.059554e-10 | 2.014461e-7 |
| d | -2.769002e4 | 5.285977e2 | -2.919275 | 0 | 0 |
| e | 8.463095 | -3.000307e-3 | 7.213445e-5 | 0 | 0 |
| f | 1.148738e-3 | -9.352312e-5 | 5.161592e-7 | 0 | 0 |
| g | -7.448624e-4 | -6.329325e-6 | 0 | -1.407339e-10 | 1.830014e-4 |
| h | -2.782082e4 | 5.285977e2 | -8.473231e-1 | 0 | 0 |

Simple alternative (Frost 1991 RiMG 25, form log10 fO2 = A/T + B + C (P-1)/T, P in bar; constants 2ND-HAND from PetroBuffer 0.0.1 `buffers.py`): IW (-27489.0, 6.702, 0.055); FMQ above 846 K (-25096.3, 8.735, 0.110); NNO (-24930.0, 9.360, 0.046). At 1 bar it agrees with H21 within 0.003 to 0.024 log units for 1500 to 2200 K (COMPUTED). The pressure term is a solid-volume correction; it is not reliable at GPa for melts.

QFM minus IW (COMPUTED, `G3_buffers_and_gas_ratios.py`; Frost both buffers; value with H21 IW in brackets where different by more than 0.02):

| T (K) | 1 bar | 1 GPa | 3 GPa | 5 GPa |
|---|---|---|---|---|
| 1500 | 3.63 | 4.00 (3.97) | 4.73 (4.67) | 5.46 (5.38) |
| 1800 | 3.36 | 3.67 (3.65) | 4.28 (4.24) | 4.89 (4.83) |
| 2000 | 3.23 | 3.50 | 4.05 (4.02) | 4.60 (4.56) |
| 2200 | 3.12 | 3.37 | 3.87 | 4.37 (4.34) |

NNO minus IW is 3.8 to 4.4 at 1 bar (1500 to 2200 K).

Caution on the IW scale: Sossi's gas ratios at dIW +0.5 and 2173 K (f_H2O/f_H2 = 1.05, f_CO2/f_CO = 0.21) correspond to dIW 0.0 on the H21 scale when computed from NASA Glenn gas data (f_H2O/f_H2 = 1.08, f_CO2/f_CO = 0.208 at dIW 0.0; 1.92 and 0.369 at +0.5). So Sossi's IW sits about 0.5 log units lower in fO2 than H21's at 2173 K. H21 itself says JANAF-based IW is 0.4 to 0.7 too reducing at 2000 K, and Maurice et al. 2024 found a 0.5 offset between the O'Neill and Eggins IW (used by Bower 2022) and the IW they use. Any generator that stores dIW must fix the reference (H21 Table 1) and must not mix values quoted in the literature without converting. This applies to F item 5 and work item F2 in `F-boss-note-check.md` (Earth's +0.5 is a calibration target on Sossi's scale, about 0.0 on H21's).

### F2 Melt iron redox to dIW

Sossi 2020 Eq. 4 (peridotite liquid, 1 bar, 2173 K, 11 experiments, dIW -1.5 to +6.4): log10(X_FeO1.5 / X_FeO) = (0.252 +/- 0.015) dIW - (1.530 +/- 0.053). Fe3+/Fe-total measured 0.014 +/- 0.009 at dIW -1.5 to 0.440 +/- 0.008 at +6.4. Inverse (COMPUTED, `G3_fe3_to_diw.py`): dIW = (log10(f/(1-f)) + 1.530)/0.252.

| Fe3+/Fe_total | 0.014 | 0.02 | 0.037 | 0.05 | 0.10 | 0.20 | 0.44 |
|---|---|---|---|---|---|---|---|
| dIW (Sossi scale) | -1.26 | -0.64 | +0.45 | +1.00 | +2.28 | +3.68 | +5.66 |

So a stored mantle Fe3+/Fe-total (a part-B quantity) gives dIW directly, with the stoichiometric slope (about 1/4) and one decade of Fe3+/Fe2+ per 4 dIW units. Magma-ocean Fe3+/Fe-total from metal-silicate equilibration: Earth 0.034 to 0.10 (28 to 53 GPa, 2300 to 4100 K), Mars 0.026 to 0.038 (Hirschmann 2022, GCA 328:221, READ abstract); these map to dIW +0.3 to +2.3 and -0.2 to +0.5 (COMPUTED, indicative, pressure outside calibration).

### F3 Gas speciation versus dIW (ideal gas, COMPUTED, NASA Glenn 9-coefficient fits, 2ND-HAND from atmodeller 1.0.2 data file, subset in `G3_nasa_glenn_subset.csv`)

log10(f_CO2/f_CO) = logK1(T) + 0.5 log10 fO2; log10(f_H2O/f_H2) = logK2(T) + 0.5 log10 fO2; log10 fO2 = IW(T) + dIW.

| T (K) | logK1 (CO + 1/2 O2 = CO2) | logK2 (H2 + 1/2 O2 = H2O) | dIW where CO2 = CO | dIW where H2O = H2 |
|---|---|---|---|---|
| 1500 | 5.311 | 5.724 | +1.00 | +0.17 |
| 1800 | 3.687 | 4.269 | +1.20 | +0.03 |
| 2173 | 2.303 | 3.018 | +1.37 | -0.06 |
| 2200 | 2.221 | 2.943 | +1.38 | -0.07 |

CO2/CO and H2O/H2 at 2173 K: dIW -4.7: 9.3e-4, 4.8e-3; -3: 6.6e-3, 0.034; -2: 0.021, 0.108; -1.4: 0.041, 0.215; 0: 0.208, 1.08; +0.5: 0.369, 1.92; +2: 2.08, 10.8; +3: 6.6, 34. (At 1500 K the same dIW values give 1.4e-3, 3.7e-3 at -4.7; 0.063, 0.163 at -1.4; 3.2, 8.2 at +2.) Check against Deng 2020 at 1800 K: dIW -2 gives H2O/H2 about 0.1 (Deng: 10/70 = 0.14), dIW +2 gives CO2/H2O in the gas of 15/70 (Deng) with H/C = 0.5 mass. Gas CH4 is small at these T and P <= 100 bar (G3.4c).

### F4 Mass action rule for a generator (2 reservoirs)

Per volatile element e: bulk mass = melt mass x dissolved fraction + atmosphere mass; atmosphere mass of species s = 4 pi R^2 p_s mu_s / (g mu_atm) (Lichtenberg Eq. 10; Maurice Appendix C). Gas speciation: reactions CO/CO2, H2/H2O, CO + 3H2 = CH4 + H2O from F3 at the surface T, with fO2 = IW(T,P) + stored dIW. Dissolved amounts: Henry-type laws F5 per species. Solve for p_s with mass balance (Maurice 2024 and Sossi 2020 do this; Maurice reports surface dIW varies by about 3 log units over mass and magma-ocean-fraction space for each fixed dIW at the magma-ocean bottom).

### F5 Solubility laws the generator can use (form: fugacity, T, P to dissolved mass fraction)

Evaluations in `G3_out_solubility_eval.txt`. All masses are of the named species. Calibration ranges are what the sources state.

| species | form | constants | calibration | source | label |
|---|---|---|---|---|---|
| H2O | Henry power law X = k f^0.5 (ppm, f in bar) | k = 647 (peridotite melt) | 2173 K, 1 bar, dIW -1.9 to +6.0 | as coded in atmodeller | 2ND-HAND |
| H2O | X = k f^0.5 (ppm, f in bar) | k = 965 (MORB) | 1200 C, 200-717 bar | refit by Sossi, atmodeller | 2ND-HAND |
| CO2 (carbonate) | X_CO3 = X0 f_CO2 exp(-DeltaV (P-1)/(R T0)); wt% CO2 = 4400 X/(36.594 - 44 X) | X0 = 3.817e-7 (SiO2 > 48.9 wt%) or 8.697e-6 - 1.697e-7 SiO2 (wt%) below; DeltaV = 23 cm3/mol; R = 83.15 cm3 bar/(mol K); T0 = 1473.15 K | basalt, 1200 C, below about 2 kbar; 45 ppm CO2 at f = P = 100 bar, 381 ppm at 1000 bar (COMPUTED) | Dixon 1997 Am. Min. 82:368 (Eqs. 1, 8), via VESIcal 1.2.12 | 2ND-HAND |
| CO2 (carbonate) | ln[CO3 ppm] = sum(d_i x_i) + a ln f_CO2 + b NBO/O + B + C P/T (P bar, T K) | anhydrous set: d = (2.3, 3.8, -16.3, 20.1) for x = (H2O, Al2O3/(CaO+K2O+Na2O), FeO+2Fe2O3+MgO, Na2O+K2O), a = 1.0, b = 15.8, B = -5.3, C = 0.14 | 2012 calibration (range not read) | Iacono-Marziano et al. 2012, via VESIcal 1.2.12 | 2ND-HAND |
| CO2 | S_CO2 = 1.1 +/- 0.15 x P_CO2 (komatiite estimate) | units not printed in Gaillard 2022 Eq. 7 (ppm per bar assumed) | extrapolated to magma-ocean T | Gaillard 2022 supplement, from Iacono-Marziano 2012 | READ (units unverified) |
| CO2 (generic) | Henry fit X (ppmw) = alpha p^(1/beta), p in Pa | alpha = 1.937e-9, beta = 0.714; 310 ppmw at 1000 bar, 12 ppmw at 100 bar (COMPUTED); within a factor 4 of Dixon at 100 bar, 20 percent at 1000 bar | compilation fit to Mysen, Stolper and Holloway, Pan, Blank, Dixon | Lichtenberg 2021 Table 1, Eq. 9 | READ |
| CO | log10 CO (ppm) = -0.738 + 0.876 log10 f_CO - 5.44e-5 P (bar) | as written; 68 ppm at f = P = 1 kbar, 130 ppm at f = 1e4 bar, P = 12 kbar (COMPUTED) | 1.2 GPa, 1400 C, dIW -3.65 to +1.46, Martian and terrestrial basalt; pressure term fitted with data at 1 to 1.2 GPa; do not extrapolate to 3 GPa (gives 14 ppm at 30 kbar) | Armstrong et al. 2015 GCA 171:283 Eq. 10, via atmodeller | 2ND-HAND |
| CO | log10 CO (wt%) = -5.20 + 0.8 log10 f_CO (MORB); -4.08 + 0.52 log10 f_CO (rhyolite) | as written; 100 ppm at f_CO = 1e4 bar (MORB) | 3 GPa, 1500 C, graphite plus CO-CO2 fluid | Yoshioka et al. 2019 GCA, via atmodeller | 2ND-HAND |
| CO | Henry fit alpha = 1.6e-7, beta = 1.0 (ppmw per Pa) | 1.6 ppmw at 100 bar | from Yoshioka 2019 | Lichtenberg 2021 Table 1 | READ |
| CH4 | X = alpha p, alpha = 9.937e-8 ppmw/Pa (beta = 1) | 0.99 ppmw at 100 bar | Fe-free basalt, Ardia 2013 | Lichtenberg 2021 Table 1 | READ |
| CH4 | ppm = exp(4.93 - 1.93 (P_GPa - 1e-4)) f_CH4 (f in GPa) | as coded; gives 20 ppm at f = P = 1e4 bar; units not verified against the paper | haplobasalt, 0.7 to 3 GPa, 1400 to 1450 C | Ardia et al. 2013 GCA 114:52 Eqs. 7a, 8, via atmodeller | 2ND-HAND (flag: check units) |
| H2 | X = alpha p, alpha = 2.572e-6 ppmw/Pa | 26 ppmw at 100 bar | Hirschmann 2012 | Lichtenberg 2021 Table 1 | READ |
| N2 | X = alpha p^(1/beta); alpha = 7.0e-5, beta = 1.8 (oxidised); reduced set alpha = 74.16, beta = 4.582 | | IW-8.3 to IW+8.7 in Libourel 2003 | Lichtenberg 2021 Table 1; fO2 form of Dasgupta 2022 in atmodeller | READ |
| H2O | X = alpha p^(1/beta), alpha = 1.033, beta = 1.747 (ppmw, p in Pa); 1.05e4 ppmw at 100 bar | | compilation | Lichtenberg 2021 Table 1 | READ |

The three-regime Henry constants for CO, CO2 and CH4 under "oxidised / reduced / very reduced" conditions that Maurice et al. 2024 take from Hirschmann 2016 (Am. Min. 101:540, Table 2) are the likely origin of Boss's three-regime scheme. Maurice Table A.1 gives the calibration domains: CO and CH4 (Hirschmann 2016) 1400-1800 C, 350-32000 bar, IW-3.6 to IW+2.3 (basalt); CO2 (Hirschmann 2016) 1200-1625 C, 100-25000 bar, IW-1 to FMQ+2; CO2 (Eguchi and Dasgupta 2018) 950-1600 C, 500-30000 bar, FMQ-6 to FMQ+1.5 (foidite to rhyolite); CO (Yoshioka 2019) 1200-1600 C, 2000-30000 bar. The constants themselves were not reached (Hirschmann 2016 Table 2, Eguchi and Dasgupta 2018).

Physical versus schematic: physical (mass action or Henry): carbonate proportional to fCO2 with exp(-DeltaV (P-1)/RT); CO dissolved as a neutral molecule with exponent 0.5 to 1.0 on fCO; CH4 proportional to fCH4; H2 proportional to fH2; H2O proportional to f_H2O^0.5; the gas-gas ratios F3. Schematic or empirical: the NBO/T factor (really exponential in NBO/O with extra cation terms), the exp(-dIW) form, the fH2^2 X_H2O form, and any extrapolation of the laws above to 2000 K and above (none is calibrated there except the H2O law at 2173 K).

### F6 Carbon inventory check from the readable papers

Sossi 2020 (BSE, H = 0.01 wt%, C = 0.01 wt%, N = 0.0002 wt% in his Fig. 3 model, Earth calibration only): with existing laws about 99 percent of H2O and about 35 percent of CO2 are dissolved at 2173 K, N2 almost all in the gas. Gaillard 2022: carbon dissolved in the equilibrium magma ocean is 2 to 70 ppm C depending on fO2 and ocean size (depleted-mantle value 37 +/- 12 ppm C); graphite saturation below about IW-3; a total of 40 to 220 bar atmosphere. Earth's BSE C/H = 1.3 +/- 0.3, C/N = 40 +/- 8, C/S = 0.49 +/- 0.14 (Hirschmann 2016 abstract, READ): calibration targets only.

### F7 Lava-world vapour (dry melt, no H, C, N, S)

Use reaction-equilibrium vapour pressures, not an on-the-fly Gibbs minimiser: log10 P_i = [sum nu_ij mu_j + nu_O2 mu0_O2 - mu0_i]/(RT ln 10) + nu_O2 log10 fO2, with melt oxide activities from a solution model (VapoRock Eq. 3). Ideal-mixing oxide activities (older MAGMA) differ by factors 3 to 10 and up to 300 for Si and K, so a generator that wants correct Na/K/SiO ratios needs tabulated outputs of a MELTS-based run, not ideal activities. Both published codes are copyleft (VapoRock AGPL-3.0, LavAtmos GPL-3.0): do not vendor the code; re-implement Eq. 3 with JANAF data, or run the codes offline and store the output tables (licensing of outputs is a decision for Boss).

fO2 during congruent vaporisation of BSE (Schaefer and Fegley 2009, used by VapoRock): dIW +4.5 at 1500 K falling to +1 at 3000 K (VapoRock Section 4.3, READ). VapoRock Table 3 (1900 K, relative to Na = 1, dIW +2.6 to +3.4): K/Na 0.126 (tholeiite), 0.247 (alkali basalt), 0.113 (komatiite), 1.53 (dunite), 0.435 (Bishop Tuff); Fe/Na 0.0138, 0.0076, 0.0173, 0.0087, 0.0052; Si/Na 0.0044, 0.0014, 0.0052, 0.0003, 0.0126; Mg/Na 1.8e-4 (tholeiite), 7.6e-4 (komatiite). Total pressure over BSE melt: about 1e-7 bar at 1500 K, about 1e-3 bar at 2200 K, about 0.1 bar at 3000 K (Schaefer and Fegley 2009 Fig. 1; Ito Fig. 1).

Temperature classes (Miguel 2011, komatiite, MAGMA): Type I below 2033 K (Na, O2, O, Fe); Type II 2033 to 2588 K (SiO above Fe); Type III 2588 to 2890 K (Mg above Fe); Type VI 2890 to 3168 K (SiO above O2); Type V above 3168 K (silicate oxide dominant). With MELTS-based codes the SiO-over-Na crossover moves down to 2700 to 2800 K for BSE.

Mean molar mass of the vapour: 25 g/mol near 2000 K to 33 g/mol near 3000 K for BSE (F7 computation); use 34 g/mol for fractionated (Na-free) vapour (Kite 2016). These feed scale height and escape in the atmosphere lane.

### F8 Stored quantities and ΔIW parametrisation

Stored per rocky planet (new): mantle dIW at the magma-ocean surface (H21 scale, reference T stated), mantle Fe3+/Fe-total (derived through F2), pressure at the base of the magma ocean or melt fraction (Maurice's Phi), H, C, N, S inventories as mass fractions of the silicate (Earth calibration only: Sossi H 0.01, C 0.01, N 0.0002 wt%), surface temperature; from part B: mantle oxides (SiO2, MgO, FeO, Fe2O3 if kept, Al2O3, CaO, Na2O, K2O, TiO2, Cr2O3) and the melt NBO/O.

Is there a published distribution of dIW per planet? No. What exists:

| piece | figure | source | label |
|---|---|---|---|
| Solar System rocky bodies at core formation | span of about 6 dex in dIW (Mercury and enstatite chondrites lowest); Earth -1 to -2 depending on the FeO activity-coefficient ratio (1 to 4) | Doyle et al. 2019 (arXiv 1910.12989) | READ |
| Polluted white dwarfs (rocks accreted, dIW as at core formation) | -1.17, -1.37, -1.14, -0.59, < -1.21, -0.63 (six stars; the IAU proceedings version cites 15 and 16 stars in later papers, SNIPPET); most oxidised -0.6 would form no Fe core. Detection limit about -8 to -0.9 | Doyle 2019 Table S1 | READ |
| Magma-ocean surface dIW from metal-silicate equilibrium at the base (dIW -2, -1.5, -2) | Earth +2 (base 55 GPa), Mars -0.3 (14 GPa), Moon -2 (5 GPa); Earth falls by about 1.5 log units if the base is at 25 GPa instead of 55 GPa; 1-sigma 0.5 log unit | Deng et al. 2020 | READ |
| Equilibration dIW range proposed for Earth accretion | about -4 (Rubie 2011) to about -1 | quoted in Dalou 2019 introduction | 2ND-HAND |
| Outgassing of meteoritic mixtures | CI-like about QFM -1, EH-like about IW -1 | Schaefer and Fegley 2017 | READ |
| Crystal-mineralogy contribution to mantle fO2 spread across planets of known host star composition (Hypatia, N about 1000 to 1250 stars) at fixed Fe3+/Fe-total | sigma(log fO2) about 0.3 dex (1 GPa), 0.4 dex (4 GPa); dFMQ changes by about 4 log units when Fe3+/Fe-total goes from 1 to 9 percent; correlation of dFMQ with Mg/Si slope 1.4 to 4.0 per unit Mg/Si; "Fe3+/Fe-total is unconstrained for exoplanets" | Guimond et al. 2023 (arXiv 2308.09505) | READ |
| Putirka and Rarick 2019 (arXiv 1907.05506) | does not estimate fO2; uses Fe as a "combined proxy for fO2 and core-formation efficiency" | | READ |
| Hirschmann 2022 | model for Fe3+/Fe-total vs T, P, composition, fO2; no distribution | | READ |

Parametrisation proposal (a parametrisation, not a published distribution; the shape and bounds are my choice and must be labelled so in the design): draw the dIW at the base of the magma ocean from a uniform distribution on about -4 to -1 (the Rubie to Siebert range, which brackets Deng's -2, -1.5, -2 and Maurice's tested -5, -2, 0), then add a pressure-dependent upward shift from Deng's profile (about +4 at 55 GPa, about +2.5 at 25 GPa, about +1.2 at 14 GPa, 0 at 5 GPa, read from the four points above; a fit needs Deng's Fig. 3 data) to get the surface dIW; then add a Gaussian scatter of sigma 0.3 to 0.4 dex for composition (Guimond). Allow an oxidised tail to about IW +3.5 (Earth's upper mantle today, Sossi and Deng) for crystallised mantles after hydrogen loss (Gaillard 2022 on Venus). Solar System values are calibration targets only (Earth +0.5 to +2.3 for the magma ocean on different scales, Mars -0.3 to +0.5, Moon -2).

Consistency tests the TODO tree can assert: (1) for any stored dIW and T, f_CO2/f_CO and f_H2O/f_H2 equal the F3 table values within 1 percent; (2) element mass balance for H, C, N, S between melt, atmosphere and graphite closes to 1e-9; (3) F2 round trip: Fe3+/Fe-total to dIW and back; (4) QFM minus IW is 3.1 to 3.7 at 1 bar for 1500 to 2200 K.

## Still unreachable (Boss must research)

Ranked by effect on the generator.

1. Hirschmann 2016 Am. Min. 101:540, Table 2: Henry constants for CO, CO2, CH4 in "oxidised, reduced, very reduced" melt. Only the first page (abstract) opened (`msaweb.org/MSA/AmMin/TOC/2016/Abstracts/AM101P0540.pdf`). Would give the constants behind Boss's three-regime scheme.
2. Eguchi and Dasgupta 2018 Chem. Geol. 487:23: CO2 solubility from fluid to graphite/diamond saturation, 950-1600 C, 0.05-3 GPa, FMQ-6 to FMQ+1.5. Constants and T-dependence needed for magma-ocean T; ScienceDirect refused.
3. Armstrong et al. 2015 GCA 171:283 and Ardia et al. 2013 GCA 114:52: confirm Eq. 10 (Armstrong) and Eqs. 7a, 8 with units (Ardia, my evaluation of the coded form gives only 20 ppm CH4 at 1e4 bar, which looks low). ScienceDirect refused.
4. Duncan, Dasgupta and Tsuno 2017 EPSL 466:115: CO2 at graphite saturation versus NBO/T (0.4 to 0.9 wt% at IW+1.6, 1 GPa, 1600 C). Would test the NBO/T claim in G3.3a.
6. Deng et al. 2020 Supplementary Notes and Fig. 3 data: surface dIW versus magma-ocean base pressure, to fit the shift in F8.
7. Hirschmann 2012 EPSL 341-344:48 (magma ocean and early atmosphere), Schaefer and Fegley 2010 Icarus 208:438, Schaefer and Fegley 2004 Icarus 169:216 (MAGMA code and Io vaporisation), Dixon 1997 original (Am. Min. 82:368; msaweb path 404). None has an arXiv version that was found.
9. VapoRock journal reference (ApJ volume and page) and the VapoRock default fO2 grid. Fig. 4 partial pressures of BSE vapour are only in graphics: a data table (or a run of VapoRock) would replace my two-point mean molar mass.
10. G3.2b: the study that gives dIW -1.4 to -4.7 for carbon, if it exists.

## Counts

Rows: 26. Confirmed 5, partly 16, contradicted 3, unverifiable 2.

- CONFIRMED 5: G3.1a, G3.4b, G3.4d, G3.5b, G3.6e.
- PARTLY 16: G3.1c, G3.1d, G3.2a, G3.2c, G3.2d, G3.3a, G3.3b, G3.3c, G3.4a, G3.4c, G3.4e, G3.5a, G3.5c, G3.6a, G3.6b, G3.7.
- CONTRADICTED 3: G3.1b, G3.6c, G3.6d.
- UNVERIFIABLE 2: G3.2b, G3.4f.

Verification labels in this file: READ about 60 items, 2ND-HAND (open-source code constants) 14, COMPUTED 12, SNIPPET 6, RECALLED 1 (VO/FeH/CrH as hot-Jupiter opacities), NOT REACHED 8 sources (list above).

Files in `/mnt/project-files/design-drafts/parts/G-files/`: `G3_buffers_and_gas_ratios.py` and `G3_out_buffers_and_gas_ratios.txt`; `G3_ch4_equilibrium.py` and output; `G3_mean_molar_mass.py` and output; `G3_solubility_eval.py` and output; `G3_fe3_to_diw.py` and output; `G3_nasa_glenn_subset.csv`.
