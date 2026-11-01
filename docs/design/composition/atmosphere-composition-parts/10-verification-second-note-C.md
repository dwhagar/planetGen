# 10 Verification of the second owner note, part C (dust and haze, comet table, radioisotopes, stellar metallicity, giant-planet metallicity, licences)

Every claim in the note was treated as unverified data. Labels follow _worker-rules.md: READ (original text), 2ND-HAND (read in a reachable paper or code that cites it), DATA (data file in a package), COMPUTED (re-derived by me with Python from stated inputs), SNIPPET (search summary only, not relied on), RECALLED (memory, unverified, listed in section 9), NOT REACHED. Verdicts: CONFIRMED, PARTLY, CONTRADICTED, UNSUPPORTED (searched, nothing found either way), NOT REACHED. No repo file was edited. Results already in 08-verification-of-boss-added-research-B.md (Mars dust omega0 and k, Titan DISR laws, comet ratio ranges, licence identities of Mainz, Heays, Bremen) were not redone; they are cited where they matter.

Network note: the widened proxy let me reach arXiv PDFs, Crossref, Europe PMC, NTRS (API), Zenodo, IOPscience, nature.com (via curl; WebFetch cannot resolve it), PyPI, the Leiden and NMSU sites. Refused or dead: jpldataeval.jpl.nasa.gov, phidrates.space.swri.edu (no connection), agupubs, sciencedirect, repository.si.edu (403), www.geol.umd.edu PDF (404).

---

## 0. Two corrections to earlier files found on the way

1. [R0768] is Cheng et al., "Haze in Pluto's Atmosphere" (Icarus 290, 112; READ header and text). Part 08B section 2.1 attributes it to "Gao and Stam". The "too small backscatter lobe by factor about 3" statement is Cheng et al. (their aggregate model). [R1510], "Constraints on the Microphysics of Pluto's Photochemical Haze from New Horizons Observations", is [R1510] (READ: charge-to-radius ratio of about 30 e-/um, aggregates grow to 0.1 to 0.2 um).
2. [R3882], the id given for [R4352] in the note, is Sallum et al., "Accreting protoplanets in the LkCa 15 transition disc" (READ header). The real paper is [R4352].

---

## 1. Mars dust (claim 1)

Primary text read in full: Lemmon, Wolff, [R2644]. Abstracts read (Crossref): [R0809, R2645, R4742, R4743] (JGR 114, E00D04), [R0327, R3172]. [R2800] full text (HAL hal-01120102) is 2ND-HAND for the dust papers

| id | Sub-claim | Evidence | Label | Verdict |
|---|---|---|---|---|
| 1a | Modified gamma size distribution | Dust: gamma distribution, r_eff 1.5 um, effective variance 0.3, T-matrix cylinders of axial ratio 1 ([R2800] s2.2). [R2644] quotes variances 0.2 to 0.5 for Viking, Pathfinder, MER | 2ND-HAND; READ | CONFIRMED |
| 1b | r_eff 1.4 to 1.6 um quiescent | [R2645]: 1.47 +/- 0.21 (Gusev), 1.52 +/- 0.18 um (Meridiani). [R4742]: canonical 1.6 to 1.7 um. [R0809]: 1.5 +/- 0.1 um. [R2644]: 1.4 um early in the mission | READ | CONFIRMED |
| 1c | Above 2.1 um in regional and global storms | [R2644] abstract and text: "1 um effective radius during northern summer and a 2 um effective radius at the onset of a dust lifting event", "over 2 um with the onset of a dust storm" with the caveat that ice hazes depress the ratio. [R0809]: 1.8 to 2.5 um in the southern hemisphere at the peak of the 2001 global storm. [R4742]: departures from 1.6 to 1.7 um during 2001A. The value 2.1 um is one of three reference curves (0.7, 1.4, 2.1 um) drawn in Lemmon Fig. 10C, not a measured threshold | READ | PARTLY (about 2 to 2.5 um at storm onset; no sharp 2.1 cut) |
| 1d | IR(9 um) to 880 nm opacity ratio about 0.3 (northern summer, r_eff 0.7 um) to 0.5 (major events) | [R2645]: visible/9 um = 2.0 +/- 0.2, i.e. IR/visible 0.50 +/- 0.05 (r about 1.5 um). [R0809]: visible/IR 1.7 at r_eff 1.5 +/- 0.1 um (IR/vis 0.59) and 2.7 at northern Ls 50 to 200 deg with r_eff 1.0 +/- 0.2 um (IR/vis 0.37); COMPUTED reciprocals. [R2800]: Qext(0.67 um)/Qext(9 um) = 1.84 at r_eff 1.5 um, nu_eff 0.3 (IR/vis 0.54). [R2644] text: the northern-summer size is about 1 um, not 0.7 um; the 0.7 um is a plotted reference line; the plotted ratio is Mini-TES 9 um over MER 880 nm and "increases significantly with the onset of dust storms". [R4742]: TES PDS 9 um depth is an absorption depth about 30 percent below extinction | READ; 2ND-HAND | PARTLY. Supported: 0.37 at r 1.0 um, 0.5 to 0.6 at 1.4 to 1.5 um, rising with storms. "0.3 at 0.7 um" has no source; the curve it came from was not reached |
| 1e | omega0 "0.60 to 0.70 at 300 nm to 0.92 to 0.97 at 800 to 2000 nm" | [R4743] covers CRISM 440 to 2920 nm with MER; Madeleine says the MER-MRO analyses span 0.26 to 3 um. Channel means (Madeleine Table 2): 0.1 to 0.5 um 0.665 (Ockert-Bell) and 0.738 (Wolff); 0.5 to 5 um 0.928 and 0.955. Near-IR 1.0 to 2.6 um 0.97 (Vincendon, earlier result). [R0809] TES solar band 0.92 to 0.94. [R4743] spectral table itself NOT REACHED, so no value exactly at 300 nm | 2ND-HAND | PARTLY (800 to 2000 nm value fits 0.955 to 0.97; the 300 nm value is plausible inside a band mean of 0.67 to 0.74 but unconfirmed) |
| 1f | UV-blue absorption by octahedral Fe3+ charge transfer in ferric oxides (hematite, nanophase) in silicate matrices | No readable text states this for airborne dust as worded. [R0327] (READ abstract): Pathfinder soil and dust "consistent with ... poorly crystalline or nanophase ferric oxide(s)". [R3172] (READ abstract): "structured absorption edge due to crystal field transitions ... from weak absorption in the near-IR to intense absorption in the near-UV is characteristic of the ferric oxides" (crystal-field wording). NTRS 19940020405 (READ abstract): Viking-era analysis needed an extra 1 percent visible-UV absorber. [R3560]: no readable text reached | READ (abstracts) | PARTLY. Fe3+ in nanophase ferric oxide or palagonite is the standard attribution; "octahedral", "charge transfer" and "silicate matrices" are plausible (Sherman, Clancy) but not stated for the airborne dust in any text read |
| 1g | Early-morning brightness spikes from thin water-ice shells alter g(lambda) | [R2644] full text searched for ice, morning, coat, shell, asymmetry, g. Only: early-morning imaging "was especially rare"; the Mini-TES/880 nm ratio rise at Opportunity near Ls 30 to 130 deg "may indicate a sensitivity to ice-coated dust"; ice hazes depress the ratio; cirriform clouds at Ls 20 to 136 deg; possible water-ice haze. No spikes, no shell model, nothing on g | READ | UNSUPPORTED as worded (only the ice-coated-dust speculation exists) |

Generator-relevant facts read in [R2644]: 1.5 to 1.65 um particles with variance 0.2 to 0.5 have an extinction maximum near or beyond 1000 nm and 7 to 11 percent less opacity at 440 nm than at 880 nm; measured MER ratios correspond to an Angstrom exponent near zero (READ). Dust is grey in extinction in the visible; the blue/UV drop is absorption (omega0), not extinction slope.

---

## 2. Haze (claim 2)

Read: [R0156], "The Pale Orange Dot", Astrobiology 16, 873 ([R0156], full text); abstracts of [R4415] (PNAS 103, 18035, Crossref), [R1097] (Astrobiology 9, 447, Europe PMC), [R4883] (Nat. Geosci. 5, 359, nature.com page); [R1510] Pluto [R0768, R1510] searched for wavelength laws.

| id | Sub-claim | Evidence | Label | Verdict |
|---|---|---|---|---|
| 2a | Organic haze needs CH4/CO2 above about 0.1 | [R0156]: haze "would have formed when the ratio of CH4 to CO2 exceeded about 0.1 [R4415]". Their model (surface ratio): CH4/CO2 = 0.1 gives "insignificant haze ... spectrally indistinguishable from an atmosphere with no haze"; 0.2 gives the "thick" haze used for all hazy cases. Haze formation observed in the lab "down to C/O ratios as low as 0.1 [R1097, R4415]"; the O from CO2 photolysis is incorporated into carbonyl and carboxyl groups, so C/O below 1 makes oxygen-bearing haze (succinic acid by mass spectrometry), unlike Titan with C/O >> 1. [R1097] abstract: CO2/CH4 = 10 (CH4/CO2 = 0.1) with H2/CO2 up to 15; H2 gives an order-of-magnitude drop in haze production with no change in composition. [R4415] abstract: haze from a CH4/N2 and early-Earth mixture, early-Earth production about 1e14 g/yr (no threshold number in the abstract). [R4883] abstract: bistable organic haze/haze-free states driven by biological CH4; the numeric threshold is not in the abstract and the paper is paywalled | READ (Arney text, abstracts); 2ND-HAND for Trainer's 0.1 | CONFIRMED as an onset (0.1), with the corrections that a thick haze needs about 0.2 and that added H2 cuts the rate about 10 times. "C/O" in the lab statement is not the C/O of the CH4+CO2 gas (that is (1+r)/2 for r = CH4/CO2, COMPUTED, never below 0.5); the quoted 0.1 appears to be CH4 carbon over CO2 oxygen (r/2 = 0.1 at r = 0.2), my reading, flagged |
| 2b | tau_haze proportional to lambda^-alpha, alpha 2.0 to 4.0 in the visible, flattening in the IR "as fractal aggregate monomers coalesce" | Titan DISR (READ in 08B, Tomasko team slides): cumulative tau = 1.012e7 lambda^-2.339 (above 80 km), 2.029e4 lambda^-1.409 (30 to 80 km), 6.27e2 lambda^-0.9706 (below 30 km), lambda in nm. COMPUTED whole-column exponents from the 08B totals: 1.72 (500 to 1500 nm) and 1.77 (355 to 934 nm). So Titan never reaches 3 or 4 and only the top layer (2.34) falls in the note's lower bound. [R0156] (READ): fractal particles (Df 1.5 to 2.4, monomers 0.06 um) give "more extinction in the ultraviolet but less in the visible and near infrared" than equal-mass spheres and are more forward scattering; no single exponent is quoted. Pluto: haze is blue (peak I/F 0.7 to 0.8 in the blue filter), consistent with Rayleigh scattering by r about 0.01 um monomers but with a forward-scattering phase function that needs r > 0.1 um, fitted by aggregates of 0.01 um spheres, Df 2, bulk radius 0.15 um ([R0768], READ; [R1510] gives no exponent. A Rayleigh sphere gives 4 only when scattering dominates and k is flat; for an aggregate above the Rayleigh-Gans knee the exponent drops by about Df (my derivation, not a source, flagged). Aggregates keep their monomers (Titan: about 3000 to 4000 monomers of 0.04 to 0.05 um, 08B); coating that rounds particles is discussed only for the lowest layers [R0768] | READ; COMPUTED | PARTLY. "2 to 4" is correct only for the smallest particles (top of Titan's haze 2.34, Rayleigh limit 4); aggregate columns give 1 to 2.3 (Titan 0.97 to 2.34, whole column about 1.7). "Monomers coalesce" is not a source-supported mechanism; flattening comes from aggregate size exceeding the wavelength |

[R2009] ([R2009], Nature Astronomy 2, 303; 08B READ): production rates spread over about 250 times for sub-Neptune analogues at 300 to 600 K and 100 to 10000 times solar metallicity. [R1513] (Nature Astronomy 4, 951, [R1513]): silicate/hydrocarbon crossover near Teq 950 K (as used in 09 row n6; not re-read here).

---

## 3. Comet table (claim 3)

### 3.1 Photochemical lifetimes at 1 AU

Sources reached: sbpy 0.6.0 (PyPI sdist, file sbpy/activity/gas/data/__init__.py, DATA, quoting [R2040] = HM15, [R0934] = C94) Leiden photo-rate table ([R1873] data; page photo_rates.html, "Solar" column times the stated factor 37700 for 1 AU, 91.2 to 200 nm only, quiet-Sun spectral shape, READ; lifetimes COMPUTED by me); [R1163] ([R1163], NH3 rate, READ). HM15 itself and C94 full tables NOT REACHED (Elsevier and Wiley refused; NTRS record of HM15 is metadata only).

| Species | Note (s) | Best reached value | Source and label | Verdict |
|---|---|---|---|---|
| H2O | 8.0e4 | HM15 quiet 8.294e4, active 4.939e4 (ratio 0.595); [R0837] intermediate 5.2e4; quiet 8.5e4; Leiden solar 1.15e5 | DATA sbpy; READ 2307.11486; COMPUTED | CONFIRMED (quiet Sun; active Sun is 0.6 times) |
| CO | 1.5e6 | CE83 1.5e6 (sbpy); Leiden solar 2.6e6 (dissociation); HM15 rate range 0.075 to 0.188e-5 /s (SNIPPET, implies about 5e5 to 1.3e6 s) | DATA; COMPUTED; SNIPPET | CONFIRMED as the CE83 value; HM15 may be shorter |
| CO2 | 5.0e5 | CE83 5.0e5; HM15 quiet 4.948e5, active 2.101e5 (ratio 0.425); Leiden solar 1.22e6 | DATA; COMPUTED | CONFIRMED (quiet) |
| HCN | 7.0e4 | C94 6.7e4; Leiden solar 1.26e5; HM15 rate 1.31 to 3.24e-5 /s (SNIPPET) | DATA; COMPUTED | CONFIRMED |
| CH4 | 1.2e5 | Leiden solar (dissociation 4.83e-6 /s; ionisation adds 3.5e-7) 2.07e5 (A+ class, 20 percent). No tabulated Huebner value reached | COMPUTED | UNSUPPORTED (factor 1.7 from the only value reached; Huebner/Crovisier value NOT REACHED) |
| C2H6 | 1.0e5 | Leiden solar 1.67e5 (dissociation 5.99e-6 /s, A+). Huebner value NOT REACHED | COMPUTED | UNSUPPORTED (factor 1.7) |
| NH3 | 2.5e4 | k = 1.8e-4 /s, "average rate for a quiet and an active Sun from [R2041] and HM15", giving 5.6e3 s; [R1873] gives 1.5e-4 /s (6.7e3 s); Leiden solar COMPUTED 1.538e-4 /s, 6.5e3 s | **Superseded by part 10E:** recommended 1.5e-4 /s (6.7e3 s), range 1.3e-4 to 1.8e-4 /s; the Leiden Solar column is NOT limited to 91.2-200 nm for NH3 (it uses a 67-420 nm solar file and the cross section to 226 nm), so the statement above that it covers only 91.2 to 200 nm is wrong for this species | READ 2307.11486; COMPUTED | CONTRADICTED (note about 4 times too long) |

Leiden-derived values run 1.3 to 2.5 times longer than HM15 for H2O, CO2 and HCN because the table covers 91.2 to 200 nm only; use them as order-of-magnitude fallbacks (spread about a factor 2). Other values read: C94 CH3OH 7.7e4 s, H2CO 5.0e3 s; HM15 CN 3.15e5 (quiet) and 1.35e5 s (active), from H92 in sbpy.

### 3.2 Photolysis products and diagnostic bands

Products at the threshold wavelength (Leiden cross-section page, READ): H2O to H + OH; CO to C + O; CO2 to CO + O; CH4 to CH3 + H; C2H6 to C2H5 + H; HCN to H + CN; NH3 listed as NH + H2 at threshold (301 nm) but NH2 + H is the dominant channel in the solar field ([R1163] citing [R1873], READ). The site warns that branching at shorter wavelengths (for example CH4 to CH2 + H2) is often dominant; wavelength-resolved branching is on its H2O and NH3 pages only.

| Band or line | Note | Check | Label | Verdict |
|---|---|---|---|---|
| OH 18 cm | yes | 1.667 GHz is 18.0 cm | COMPUTED | CONFIRMED |
| OH 309 nm | yes | "OH (0-0) band around 308 nm" (Biver et al., Comets III, [R0410]) | READ | CONFIRMED (308 to 309 nm) |
| CO 4.67 um | yes | No read text gives the band wavelength; 2143 cm-1 (RECALLED) converts to 4.666 um | COMPUTED from RECALLED | NOT REACHED (consistent with memory) |
| CO 230 GHz | yes | CO 2-1 at 230.538 GHz ([R1088] Table) | READ | CONFIRMED |
| CO2 4.26 um | yes | "very strong band near 4.26 um" (Comets III) | READ | CONFIRMED |
| CH4 3.30 um nu3 | yes | "strong IR fundamental band (nu3) near 3.3 um"; needs |d(delta)/dt| above about 15 km/s to avoid telluric lines | READ | CONFIRMED |
| C2H6 3.35 um nu7 | yes | "nu7 band near 3.35 um" | READ | CONFIRMED |
| NH3 "1.15 mm" | note queries it | The 572 GHz line is 0.524 mm (COMPUTED); 572 GHz is "only accessible from space" ([R1163] citing [R0407]). 1.15 mm is wrong (1.15 looks like CO 1-0 at 115 GHz = 2.6 mm). Ground radio uses inversion lines JK = 11-11 to 55-55 near 23.7 GHz (12.7 mm; frequency RECALLED). Near-IR NH3 at 2.9 to 3.0 um (RECALLED; Comets III says the 3.0 um region is dominated by H2O) | READ; COMPUTED | CONTRADICTED (1.15 mm) |
| HCN 88.6 GHz | yes | 88.632 GHz J = 1-0 (not found in a read text) | RECALLED | NOT REACHED |
| CN violet 388 nm | yes | 3883 A (not found in a read text) | RECALLED | NOT REACHED |

### 3.3 Abundance ranges (percent of water)

Compared with 08B (Comets III Table 3 READ, plus [R0859, R4705]):

| Species | Note | Found | Verdict |
|---|---|---|---|
| CO | 0.5 to 30 | IR 0.3 to 26, radio < 1.23 to 35 (limits), median 3.3, Hale-Bopp 23 to 26 | PARTLY (bottom 0.3) |
| CO2 | 2 to 20 | 2 to 30 (IR 4 to 30), median 12 | PARTLY (top 30) |
| CH4 | 0.2 to 1.5 | 0.15 to 2.7 (OCC mean 0.88) | PARTLY |
| C2H6 | 0.1 to 1.5 | IR 0.1 to 2.7 (OCC mean 0.63, JFC 0.34) | PARTLY (top 2.7) |
| NH3 | 0.2 to 1.5 | radio 0.18 to 0.60, IR 0.1 to 3.6 | PARTLY (top 3.6; radio-only top 0.6) |
| HCN | 0.05 to 0.5 | radio 0.05 to 0.25, IR 0.03 to 0.5 | CONFIRMED (0.03 to 0.5) |

---

## 4. Radioisotope reservoir table (claim 4)

Sources read: [R4890], Planetary Science Journal, Feb 2022 ([R1003]; IOPscience full text incl. Table 1, READ; CI values from [R0270, R2729]; Moon, Mars, Vesta, Earth from the paper's own estimates); [R2032] ([R2032: Table 3], READ; adopts McDonough and Sun 1995 for the [R4687] for oceanic crust: 2ND-HAND); [R3483] (Science 333, 1850; abstract READ via experts.arizona.edu; K/Th lines 360 and 5500 READ in the JHU/APL highlight page); [R4266] (JGR 92, E537; abstract READ via Crossref); [R0138] (SNIPPET: 13,800 +/- 2600, MORB source 19,000 +/- 2600, OIB 11,900, crust 13,000). Ratios below are mass ratios (g/g), COMPUTED from the stated concentrations.

| Reservoir | K/U note (x1000) | Reached value | Th/U note | Reached value | Verdict (K/U; Th/U) |
|---|---|---|---|---|---|
| CI chondrites | 60 to 70 | 71.4 ([R0270]: K 550 ppm, U 7.7 ppb), 66.1 ([R2729]: 539 ppm, 8.16 ppb), 68.5 ([R0169] meteoritic eps, COMPUTED); CM 36.4, CV 18.2, H 65 | 3.6 to 3.8 | 3.68 (Barrat), 3.65 (Lodders), 3.88 (Asplund) | CONFIRMED (top 71); PARTLY (3.65 to 3.88) |
| Bulk silicate Earth | 10 to 12 | 9.7 +/- 2.1 ([R1003], from flat-pattern basalts); 13.8 [R0138], 14.05 [R3417], 12.1, 11.8 [R3363] quoted in Dauphas; 14.0 (McDonough and Sun 1995: K 280 ppm, U 20 ppb) | 3.7 to 3.9 | 3.95 (McDonough and Sun: Th 79 ppb), 3.68 (Dauphas, set by CI) | PARTLY (literature spans 9.7 to 14.05); CONFIRMED |
| Oceanic crust (MORB) | 10 to 12 | 10.0 ([R4687] via [R2032]: K 0.07 percent, U 0.07 ppm); flat MORB median 10.6 +/- 0.5 (Dauphas); MORB source 19.0 (Arevalo, SNIPPET) | 3.7 to 4.0 | 3.0 (Th 0.21, U 0.07 ppm); depleted mantle 2.75 (Th 0.022, U 0.008 ppm) | CONFIRMED; CONTRADICTED (about 3.0) |
| Continental crust (not in note) | | upper 8.6 (K 2.32 percent, U 2.7), bulk 11.5 to 11.8 (part 07) | | upper 3.9, bulk 4.3 to 4.5 (part 07) | |
| Moon | 2.5 to 3.0 | 2.45 +/- 0.14 (bulk silicate Moon, Dauphas: K 44 ppm, U 18 ppb) | 3.6 to 3.8 | 3.67, but this equals the CI ratio the paper assumes; no independent lunar measurement reached | PARTLY (best 2.4 to 2.6); UNSUPPORTED |
| Mars | 16 to 20 | 17.3 +/- 1.8 (bulk silicate Mars from SNC meteorites, Dauphas); K/Th 5300 mean crust times Th/U 3.67 gives 19.5 (COMPUTED, depends on assumed Th/U) | 3.7 to 3.9 | 3.68 assumed in Dauphas Table 1; no U measurement of Mars crust reached | CONFIRMED; UNSUPPORTED |
| Venus | 10 to 13 ("Earth-like, confirmed by Venera gamma-ray data") | Vega 1: K 0.45 +/- 0.22 percent, U 0.64 +/- 0.47 ppm, Th 1.5 +/- 1.2 ppm gives K/U 7.0 (+/-6.2), Th/U 2.3 (+/-2.6); Vega 2: K 0.40 +/- 0.20, U 0.68 +/- 0.39, Th 2.0 +/- 1.0 gives 5.9 (+/-4.5), 2.9 (+/-2.2) (READ abstract, COMPUTED). Venera 9 and 10 (SNIPPET values K 0.47 and 0.30 percent, U 0.60 and 0.46 ppm, Th 3.65 and 0.70 ppm): K/U 7.8 and 6.5, Th/U 6.1 and 1.5. Venera 8 (SNIPPET, K 4.0 percent may be K2O): 18.2 (+/-8), Th/U 2.95 | 3.7 to 3.9 | individual landers 1.5 to 6.1, each with 30 to 100 percent errors; unconstrained | CONTRADICTED (landers give about 6 to 8, with one outlier; Surkov states only "similar to basic rocks ... tholeiitic basalts and gabbros", no Earth K/U claim); UNSUPPORTED |
| Mercury (surface, northern hemisphere) | 8 to 18 | K 1150 +/- 220 ppm, Th 220 +/- 60 ppb, U 90 +/- 20 ppb gives K/U 12.8 +/- 3.7, K/Th 5230 +/- 1740 (cluster line K/Th = 5500 with Mars, Venus, Earth; Moon on 360) (READ, COMPUTED); Dauphas Table 1 K/U 14.3 +/- 4.1 | 3.5 to 4.0 | 2.44 +/- 0.86 (COMPUTED from the abstract) | CONFIRMED (surface, not bulk); PARTLY (central value 2.4, note's 3.5 is 1.2 sigma above) |

Also: [R3484] (JGR Planets, [R3484], Crossref abstract READ): surface K varies 300 to 2400 ppm, Th varies little, so K/Th varies widely across the surface; the 2011 average is therefore a northern-hemisphere mean, not a bulk value.

Summary of corrections for the generator: K/U by body class (mass ratio) CI 66 to 71k; BSE 9.7 to 14k; oceanic crust 10 to 10.6k; continental crust 8.6 to 11.8k; Moon 2.4 to 2.6k; Mars 17 to 19.5k; Mercury surface 12.8 to 14.3k (+/- 3.7 to 4.1k); Venus surface 6 to 8k (+/- 30 to 100 percent). Th/U: CI 3.65 to 3.88; BSE 3.7 to 3.95; upper crust 3.9 to 4.0, bulk crust 4.3 to 4.5; MORB 3.0, depleted mantle 2.75; Mercury surface 2.4 +/- 0.9; Venus unconstrained.

---

## 5. Stellar abundances (claim 5)

Read: [R0478, R0479] ([R0479], MNRAS; 53 thin-disc solar twins, HARPS, Th II 4019 A, Table 4 slopes); [R1446] ([R1446], ARNPS 68, 237); [R0309] ([R0309]; 30 elements, 79 solar twins; no Th or U); [R4469].

| id | Sub-claim | Evidence | Label | Verdict |
|---|---|---|---|---|
| 5a | Sub-solar metallicity has elevated [Th/Fe] | Botelho Table 4 (READ): [Th/Fe]_ZAMS = -0.338 +/- 0.049 [Fe/H] + 0.082 (rms 0.040), [Th/H]_ZAMS slope +0.671 +/- 0.049, but only over the narrow solar-twin range; mean [Th/Fe]_ZAMS +0.086 +/- 0.008 (rms 0.047), no correlation with age. Halo: median [Eu/Fe] about +0.3 at [Fe/H] about -1.5, r-I stars ([Eu/Fe] +0.3 to +1.0, about 130 known, about 15 percent of metal-poor stars) and about 30 known r-II stars ([Eu/Fe] > +1.0); Th follows Eu not Nd (Botelho; Nd about 60 percent s-process) | READ | PARTLY (small negative slope in thin-disc twins; large only in r-process-enhanced halo stars; Th tracks r-process timing, as claimed) |
| 5b | Same for [U/Fe] | Frebel: uranium "very challenging to detect", one weak blended optical line; "only five stars have reliable abundance measurements", all r-process-enhanced metal-poor stars (CS 31082-001, HE 1523-0901, others). [R4469] notes variation in U is not yet known for metal-rich planet hosts. Disc-star [U/Fe] vs [Fe/H]: no measurement exists in the texts read | READ | UNSUPPORTED |
| 5c | Hence more radiogenic heat per unit silicate mass | Botelho: [Th/Si]_ZAMS (the Th-per-mantle-mass proxy, citing [R4469]) falls with age: slope -0.0121 +/- 0.0018 dex/Gyr, from +0.149 dex at zero age for the youngest stars to +0.045 dex for 8.6 Gyr-old stars, so older, more metal-poor thin-disc stars have less Th per silicon. [R4469]: Th between 59 and 251 percent of solar in the solar-analogue hosts studied (READ in part 07). Heat also depends on decay since formation; K-40 and U-235 have decayed most in old systems | READ | CONTRADICTED in the thin-disc sample read; untested for thick-disc and halo hosts |

---

## 6. [R4352] (claim 6)

Primary text, [R4352] (ApJ 831, 64), READ in full. The note's id 1511.07456 is a different paper (section 0).

| Item | Note | Paper | Verdict |
|---|---|---|---|
| Z_planet/Z_star | proportional to M^-alpha, alpha 0.45 to 0.55 | Fit (9.7 +/- 1.28) M_J^(-0.45 +/- 0.09), planet mass in Jupiter masses (Fig. 11) | PARTLY (0.45 correct; the stated band 0.45 to 0.55 is not the paper's 1-sigma 0.36 to 0.54) |
| Heavy-element mass | not given | M_Z = (57.9 +/- 7.03) M_J^(0.61 +/- 0.08) Earth masses; "roughly M_Z proportional to sqrt(M)", Z_planet proportional to 1/sqrt(M) | CONFIRMED (the caller's guess) |
| Scatter | "drivers listed" | Intrinsic spread of M_Z about the fit: factor 10^sigma = 1.82 +/- 0.09 (1 sigma in log); "accurate to a factor of a few". Named drivers: observational uncertainty, migration history, stochastic formation. No relation to orbital properties or stellar mass. M_Z vs [Fe/H]: Kendall tau 0.088, p = 0.38. Z_p/Z_star vs mass: tau 0.44, p = 1.3e-5 | PARTLY |
| "Below 30 Earth masses accrete atmospheres enriched 10 to 50 times solar" | | Sample is 47 planets with 20 M_E < M < 20 M_J, so 20 to 30 M_E is the extreme low end; enrichment is relative to the host star's metallicity, and the paper does not discuss atmospheres below 20 M_E. Table 1 (M < 0.1 M_J, Z_p/Z_star): K2-24 c 11.3, K2-24 b 18.4, Kepler-30 d 21.4, HAT-P-11 b 27.8, CoRoT-24 c 29.2, K2-27 b 43.9, GJ 436 b 132. Fit values (COMPUTED): 34 at 20 M_E, 28 at 30 M_E, 16 at 100 M_E, 9.7 at 1 M_J; M_Z/M 0.54, 0.46, 0.29, 0.18 | PARTLY (about 11 to 130, fit 28 to 34; metals are in the envelope: "metal-enriched giant planet atmospheres should be the rule") |

---

## 7. Replacement list R1 to R15 versus part 09 section 4 (claim 7)

Verdict PARTLY (only the R15 numbering clash can be shown). Limit: I was given the note's names for six items, not its full list; the note's own numbers are known to me only for R15. The rest are matched by topic. Our list is part 09 section 4.

| Ours (part 09 s4) | Rows it feeds | Note item as named to me | Match |
|---|---|---|---|
| R1 envelope occurrence by radius, mass, insolation | b4, e3, e4, e5 | Radius valley and envelope occurrence | same topic |
| R2 volatile delivery vs formation distance | b1 to b6 | Volatile delivery vs formation distance | same topic |
| R5 giant-planet C/O, C-rich host fraction | f7 | Giant C/O | same topic |
| R9 comet coma classes and ratios | k1 to k5, k7 | Comet taxonomies | same topic |
| R11 trace gases | m2 to m6, m11 | Trace gases | same topic |
| R15 brown-dwarf type-Teff relation, T/Y boundary near 500 to 700 K, hot end above 2200 K | g1, g4 | R15 "Brown dwarf and exoplanet cloud regimes" | MISMATCH: that is our R12 (clouds, hazes: rows n6, n7), not R15 |
| R12 clouds and hazes: [R1513] crossover width, [R2009] yield, H2SO4-H2O curve, perchlorate/carbonate decomposition | n6, n7 | (no note item for H2SO4, perchlorate) | partial |
| R3 fO2 vs mass; R4 shoreline constants; R6 hot-Jupiter thermochemistry; R7 lava-world vapour; R8 frost atmospheres; R10 WD and pulsar planets; R13 radioisotopes; R14 initial rotation / UV | a1 to a4, d6, h2 to h4, i6, j1 to j5, l1 to l3, o2 to o6, p6 | not named in what I was given | unknown |

Which note items would actually source our DESIGN rows:
- "Brown dwarf and exoplanet cloud regimes" feeds n6 (silicate versus haze crossover, DESIGN R12) only through papers the note must name: [R1513, R1515] (already used), [R1512] (ApJ 855, 86, [R1512], abstract READ: f_sed from 0.5 to above 1 depending on nucleation/surface energy and K_zz, not gravity), [R3592] (ApJ, [R3592], abstract READ: bin-scheme size distributions for silicate and titanium clouds on hot Jupiters, cloud mass anticorrelated with Teq, bimodal sizes), [R0024] (ApJ 556, 872, f_sed formalism). These give particle sizes and vertical structure, not weights, and give nothing for g1 or g4 (type-Teff), which need an observational compilation (candidates RECALLED: [R1370, R2374]; verify).
- "Comet taxonomies" can source k1 to k5 only with a per-comet table (Harrington [R3535], Dello [R0048, R1063]). The class fractions in k1 to k3 remain SNIPPET until that table is read.
- "Giant C/O" feeds f7 with [R0535, R0536]; the [R4352] scatter in section 6 gives the Z_p/Z_star spread for the metallicity bins.
- "Trace gases", "Radius valley", "Volatile delivery" map one-to-one to R11, R1, R2 and add nothing without named papers ([R1478]; Owen and Wu 2017; [R4335] for the first two groups; [R1137, R1813] for trace gases are in our own list).

---

## 8. Licences (claim 8): what the pages actually say

| Item | Note | Page text found | Label | Verdict |
|---|---|---|---|---|
| JPL Publication 19-5 (Evaluation 19) | "public domain (US Government work)" | NTRS record 20210006316 (title "Chemical kinetics and photochemical data for use in atmospheric studies; evaluation number 19", JPL, 2020-05-01, 14 authors, no downloadable file): rights field determinationType OTHER, belongsToPublisher true, belongsToUsGov false, approver comment "Items in TRS are protected by copyright, but are furnished with U.S. government purpose use rights." The abstract says the evaluation is at jpldataeval.jpl.nasa.gov (refused, so its own terms are NOT REACHED). Evaluation 18 (NTRS 20160006735) is PUBLIC_USE_PERMITTED; Evaluation 15 (20090033862) PUBLIC_USE_PERMITTED | READ (API JSON) | CONTRADICTED for 19-5 as "public domain"; the NTRS entry says copyright-protected, government-purpose rights |
| [R2040] (not in the note's list but cited) | | NTRS has a CHORUS metadata-only accepted-manuscript record (id 23332861698053, funder NNX08AF03G, Elsevier BV, no file, no rights field); phidrates.space.swri.edu unreachable | READ | terms NOT REACHED |
| NMSU Karkoschka files | "proprietary, distribute only derived tables" | The data page saturn_infrared.html says only "Source: Karkoschka, E., 1994 ... Icarus, 111, 174-192" and offers a download; infraredsaturn.txt has only a column header; the directory README and documentation files are about page building; the site footer "Copyright" links to the JPL image policy. No terms text on any of these | READ | UNSUPPORTED (no licence stated; "proprietary" is not on the pages) |
| Leiden photodissociation database | "academic citation request / open access" | index page, exact text: "When using material from this site, please cite one of the overview papers. In cases of specific molecules, please cite the original source of the cross sections as appropriate." downloads.html and references.html: no licence words. Nothing says "open access" | READ | PARTLY (citation request only) |
| Bremen (Serdyuchenko O3) | "proprietary/unspecified" | 08B: no licence or terms on the data page (only a data privacy link) | READ (08B) | PARTLY ("unspecified" is right, "proprietary" is not shown) |
| Mainz atlas | CC BY 4.0 | Zenodo record 6951 metadata (API): license cc-by-4.0, access open, publication date 2013-06-24, creators Keller-Rudek, Moortgat, Sander, "frozen snapshot of all cross section and quantum yield data files"; points to the live atlas for the up-to-date version. The live site states only a citation request (08B) | READ | CONFIRMED for the 2013 snapshot; the live site has no licence |

Rule for planetGen stays as in 08B: vendor only what carries a licence line; otherwise read user-supplied copies.

---

## 9. What this changes in the generator items

GEN.414 (aerosols):
- Dust: draw r_eff from 1.0 um (northern summer) through 1.4 to 1.6 um (quiescent) to 2 to 2.5 um (storm onset); IR(9 um)/visible(0.88 um) opacity ratio 0.37 at 1.0 um, 0.5 to 0.6 at 1.5 um [R0809, R2645, R2800]; visible Angstrom exponent near 0; omega0 about 0.955 to 0.975 red/NIR, band mean 0.67 to 0.74 for 0.1 to 0.5 um. Do not use "0.3 at 0.7 um" or the 2.1 um cut or the ice-shell mechanism (no source). Attribute UV absorption to Fe3+ in nanophase ferric oxide or palagonite with a hedge.
- Haze: haze rows keep the onset CH4/CO2 = 0.1, treat 0.2 as thick [R0156], and cut production about 10 times when H2 is abundant [R1097]. Visible exponent: use Titan's 0.97 to 2.34 for aggregates and 3 to 4 only for tiny spheres; do not apply 2 to 4 to fractal hazes. Replace the mechanism "monomers coalesce" by "aggregate size exceeds the wavelength".
- Giant-planet envelope metallicity (row f): Z_p/Z_star = 9.7 M_J^-0.45 with spread factor 1.82 (valid 20 M_E to 20 M_J).

Comet rows (k): widen NH3 to 0.1 to 3.6 (IR) and note radio 0.18 to 0.60; CO2 to 30; CH4 and C2H6 to 2.7; coma lifetimes: use H2O 8.3e4 s (quiet) to 4.9e4 s (active), CO 1.5e6, CO2 4.9e5 to 2.1e5, HCN 6.7e4, NH3 5.6e3 to 6.7e3 (not 2.5e4), CH4 and C2H6 about 2e5 and 1.7e5 with a factor-2 band until Huebner's table is read.

Ground radioisotope rows (o): o3 K/U range LU[8.9e3, 1.4e4] is right for the Earth; add body-class values from the table in section 4 (Moon 2.4k, Mars 17 to 19.5k, Venus 6 to 8k, Mercury 12.8k +/- 3.7k, CI 66 to 71k). o4 Th/U U[3.8, 4.5] excludes MORB (3.0), depleted mantle (2.75) and Mercury (2.4); keep it for BSE and continental crust and add a lower tail for mantle-derived reservoirs. o1 and o2 remain Th-only; there is no disc-star U measurement; the "sub-solar metallicity means more heat per silicate mass" claim should not be added (5c).

Licences: do not mark JPL 19-5 as public domain; treat JPL 19-5, Huebner and Mukherjee, NMSU and Bremen as "no licence found, read user-supplied copies"; Mainz Zenodo snapshot is CC BY 4.0; Leiden is cite-on-use.

## 10. NOT REACHED

1. [R0934, R2040] full tables (quiet and active Sun, all species); CH4 and C2H6 lifetimes from either.
2. [R4743] omega0(lambda), (n, k) and Fig. 13 ratio-vs-size curves; [R4744] MARCI UV omega0; [R3560]; any text stating the airborne-dust Fe3+ mechanism.
3. [R4883] numeric haze threshold; [R4415] text; Pluto colour slope.
4. [R0138] text; measured Th/U for lunar and martian rocks; Venera 8 to 10 original papers (values SNIPPET); [R4687] (via [R2032] only).
5. Terms: jpldataeval.jpl.nasa.gov, phidrates.space.swri.edu, any Bremen terms.
6. Disc-star uranium abundances; thick-disc and halo [Th/Si] vs [Fe/H].
7. The note's full R1 to R15 list (only six names were provided).

## 11. Recalled / to verify

CO 2143 cm-1 (4.67 um), CO2 2349 cm-1, CH4 3019 cm-1, C2H6 2985 cm-1, NH3 3337 and 3444 cm-1 (wavenumbers); NH3 inversion lines near 23.7 GHz; HCN 88.632 GHz; CN 3883 A; [R1370, R2374] as type-Teff sources. Venera 8, 9, 10 values and Arevalo/Taylor values are SNIPPET.

## 12. Verdict counts (61 sub-claim rows; the per-row verdicts above are authoritative)

CONFIRMED 21: 1a, 1b, 2a, 3.1 (H2O, CO, CO2, HCN), 3.2 (OH 18 cm, OH 309 nm, CO 230 GHz, CO2 4.26 um, CH4 3.30 um, C2H6 3.35 um), 3.3 HCN, section 4 (CI K/U, BSE Th/U, MORB K/U, Mars K/U, Mercury K/U), 6 (M_Z fit), 8 (Mainz).
PARTLY 21: 1c, 1d, 1e, 1f, 2b, 3.3 (CO, CO2, CH4, C2H6, NH3), section 4 (CI Th/U, BSE K/U, Moon K/U, Mercury Th/U), 5a, 6 (alpha, scatter, low-mass enrichment), 7, 8 (Leiden, Bremen).
CONTRADICTED 7: 3.1 NH3 lifetime, 3.2 NH3 "1.15 mm", section 4 (MORB Th/U, Venus K/U), 5c, 6 (arXiv id), 8 (JPL 19-5 public domain).
UNSUPPORTED 8: 1g, 3.1 (CH4, C2H6), section 4 (Moon Th/U, Mars Th/U, Venus Th/U), 5b, 8 (NMSU).
NOT REACHED 4: CO 4.67 um, HCN 88.6 GHz, CN 388 nm (all RECALLED only), Huebner and Mukherjee terms.
