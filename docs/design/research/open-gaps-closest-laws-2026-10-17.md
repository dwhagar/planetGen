# Open gaps and the closest laws that cover each (2026-10-17)

Filed by the TODO thread on 2026-10-17 at about 1:16 PM EDT (as this thread's clock reads); the notes on the items are in docs/TODO.md and the summary is section 4.11 of [gap-filling-plan.md](gap-filling-plan.md).

Status: research result, read-only. Written 2026-10-17, about 1:10 PM EDT as this thread's clock reads (the lane clock can drift). No code, no TODO.md, issue or docs edits: the TODO thread files it. Answers Boss's question of 12:37 PM EDT: a list of every outstanding item and, for each, the nearest law, model or standard calculation, where to look for it, and how good the stand-in is.

Inputs: the gap-filling plan and its registry, the boss-papers report ("still open" list and the contradictions table), the Gemini batch 1 check, GEN.531 (R1 to R15), the papers-to-fetch list, the tunable-gaps and stellar-evolution-tables reports, and the TODO items GEN.450, 453, 463, 467, 469, 493, 494, 502 and 529. Filled items are not listed as open; they are in section 4.

## 0. How to read the tables

**Columns.** Each open item has: what is missing; the nearest law or model with its inputs, in words (the equation form is given when it is short); where to look (a book or review by topic, then a free source or tool Boss can use); the stand-in the project uses today with its quality; and the rule that decides it.

**Rule marks** (Boss, 12:39 PM EDT on 2026-10-17): **V** = when in doubt, verified wins (a verified source decides, or will replace the stand-in when one arrives; no conflict between verified sources today). **N** = when two verified sources are in doubt, the newest wins. **M** = when two sources are both in question, "most permission take wins", which Boss confirmed by card at 1:12 PM EDT means the widest range of the two sources. **D** = Decided earlier: the six answers of "widen, keep both" (GX2 sign, GX3 slope, GR7 ratio, thick-disc height ratio, thick-to-thin mass ratio, cluster-system radius ratio) stand; Boss confirmed at 12:40 PM EDT that newest-wins applies only to later conflicts. Where a row names a conflict that is not one of the six, it says what V, N and M would each give. No row carries a "pending Boss" mark.

**Stand-in kind** (from the gap plan): D derived, B bracket, T theory prior, P placeholder. **Quality** in plain words: *good* = checked against a source this project has read and inside its stated error; *fair* = derived or bracketed from read sources with one closure check; *weak* = recalled, placeholder, or sources disagree.

**Citation marks** (section 7 lists every reference with its mark). **V-rec** = the bibliographic record (authors, year, journal or arXiv id, title) matched in Crossref, arXiv or Open Library on 2026-10-17; the text was not opened. **V-read** = this thread opened the text and read the number cited. **PT** = a prior project thread marks the number READ in the project documents. **U** = unverified. Book entries name the chapter by topic only: no chapter number was opened, so chapter locators are unverified. Arithmetic marked COMPUTED was done in this thread.

**Free tools.** "Reachable" means the site answered HTTP 200 from this container (section 6). A reachable site can still refuse a bulk download or need a registration; I did not run a calculation on each.

## 1. What most needs Boss

1. **Giant-planet visible absorption (AD1, GEN.422, atmosphere colour).** The project names [R2265] [KT10] and [R2258] [Kar94] (both paid) as the sources of the visible methane and ammonia band strengths. The free line lists do not reach the visible: the HITEMP methane list covers 0 to 13,400 cm-1, which is wavelengths longer than 746 nm [Har20] (COMPUTED from the abstract: 619 nm is 16,155 cm-1 and 727 nm is 13,755 cm-1, both outside). Fitting the bands to the giants' own albedo would tune to the Solar System (the NO FIXED CONTENT rule). This is the gap with no free route that I found, and it is the one that most limits the colour of the giant planets.
2. **[R2726] [Lod03] (BC1, GEN.468).** The free copy is behind a captcha from this container, so a person has to fetch it. The water-share fill is derived and reproduces the 0.14 to 0.39 bracket, so nothing is blocked; the paper would add the condensation temperatures and check the derivation.
3. **"Most permission take wins" is settled.** Boss chose "widest range" by card at 1:12 PM EDT. The one place it changes an outcome today is the halo scatter (Dragonfly 1.01 dex against Harmsen 0.33 dex): under M the widest, 1.01, applies (section 5, item 8).
4. **Closed since the plan, no file needed any more.** BC13 (radiation weights and dose limits) is readable from the free IAEA GSR Part 3 [IAEA14]: only the neutron weighting function is missing. AD4 for SO2 and H2S no longer needs [R0443] [Bog78] if the project computes the cross section from polarizability and anisotropy (row AD4); the paper would only check it.
5. **Three citations in the project files look wrong** and are worth correcting before anyone builds on them (section 5, items 1, 3 and 5): "[R3734]" is [R0618]; "Yu and Ho 2018" for the bulge result is probably Yu and Ho 2019; and the fO2 buffer must not be cited to [R3365].

## 2. Twelve laws that cover most of the rows

For a quick approximation by hand, these are the models behind the tables. Inputs are what the generator already stores or can derive.

| # | Law or model | Form in words | Inputs | First place to look | Free tool |
|---|---|---|---|---|---|
| L1 | Hydrostatic equilibrium with an equation of state (stars, planets, neutron stars) | pressure gradient balances gravity; for a neutron star add the general-relativistic terms (Tolman-Oppenheimer-Volkoff) | mass, equation of state table | [KWW12] on stellar structure, [Poi00] and [SD08] on planetary interiors, [LP16] on neutron stars | CompOSE tables [Typ22], any ODE integrator; BurnMan [Cot14] for mantles |
| L2 | Virial theorem | 2 x kinetic + potential = 0, so sigma squared is about G M / (eta r); for a polytrope of index 3/2 the total energy is -(3/7) G M squared / R | mass, size | [BT08] on the virial theorem and dynamical friction | pencil and paper |
| L3 | Kennicutt-Schmidt law with gas continuity | star-formation rate per area = A x (gas per area) to the power N, with N = 1.4 +- 0.15 in [Ken98]; gas changes by inflow minus star formation | gas profile, inflow rate | [KE12] review, [Big08], [Pag09] on infall models | TNG data [Nel19] for inflow histories |
| L4 | Halo mass function, subhalo mass function and stellar-to-halo relation | number of haloes per mass from collapse theory ([ST99], [Tin08]); subhaloes follow a universal function of mass ratio [Gio08]; stars from [Mos13] | cosmology, mass | [MBW10] on halo statistics | Colossus [Die18], HMFcalc [Mur13] |
| L5 | Fuel consumption theorem | stars burn a fixed fraction of fuel per unit luminosity: number of stars in a phase = B x L_total x time spent, B = 1.7 +- 0.4e-11 per solar luminosity per year in [Bu97] (V-read) | luminosity, phase duration | [SC05], [Bu89] | MIST tracks [Cho16] |
| L6 | Rayleigh scattering with the King factor | cross section = 24 pi cubed / (wavelength to the 4th x N squared) x ((n squared - 1)/(n squared + 2)) squared x F_K; F_K = (6 + 3 delta)/(6 - 7 delta) in terms of the depolarization ratio delta [You81], equivalently F_K = 1 + 2 gamma squared / (9 alpha squared) for polarizability alpha and anisotropy gamma (COMPUTED algebra, this thread) | refractivity or polarizability, anisotropy | [BH83], [Bat84], [Buc95] | NIST CCCBDB for polarizability (page contents not opened, U) |
| L7 | Beer-Lambert absorption with correlated-k | transmission = exp(-column x cross section); correlated-k sorts line strengths to avoid line by line | line list, temperature, pressure | [GY89], [LO91], [Sea10] | HITRAN [Gor22], ExoMol [Ten20], petitRADTRANS (reachable) |
| L8 | Clausius-Clapeyron and Gibbs minimisation | vapour pressure rises as exp(-latent heat / (R T)); equilibrium composition is the minimum of total Gibbs energy | species, temperature, pressure | [Pie10], [Sea10], [LF02] | NIST WebBook, CoolProp [Bel14], MELTS, easyCHEM |
| L9 | Oxygen fugacity from the iron-wuestite identity | delta IW = 2 log10(a_FeO / a_Fe) | FeO and Fe activities | [Hir21] for the IW calibration, [Fro91] | MELTS; GEOROC for test compositions |
| L10 | Energy balance for sublimation | absorbed sunlight = thermal emission + latent heat x sublimation rate | distance, albedo, vapour pressure | [CA79], [FS09] | NIST WebBook phase data |
| L11 | Stoermer cutoff and effective dose | vertical cutoff rigidity scales with dipole moment over distance squared [SS05]; effective dose = sum over tissues of w_T times sum over radiation types of w_R x absorbed dose | field moment, spectrum, weights | [Att86], [Tur07] | EXPACS (PARMA [Sat15]), CREME, OLTARIS, IGRF (all reachable) |
| L12 | Photolysis and Troe falloff | photolysis rate = integral of cross section x quantum yield x flux; termolecular rate = k0[M]/(1 + k0[M]/kinf) x broadening factor | cross sections, flux, rate constants | [BS05], [Troe79], [Atk04] | IUPAC datasheets, MPI-Mainz atlas [KR13], Leiden database [Hea17], NIST kinetics |

## 3. Open items by area

Counts: galaxy 16 rows, stars 11, bodies and composition 18, atmosphere and colour 14 plus 11 groups of thin rows, radiation 1 (AD5 and BC13 are listed under atmosphere and bodies), comets 4. That is 64 open lines and 11 groups of thin rows. Ids are those of the gap plan (GR, GX, BC, AT, AD) or the TODO item (GEN.nnn).

### 3.1 Galaxy

| Id and what is missing | Nearest law or model, and its inputs | Where to look | Stand-in today | Rule |
|---|---|---|---|---|
| **GR1** Stellar halo mass (or accreted fraction) for M* 3e9 to 3e10. No survey reaches haloes this faint; the six measured haloes span 3.8e10 to 8e10 [Har17] | Accreted stars = stellar mass of the merged satellites. Satellites come from the universal unevolved subhalo mass function (dN/dln mu = N0 x^-a exp(-6.283 x^3), a 0.8, N0 0.21 [Gio08]), their stars from the stellar-to-halo relation [Mos13], dark below the reionization occupation gate [BLF20], and only mergers under a mass-ratio ceiling reach the halo [RG16] | [MBW10] on mergers and halo build-up; [BT08] on dynamical friction and tidal stripping; data: Illustris and TNG [Nel15] [Nel19], Magneticum accreted fractions [Rem22]; Colossus [Die18] | B: median anchored on the six haloes, plus the D integral. Integral gives 0.071 against 0.056 measured at 5e10 but over-accretes 2 to 3 times overall and its slope (1.4 to 2.2) misses 0.79 of [DS14]. Quality: fair at 5e10, weak below 3e10 | V. The Dragonfly scatter conflict (section 5, item 8) is the one place M differs |
| **GR3 left-over** Thick-to-thin LENGTH ratio against Vc: [YD06] falls with Vc, [Com18] median 1.46 rises. (Height ratio and mass ratio are Decided widen) | Two-component disc decomposition; hard bound: thick-disc mass cannot exceed the hot-component mass, log M_h = 0.489 + 0.9161 log M* [Com21] | [Com18] and its VizieR table J/A+A/610/A5; [Com21]; [YD06]; VizieR | D bound plus a B widened draw. Quality: good for the mass ratio on the mean relation (39 percent of galaxies exceed the bound one by one), fair for lengths | Not one of the six. V cannot decide (both verified). N takes [Com18]. M keeps the wider spread |
| **GR4** Age-velocity normalisation for galaxies other than the Milky Way | Disc heating: dispersion grows as age to a power, beta_z 0.441 +- 0.007 [SH20]; age-averaged vertical dispersion tracks circular speed, sigma_z,0 = (0.248 +- 0.038) V(2.2 h_R) falling as exp(-R/2h_R) [Mar13]; for flat star formation the age average of (t/tau)^beta is 1/(1+beta) | [BT08] on disc heating; [AB09] Table 6 for the age dependence of sigma_z/sigma_R (0.41 at 1 Gyr to 0.57 at 10 Gyr); [SH20]. Age-resolved kinematics of other galaxies: integral-field surveys (MaNGA, CALIFA), not opened (U) | D rule. Milky Way closure is 20 percent low (16.9 against 21.1 km/s). Quality: fair, one test case | V ([AB09] wins over the constant 0.54) |
| **GR6** Start of thick-disc star formation against host mass in star-forming low-mass discs. Ages exist for the Milky Way [Hay13] [Ben14] and three S0 galaxies [Pin19a] [Pin19b]; [Sat24] gives mean ages only | Thick-disc stars form in the gas-rich turbulent phase when the halo first collapses, so the start follows the halo formation epoch from the mass-accretion history (z50 1.27 to 0.35 over halo mass 1e8 to 1e15 [Cor15]) plus a delay | [MBW10] on formation epochs; [Com21] predicts the youngest thick-disc age against M*; arXiv for the [Sat24] star-formation-history table; Colossus [Die18] for redshift to time | T/P: halo formation epoch + 0 to 1 Gyr delay (Decided placeholder). Quality: weak | V (read ages win; [Sat24] is a mean) |
| **GR7** Radius ratio of the metal-poor (blue) to the metal-rich (red) cluster system against galaxy mass. Papers disagree: more extended against equal within errors | No theory gives a number. Reasoning (inference): blue clusters follow the accreted halo and red clusters the in-situ bulge, so the ratio would rise with the accreted fraction (GR1) | [FP12], [HR14], [Kar16]; VizieR cluster tables | B: log-uniform 1.0 to the largest printed ratio, 1.5 until read. Quality: weak | D (Decided widen, keep both) |
| **GR8** Radius of the cluster system for log M* 8.5 to 10.3 (spirals and intermediate masses are poorly sampled). Ratio to the galaxy's radius is about 3.5 [HR18] [For18] against about 6 [Kar16] | Cluster-system size scales with the size of the host dark-matter halo [HR18]; the halo virial radius follows from halo mass ([Mos13]) and the spherical-collapse overdensity ([MBW10]) | [HR18]; [HR14]; [Geo10]; Colossus [Die18] for halo radius from mass | B: interpolation with 0.3 dex (P). Quality: weak to fair | D (the 3.5 against 6 conflict is one of the six) |
| **GR9b** Cluster ages against host mass. Milky Way two-track data only | Cluster formation follows the host's star-formation history ([Beh13] average histories) with efficiency rising with gas pressure; old track from the infall times of accreted satellites (GR10) | [Pfe18] E-MOSAICS; [For18]; [Kru14] | B: 1 Gyr per track (P). Quality: weak | V |
| **GR9c** Nuclear cluster and nuclear disc ages and dispersions | Virial theorem: sigma squared about G M / (eta r_h); nuclear cluster mass log M = 0.48 log(M*/1e9) + 6.51, occupation about 80 percent near 1e9 [Neu20]; old base population plus recurrent young bursts above 1e9; in-spiral by dynamical friction [TOS75] [BT08] | [Neu20] (free review tables); [Geo16]; [BT08] | B: virial dispersion, young mass fraction uniform 0 to 0.3 above 1e9 (P). Quality: fair for dispersion, weak for ages | V (slope conflict is GX3) |
| **GR10** Stellar halo age against halo mass as a simulation fit (six measured haloes only) | Mean halo assembly history, z50 = 1.27, 1.15, 1.02, 0.88, 0.74, 0.60, 0.47, 0.35 at log M 8 to 15 [Cor15], then infall times of the accreted satellites | [CC10] (Aquarius, Milky Way mass only); [Cor15] (the authors' code, not opened, U); Colossus [Die18] | D analytic mean; the simulation fit is open. Quality: fair (mean, no scatter) | V |
| **GR11** Dwarf stellar-to-halo relation and scatter for 1e8.5 to 1e10.5; abundance-matching fits stop at 1e10.5 | Abundance matching with the steeper low-mass slope: [Beh13] Eq. 3 with slope 1.92 [GK14]; scatter bracket 0.19 dex [Nad20] to 0.5 dex [Fit17] | [GK14]; [GK17]; [Nad20]; Local Group dwarf masses through VizieR | B: median extension flagged, scatter uniform 0.19 to 0.5. Quality: weak below 1e9 | Already the widest bracket (M). N would take [Nad20], newer than [Fit17] |
| **GR14** Accreted outer envelope above M* 1e11. No method isolates an accreted halo in bulge-dominated galaxies; [DS14] reaches 11.4 only as an outer-component proxy | Ex-situ fraction rises toward 1: Illustris-1 medians 0.12 at 4.7e10, 0.19 at 1e11, 0.34 at 2e11, 0.54 at 4.7e11, 0.72 at 1e12 [RG16]; Magneticum has its own curve [Rem22] | [RG16] Fig 3 (digitised table in the Gemini-check files); [Rem22]; TNG data [Nel19] | B: logistic, centre bracket 5e10 to 1e12 (P). Quality: weak | M (two simulation curves, both in question: widest centre) |
| **GX1** Halo mass given M* (the conditional, not the inverse of the median). A derivation to code, not a data gap | Bayes: p(M_h given M*) proportional to p(M* given M_h) times the halo mass function; the scatter of [Mos13] shifts the median (Eddington bias) | [ST99]; [Tin08]; [MBW10]; Colossus [Die18] and HMFcalc [Mur13] | D. Quality: good once coded | V |
| **GX2** Arm pitch against bulge fraction: the sign differs between [Har17b] and the Yu and Ho result | Density-wave theory: pitch angle is set by the shear rate and disc stability | [BT08] on spiral structure; [YH19]; VizieR | B: sign 50/50, magnitude from each paper. Quality: weak | D. Citation fix: the project's "Yu and Ho 2018" is probably [YH19]; [YH18] is about wavelength dependence |
| **GX3** Nuclear cluster mass slope: 0.48 [Neu20] against 1.0 to 1.36 [Geo16] | Two formation channels, in-spiral of clusters by dynamical friction [TOS75] and in-situ star formation (discussed in [Neu20]); which channel gives which slope is not settled | [Neu20]; [Geo16] | B: uniform 0.48 to 1.36. Quality: weak | D |
| **GX4** Environment (field, group, cluster) for the 0.3 dex cluster-count factor | Halo mass function gives the group probability of a halo | [MBW10]; Colossus [Die18] | T: environment drawn from halo mass, no new concept | D (Decided) |
| **GX5** Halo metallicity baseline (-0.6 against -0.9) and its law outside halo mass 5e8 to 1e10 | Halo [Fe/H] = mass-metallicity relation of the dominant progenitor minus 0.15 [DB18]; dwarfs from the universal dwarf relation [Kir13], massive galaxies from [Tre04]; infall chemical evolution [Pag09] | [Kir13]; [Tre04]; VizieR | D with 0.2 dex (P). Quality: fair inside 5e8 to 1e10, weak outside | V |

### 3.2 Stars

| Id and what is missing | Nearest law or model, and its inputs | Where to look | Stand-in today | Rule |
|---|---|---|---|---|
| **GEN.463** Neutron star radius: the fixed range 10 to 13 km is replaced by a derived radius (the sourced radii give 11.6 to 13.1 km) | Tolman-Oppenheimer-Volkoff equations plus a cold-matter equation of state p(energy density) give R(M); constraints: causality, tidal deformability from GW170817 [Abb18], radii of PSR J0030 [Mil19] [Ril19] and J0740 [Mil21] [Ril21] | [LP16] (R(M) bands); [Rea09] (piecewise-polytrope parametrisation); [Typ22] and the CompOSE tables; a TOV integrator is about 20 lines | T/D: radius from mass and the equation-of-state band. Quality: good (bands from verified papers; the equation-of-state choice is the knob) | V |
| **GEN.467** Main-sequence lifetime: the fixed 10 Gyr per solar mass and power law -2.5 are replaced by relations from mass and metallicity | Fuel consumption: lifetime = (fraction of hydrogen burned x mass x hydrogen fraction x 0.007 c squared) / mean luminosity; 0.007 is the mass fraction released by fusing hydrogen to helium (about 0.7 percent, from the recalled atomic masses of H-1 and He-4; NIST atomic weights are reachable to check) | [KWW12]; [SC05]; analytic formulae t_MS(M, Z) in [HPT00]; tracks [Cho16] [Dot16] (MIST, reachable) | D from [HPT00] or a track grid. Quality: good if checked against a track | V |
| **GEN.469** Disc star-formation history: the local B-star share is 0.2 percent against 0.04 percent observed (the observed value stays a test, never an input) | Kennicutt-Schmidt law [Ken98] [KE12] [Big08] with gas continuity (inflow minus star formation); inflow set by halo accretion [Cor15] or taken from TNG [Nel19]; the B-star share is the stars younger than a few hundred Myr over all stars, times the IMF's B fraction | [Pag09] on infall models; [KE12]; TNG data [Nel19] | D, with the inflow as the missing input. Quality: fair | V |
| **GEN.450** Radial slope of the star-formation profile: the Schmidt index 1.12 to 1.2 gives -0.19 dex per kpc, the observed -0.28 needs about 1.65 | Same law; ways to close the gap: an index that rises at low gas density, a gas disc more extended than the stars, an outer-disc threshold | [Ken98] (1.4 +- 0.15); [Big08]; [KE12]; HI and CO radial profiles through VizieR | D, the choice to be named. Quality: fair | V ([Ken98] verified; the observed gradient is a test) |
| **GEN.529** Pre-main-sequence track: Hayashi temperature, birthline radius and energy factor are recalled ([I]), not looked up | Virial theorem and Kelvin-Helmholtz contraction: t = integral of (energy factor) G M squared / (R squared L) dR, energy factor 3/7 for a polytrope of index 3/2 and 3/4 for index 3 (COMPUTED from the polytrope result); the Hayashi line is the fully convective limit; the birthline comes from deuterium burning during accretion [Sta83] | [KWW12]; table lookup in [Bar15] (CDS tables), [Sie00] (online server), MIST pre-main-sequence tracks [Cho16] | [I] constants; the 12 Myr derived time against the 40 Myr anchor shows the factor is wrong. Fix is a table lookup (D). Quality: weak today | V |
| **GEN.493** Initial-rotation distribution and saturation time against stellar mass (UV regime of the host) | Rotation-activity relation: activity constant for Rossby number below saturation, falling as Rossby to about the -2.2 beyond [Wri11]; spin-down t^-1/2 [Sku72]; rotation tracks [Joh21] | [Wri11]; [Jac12]; [Med22]; cluster rotation periods through VizieR | T/B: tau_sat 0.16 Gyr (0.5 to 0.6 Msun), 0.6 +- 0.3 (0.2 to 0.3), 2.3 +- 1.3 (0.1 to 0.2) (PT), and a percentile u on the [Joh21] tracks. Quality: fair | V (Rossby saturation values of different authors are never mixed, 0.13 to 0.16 with Wright's times) |
| **GEN.494** Brown-dwarf spectral type against temperature; the T/Y boundary | Type is a monotone function of Teff with scatter; Teff from cooling curves | [Kir21] Table 13 (RMS 134, 79, 79 K; Y0 at 371 to 510 K); [BL93]; ATMO 2020 [Phi20] and [Bar03] grids (reachable) | D/B from [Kir21]. Quality: fair | V ([Kir21] wins over "T/Y gap 500 to 700 K", which is T8 to T9) |
| **GEN.502** Giant branch: red-giant time law against core mass, tip by [Fe/H], branch duration against the horizontal-branch ratio | Core-mass-luminosity relation of hydrogen-shell burning with time = core growth / (L / (X x energy per mass)); tip set by helium ignition [Ser17]; fuel consumption [Bu97] | [SC05]; [KWW12]; fit of [Ree25] (slope -0.83 at the base to -0.56 at the tip); [Ser17] tip 1,650 to 1,780 Lsun at Z 1e-4 | B/T: one exponent 0.78, tip 2,500, ratio 0.85 against 1.45 (measured mean 1.47 +- 0.025) [Sal04]; the 1.44 +- 0.06 once quoted is not in the paper's text. Decided (Boss, contradiction sweep card, 7:29 AM EDT 2026-10-21): the ratio is 1.45. Quality: fair (within 15 percent) | V (published numbers win over recalled) |
| **Stellar tables, leftover [R] values** Horizontal-branch colour row at [Fe/H] -1.5 (0.25 too red), extreme-HB luminosity (code 40 to 60, source 10 to 40), age shift 0.1 dex per Gyr (not a published number), AGB life 0.002, RR Lyrae edges | HB position is set by core mass and envelope mass lost on the giant branch (Reimers law [Rei75], [SC05b]); second parameter [LDZ94] | [SC05]; [Dot10]; MIST and BaSTI HB tracks (reachable) | P, tagged [R]. Quality: weak | V |
| **Evolutionary flux** 2.2e-11 per solar luminosity per year in the files | Fuel consumption theorem constant B | [Bu97] (V-read: 1.7 +- 0.4e-11), [Bu89], [RB86] (not opened) | Replace by 1.7 +- 0.4e-11. Quality: good | V and N agree (the read, newer paper wins) |
| **BC19 / T3** White-dwarf spectral type (DA fraction) against temperature | Spectral evolution: settling and convective mixing change the surface composition with cooling | [BRL97]; catalogue fits [GBB19] [Kep19]; the Montreal White Dwarf Database and VizieR (reachable) | B: binomial about a fitted fraction, downloaded at build time (Decided). Quality: fair, not yet fitted | V |

### 3.3 Bodies and composition

| Id and what is missing | Nearest law or model, and its inputs | Where to look | Stand-in today | Rule |
|---|---|---|---|---|
| **BC1** (GEN.468) Water share of solids beyond the snow line, and condensation temperatures | Oxygen budget: oxygen goes to MgO, SiO2, Al2O3, CaO, Na2O and FeO first, carbon takes CO, CH4 or CO2, what remains is water. Solar mix gives 0.17 to 0.32 (0.40 at C/O 0.35); all-methane limit 0.54. Then equilibrium condensation: ice forms when the partial pressure passes the vapour pressure (Clausius-Clapeyron) near 150 K at nebular pressure | [Lod03] (50 percent condensation temperatures; free copy behind a captcha); [WSH19]; [LF98]; abundances [Asp21] Table 2; NIST WebBook phase data, MELTS (reachable) | D (the plan's oxygen budget, COMPUTED; reproduces 0.14 to 0.39). Quality: good for the bracket, no disc chemistry | V. The "keep the switch" of 0.14 to 0.39 against about 0.5 is Decided; the derived function spans both |
| **BC2** Stellar Mg, Si, O, C against [Fe/H], with the quality cut, and thick-disc and halo stars | Galactic chemical evolution: alpha-element ratios fall as iron rises because type Ia supernovae add iron later [Pag09]. Regression [X/Fe] = a + b [Fe/H] from [Bre16] (C slope -0.275, O -0.446, Mg -0.247, Si -0.196; residual sd 0.06 to 0.13) | [Pag09]; [Bre16] (VizieR J/ApJS/225/32); Hypatia catalogue (reachable, licence unstated) | D from DATA (1592 stars, no quality cut applied; the cut must be applied before use). Quality: good for F, G, K disc stars, [Fe/H] -1.2 to +0.6 | V |
| **BC3** Mantle mineral equations of state: the files carry the 2024 Part III parameter set, not Part II | Mie-Grueneisen-Debye thermal equation of state on a Birch-Murnaghan cold curve; phase equilibria by Gibbs minimisation | [SLB05]; [SLB11] Table A1 (PT, 100 GPa density test passes); [Poi00]; [SD08]; BurnMan [Cot14] (install route not tested, U) | B: second-hand set. Quality: good (densities agree within 0.05 to 0.6 percent at 100 GPa) | V (the paper's Table A1 wins on q0, theta0, post-perovskite K0). Not a blocker |
| **BC4** Oxides and K, Th, U of the mantle: 2014 edition values (Th 0.0849, U 0.0229 ppm) unconfirmed | Bulk silicate Earth from cosmochemical constraints: refractory lithophiles at chondritic ratios; crust from a compilation | [PO14] (2014 not opened; 2003 edition read: oxides total 98.41, K 260 ppm, Th 83.4 ppb, U 21.8 ppb); [MS95]; [RG14] (crust, matches every file value); GEOROC and EarthChem (reachable) | B: 2003 values, labelled "2003 edition". Quality: good for the crust, fair for the mantle | V (2003, read, over 2014, unread). If the 2014 edition is read, two verified editions disagree and N takes 2014 |
| **BC5** Giant-planet interiors: core mass, diluted core, inflation efficiency | Hydrostatic structure with a hydrogen-helium equation of state [SCvH95] and a heavy-element mass M_Z = (57.9 +- 7.03) M^(0.61 +- 0.08), M_Z in Earth masses and M in Jupiter masses, 47 planets, scatter a factor 1.82 +- 0.09 [Tho16] (V-read) | [For07]; [SCvH95]; [Tho16]; [Cha25] (147 giants, 1.4 to 1.6 times lower than [Tho16] between 0.3 and 3 Jupiter masses) | T: Dirichlet core fraction between the two measured giants as outer limits. Quality: fair for bulk, weak for core structure | Two verified fits disagree: N takes [Cha25] (2025, newer than 2016); M keeps both, the union; V cannot decide. Not decided |
| **BC6** Moon composition against distance and origin; zero confirmed exomoons | Circumplanetary-disc ice line r_ice = 20 R_J (M_p/M_J)^(2/3) (tau/5 Myr)^(-1/3) [Sas10]; satellite mass fraction about 1e-4 of the planet [CW06] [CW02]; population synthesis [Cil18] (85 percent icy, Galilean-like in 20 percent) | [CW02]; [CW06]; [Sas10]; [Cil18] | T: ramp with distance in snow-line units, bracket 0.5 to 0.85 icy. Quality: weak | Two verified sources disagree ([Sas10] against [Cil18]), not one of the six. N takes [Cil18]; M keeps the bracket (as today) |
| **BC7** Ice densities and the ice compaction law | Pressure compaction of porous ice: density factor exp(-b P) with a floor [BN19]; solid density of an ice from its crystal cell: rho = Z M / (N_A V) | [BN19] (Table 1: ice 920, silicate 3500 kg/m3, porosity 0.60, floor 0.20, b1 -0.1 per MPa, b2 -0.11); [FS09] (paid); NIST WebBook, PubChem, COD crystal cells (all reachable) | B: [BN19] read. Quality: good for water ice, fair for other ices | V (the paper wins: reference pressure 10 MPa, not 1e6 Pa) |
| **BC8** Planet-minus-star scatter in Fe/Mg and the FeO against star and water | Planet composition follows the star's refractory ratios with a slope near 1 and giant-impact scatter | [PV26] (slope 0.94, band +1.02 to -1.07); [PV20]; [Adi21]; NASA Exoplanet Archive and DACE (reachable) | B: printed band as width. Quality: weak (few precise planets) | Two verified papers by the same group ([PV20] older, [PV26] newer): N takes [PV26]; V indifferent; M widest |
| **BC9** Water equation of state at high pressure (superionic) | Free-energy equation of state; ice-phase Gibbs surfaces; ab initio tables above about 100 GPa | [WP02] (IAPWS-95); [Fre09]; [Jou20]; AQUA [Hal20]; IAPWS and CoolProp [Bel14] (reachable) | B between ice VII, ice X and the superionic range, 5 percent. Quality: fair | V |
| **BC10** Albedo of weathered regolith and frost; grain state | Hapke radiative transfer: single-scattering albedo to reflectance, then phase integral to Bond albedo | [Hap12]; spectra from RELAB (reachable) and the USGS library [Kok17] (record only; the page was not reachable) | B: laboratory reflectance bracket, 0.05 in albedo (P). Quality: weak | V |
| **BC11** Volatile stripping when a planet is ejected | Energy-limited escape during ejection; a retained thick hydrogen atmosphere can insulate the surface [Ste99] | [Ste99]; [Pie10] | P: uniform 0 to 1 loss fraction (Decided). Quality: weak | M (no measurement; widest) and Decided |
| **BC12** Vapour-pressure validity windows and the cold end | Clausius-Clapeyron with temperature-dependent latent heat; Wagner or Antoine forms inside each fit's window | [FS09] (paid; HAL record, U); NIST WebBook phase data; [SW96]; CoolProp [Bel14] | B: each fit inside its window, Clausius-Clapeyron beyond, width doubled. Quality: fair | V |
| **BC13** Radiation weights and dose limits | Effective dose = sum over tissues of w_T times sum over radiation types of w_R times absorbed dose | **Read in this thread** from [IAEA14]: w_R 1 for photons, electrons and muons; 2 for protons and charged pions; 20 for alpha particles and heavy ions; neutrons are a continuous function of energy (blank in the extract). w_T 0.12 (red marrow, colon, lung, stomach, breast, remainder), 0.08 gonads, 0.04 (bladder, oesophagus, liver, thyroid), 0.01 (bone surface, brain, salivary glands, skin). Limits: 20 mSv per year averaged over five years, 50 mSv in any one year; the glossary cites [ICRP07] | D/B from [IAEA14]. Quality: good except the neutron function | V (IAEA copy read; ICRP 103 is its cited source, not opened) |
| **BC14** (GEN.473, 477) Belt composition kernel | Condensation zoning: rocky inside the snow line, water-rich beyond; taxon mix against heliocentric distance in one belt [DC13]; early giant-planet migration moves material [Wal11] [RI17] | BC1 for the snow line; [DC13] Table 6; [Kro14] | T: zoning boundary at the star's snow line, width 0.3 in log AU (P); 4 of 84 names (Decided). Quality: weak | M (one belt only, extrasolar ensemble absent) |
| **BC15** (GEN.474) Rock types from bulk oxides | CIPW-style norm: allocate oxides to normative minerals by mass balance; classify with the total alkali-silica diagram [LB86]; phase equilibria by Gibbs minimisation | [PA09]; [LB86]; MELTS; GEOROC for test compositions (reachable) | D, deterministic. Quality: good once a known rock reproduces its mineral list | V |
| **BC16** (GEN.475) Meteorite groups as outputs | Groups are regions in (Fe/Si, Mg/Si, oxidation state, oxygen isotopes) | [Kro14]; Meteoritical Bulletin Database (not reachable from this container) | D from BC1, BC2, BC15. Quality: fair | V |
| **BC17** (GEN.476) Al-26 and Fe-60 from the birth environment | Radioactive decay N = N0 2^(-t/half-life); supply from supernova and Wolf-Rayet yields with the odds of a nearby massive star from the embedded cluster mass function | [Lic16]; [Die06] (galactic Al-26 from gamma rays); [Rug09] (Fe-60 half-life, number not opened); NuDat (reachable) for half-lives | T: Bernoulli with cluster-size odds; evidence points to a mixture (an enriched flag). Quality: weak | V |
| **BC18 / T7** Core light elements | Metal-silicate partitioning as a function of pressure and oxygen fugacity; density deficit against pure iron | [Fis15]; [Hak18]; [DA81] (Earth only) | T: Dirichlet over S, Si, O, H (P). Quality: weak | M (cores unobservable; widest) |

### 3.4 Atmosphere and colour

| Id and what is missing | Nearest law or model, and its inputs | Where to look | Stand-in today | Rule |
|---|---|---|---|---|
| **AD1** (GEN.422) Visible absorption of CH4, NH3 and H2O in giants | Beer-Lambert absorption with a band cross section per molecule; line-by-line from a line list, or an empirical band-strength table | [KT10] and [Kar94] (paid) hold the visible methane bands; line lists [Gor22] (HITRAN), [Har20] (HITEMP CH4, longer than 746 nm), [Yur14] [Ten20] (ExoMol) do not reach 619 nm or 543 nm; theory: [GY89], [Sea10] | B: band strengths from neighbouring bands, 10 percent (P), data downloaded at build time. Quality: weak in the visible; fair beyond 750 nm | V. Fitting to the giants' own albedo is not allowed (NO FIXED CONTENT) |
| **AD2** (GEN.423) Mars dust and Titan haze optical spectra | Mie theory from the complex refractive index (n, k) gives extinction, single-scattering albedo and phase function [BH83] [HT74]; Mars dust single-scattering albedo 0.922 to 0.926 in the solar band for effective radius 1.2 to 1.8 micron [Wol09] (figures only, digitised) | [Wol09] Figs 12 to 14; tholin n and k [Khare84]; Titan haze model [Tom08] (paid); RELAB (reachable); any free Mie code | B: similar materials, 10 percent (P). Digitised values labelled "digitized from Figures 12 to 14", never "the paper's table". Quality: fair for Mars, weak for Titan | V |
| **AD3** (GEN.413, 414) Infrared opacity (correlated-k) tables; the greenhouse forcing laws are already read [Myh98] | Correlated-k from line-by-line cross sections [LO91]; analytic non-grey temperature profile [PG14] [Gui10]; mean opacities [Fre14] | [GY89]; [Pie10]; HITRAN, HITEMP, ExoMol; petitRADTRANS (reachable) | D: computed at build time, nothing stored (Decided). Quality: good once computed | V |
| **AD4** (GEN.422) Rayleigh cross section for SO2, H2S, H2, He, NH3 and water vapour (air, N2, O2, Ar, CO2 are filled) | L6: cross section from refractivity and the King factor, or from polarizability alpha and anisotropy gamma: sigma = 128 pi to the 5th / (3 wavelength to the 4th) x alpha squared x (1 + 2 gamma squared / (9 alpha squared)), alpha in volume units (COMPUTED algebra). The depolarization ratio of SO2 and H2S is the missing number in the files | [Bog78] (paid) holds the SO2 and H2S depolarization; [Kin23] origin; [You81] [Bat84] [Buc95] for the working forms; [SU05] measured 532 nm cross sections; NIST CCCBDB polarizability (page not opened, U); any free quantum-chemistry package gives alpha and gamma (not run here, U) | B: 5 percent (P). Generic error 2.5 to 4 percent for Ar and N2, 7 for CO2, 10 to 18 for CO, CH4 and N2O [SU05]. Quality: weak for these gases | V ([SU05] measured errors win over the generic 5 to 9 percent) |
| **AD5** (GEN.436, 444) Energetic particle spectra and shower runs beyond Mars (Mars surface dose is filled [Has14]) | Galactic cosmic ray spectrum [BO96] with force-field modulation [GA68]; solar proton event fluence probabilities [Xap00]; atmospheric showers by an analytic model [Sat15] | EXPACS, CREME, OLTARIS (reachable); NIST ESTAR, PSTAR and XCOM (reachable); [Att86] [Tur07] | B: 10 percent (P), downloaded at build time. Quality: fair. The Gale column is 21 g/cm2 against 16.4 in the code | V ([Has14] wins) |
| **AD6** (GEN.485, 499) Hot-giant line lists (TiO, VO, FeH) and vapour pressures | L7 plus L8: line-by-line absorption, and condensation from Gibbs equilibrium; silicate vapour pressure fits [SF09] [Kit16] | [Ten20] ExoMol (reachable) | D: downloaded at build time (Decided). Quality: good once downloaded | V |
| **AD7** (GEN.471, 491) Sulfuric acid and water: acid vapour pressure and brines | Acid pressure over a solution = mole fraction x activity coefficient x pure-acid pressure, activities from [Zel91]; water pressure from [Tab97] Eq. 1 and Table 1 (185 to 260 K, read) | [Aye80] (paid) for pure acid; [KL90] form is in the files; [Gia60]; E-AIM [WC02] (not reachable from this container) | B: 10 percent in log pressure (P). Quality: good for the water side, weak for the acid side | V. Printing faults: [Tab97] Table 2 is discontinuous at water activity 0.05 and 0.85; [Zel91] Eqs 11 and 12 have a sign and an index fault |
| **AD8** (GEN.496, 501) Photochemistry rates (JPL 19-5 not in hand) | L12: Troe falloff, Arrhenius, photolysis integrals | [Atk04] IUPAC datasheets (reachable); NIST kinetics (reachable); [KR13]; [Hea17]; [HM15]; [JPL19] (U); UMIST [Mil24] and KIDA [Wak12] are not reachable from this container | B: factor 2 (P), downloaded at build time. Quality: fair | V |
| **AT.a1, a2 (T8)** Mantle redox against mass; the later-mantle law | L9: delta IW = 2 log10(a_FeO/a_Fe); ideal-solution limit gives -2.5 to -1.6 for FeO 5 to 15 mole percent (plan 4.8); the 22 white-dwarf values give -1.25 +- 0.10 [Doy20] [Doy23] | [Hir21] (IW calibration); [Fro91]; [Sos20]; [Gai22]; [GS14]; [Hir12]; MELTS (reachable); GEOROC | D: identity with the white-dwarf zero point. Quality: fair at the small-body end; the Earth end and slope are Solar System only (weak) | V (Figure 2 of the Gaillard paper in Boss's upload, about FMQ-2.75, wins over its text's FMQ-2 +- 0.5; the report identifies it as [Gai22]). The white-dwarf scale (core formation) differs from the surface magma-ocean scale and needs a conversion test |
| **AT.b4** Water-rich share (published shares run from at most 3 percent to 87 percent by formation model) | Water delivered by migration from beyond the snow line: logistic in mass, centre 6.1 Earth masses, slope 6.3 per dex (COMPUTED by the tunable-gaps lane from NGPPS) [Ems21]; density modes [LP22] | [Ems21]; DACE and NASA Exoplanet Archive (reachable); the NGPPS Zenodo record is not reachable from this container | T: uniform 0 to 0.87 (P). Quality: weak | Two papers disagree ([R3801], at most 3 percent; [R0715], 59 percent water-rich; neither opened here, U). If verified, N takes the 2026 paper and M the union (today's bracket) |
| **AT.b5** Sub-percent water fraction | No ensemble measures it | none | P: log-uniform (Decided placeholder). Quality: weak | M and Decided |
| **AT.n7** Haze yield | Haze mass production rate = yield x methane photolysis rate (inference); laboratory rates span 250 times [Hor18]; Titan [Tra06] | [Hor18]; [Tra06]; [AM01]; [Roo22] (virga) | P: log-uniform (Decided placeholder). Quality: weak | M and Decided |
| **AT.j2, j3, j4** Mars-, Titan- and Pluto-like occurrence | These are outputs of escape and retention, not a separate frequency: atmosphere kept above a line in escape velocity and flux (cosmic shoreline) [ZC17], fitted logistically [BT25] and with slope 5.77 [MG25] | [ZC17]; [BT25]; [MG25]; [CK17] | T: derive from inventory, escape and vapour pressure; leftover prior [0, 1] flagged. Quality: weak | V (both 2025 fits beat the 13-body repo window). The two fits are the same year; the lane chose [BT25] |
| **AT.m3, m4, m11** N2O, PH3, NO2, CH3Cl presence; no measurement on any exoplanet | Abiotic yields: lightning makes NO from N2 + O2 at high temperature; photochemistry and rate constants [Atk04]; abiotic ozone [DG14] | [Atk04]; [DG14]; [BS05]; TUV and NIST kinetics (reachable) | T/X: flagged bracket [0, 1] (m4 [0, 0.05]), default off, never tuned. Quality: weak | M |

**Thin F rows** (the tunable-gaps note counts 14; 19 ids when d1 to d5 are listed singly). Each is tunable only to a few stated constraints.

| Rows | What is thin | Nearest law and the number to hold | Where to look | Rule |
|---|---|---|---|---|
| b1, b2, b3 | Water delivery payloads (dry, low, high: 0, 0.05, 0.5 weight percent, Earth-only scenarios) | Migration-delivery logistic in mass, centre 6.1 Earth masses [Ems21] | [Ems21]; [LP22] | M (Earth-only payloads) |
| d1 to d5 | Atmosphere retention classes and the K window (13 Solar System bodies) | Logistic classification: log10 f0 2.62, p 5.9, q 1.17, width 0.247 dex [BT25]; slope 5.77, intercept -4.35 [MG25] | [BT25]; [MG25]; [ZC17] | V |
| e5 | Envelope metallicity upper edge 0.6 (about 106 times solar) excludes K2-18 b | Atmospheric metallicity is not bulk envelope metallicity; retrievals give 0.6 to 162 times solar for K2-18 b [Ben19] | [Ben19]; NASA Exoplanet Archive | V |
| f7 | Giant C/O tail | Stars: median 0.47, sd 0.07, maximum 0.66, none above 0.8 [Bre16]; planets: core 0.57, 0.043 to 0.087 above 0.8 | [Bre16]; [LF02] | V (planet data sit outside the design tail 0.0013 to 0.34) |
| h4 | TiO and VO edges | Titanium kept only above about 2350 to 2640 K; VO kept at 2210 K; Gibbs equilibrium | [Ten20]; BC12 methods | V |
| i6 | Lava-world volatile fraction | 2 of 2 observed molten-dayside planets are volatile-bearing; Beta(3,1), median 0.79 | [Kit16]; [SF09] | V |
| k1, k2, k3 | Coma class shares | Counts 13, 9, 3 of 25 [HP22] (see 3.6) | [HP22] | M |
| n6 | Silicate-to-haze crossover near 950 K | Compute from the chemistry module; slope 0.0008 +- 0.00016 per K | [Gao20] | V |
| o2 | Thorium scatter in non-thin-disc hosts | 0.056 dex (solar twins) to 0.19 dex (mixed catalogues); Hypatia | Hypatia (reachable) | M (two ranges, both in question) |
| o3 | K/U and K depletion | Stellar K exists; planetary K depletion only from Solar System bodies; sodium-to-magnesium in white dwarfs as a proxy | tunable-gaps note | M |
| o6 | Al-26 present | Evidence points to a mixture (an enriched flag), not log-uniform | [Lic16]; see BC17 | V |

### 3.5 Radiation

| Id and what is missing | Nearest law or model, and its inputs | Where to look | Stand-in today | Rule |
|---|---|---|---|---|
| **GEN.453** Radiation constants: cutoff scaling, dose against column, field factor | Stoermer cutoff: vertical rigidity proportional to dipole moment over distance squared [SS05] (Earth about 14.3 GV; the code has 14.9); galactic cosmic ray spectrum [BO96] modulated by the force-field model [GA68]; shower dose against depth [Sat15] | EXPACS, CREME, OLTARIS, IGRF (reachable); [Att86] | Constants fitted to Earth and the Moon: about 500 mSv per year at the Moon and 255 at 22.5 g/cm2 are sound, 1.5 to 2.7 times low at 100 g/cm2. Quality: fair at two anchors, weak at 100 g/cm2 | V (shower simulations and [Has14] beat the flat floor 0.6) |

Also open and listed above: AD5 (energetic particle spectra and shower runs, in 3.4) and BC13 (radiation weights and limits, in 3.3, closed except the neutron weighting function).

### 3.6 Comets

| Id and what is missing | Nearest law or model, and its inputs | Where to look | Stand-in today | Rule |
|---|---|---|---|---|
| **AT.k7** (GEN.488, 495) Single-comet species N2, O2, Ar, Kr, Xe, SO2 (67P only; SO2 0.2 percent of water in Hale-Bopp [MC11]; NH3 about 1.5 percent, so the repo's 3.6 percent upper bound is not supported) | Gas trapped in amorphous water ice as it forms: trapped fraction depends on deposition temperature and partial pressure [BN85] [BN87] [BN88]; CO is enriched over N2 [NB96]; grains formed near 25 K [NB05]; the low N2 content of 67P points to a low formation temperature [Rub15] | [BB17]; [MC11]; [CII]; 67P: O2 [Bie15], Ar [Bal15], Kr and Xe [Rub18]; [BNO07] | T: ordering by volatility, Dirichlet in that order (P). Quality: weak | M (one comet) |
| **AT.k1 to k5, T16** Coma classes and abundances | Sublimation by volatility and heliocentric distance (L10): every coma beyond about 3.5 AU is CO-dominated, so class is log10(CO/CO2) at the comet's own distance; CO2/H2O geometric mean 10.9 percent (0.51 dex, N=21), CO/H2O 2.5 percent (0.65 dex, N=16) [Rob24]; medians 12 +- 2 and 3 +- 1 percent [HP22]; the two interstellar objects lie 1.8 to 3.6 sigma outside the Solar System fits | [HP22]; [DR16]; [BB17]; [CA79] for the energy balance | F: counts inside binomial intervals; per-system offset up to about 1.5 dex. Quality: fair within the Solar System, weak for other systems | M (the within-system scatter is Solar System only) |
| **T15** Cometary organic fraction (Marboeuf runs two extreme models; no probability exists) | Dust is carbon-rich in 67P [Bar17]; composition of volatiles and origin [AH12] | [Bar17]; [AH12]; [MC11] | P: organic fraction N(0.60, 0.15) clipped to [0, 1] (the lane's proposal). Quality: weak | M |
| **T19** Rock fraction of dwarf planets and Kuiper Belt objects | Compaction by self-gravity: porosity falls with size [BN19]; density gives the rock fraction | [BN19] (a 70 percent rock mass fraction fits 14 of 17 objects); dwarf planets 0.74 to 0.82, small objects 0.3 to 0.45 (COMPUTED by the lane) | P/B: dwarf planets 0.78 +- 0.12, porosity transition near 350 km. Quality: weak | Two read sources read the densities differently (rock fraction rising with size, against a constant 70 percent with porosity); dates not recorded here, so N cannot pick; M keeps the draw |

## 4. Closed since the plan (not listed as open)

| Item | What closed it | Where |
|---|---|---|
| GR2 curves; GR12 reionization spread; GR13 ultra-faint sizes; GR15 width | Digitised tables: Illustris ex-situ curve, D'Souza outer-component fractions, occupation fraction (centre 10^9.33 Msun/h, width 0.21 dex), centre against reionization redshift, reionization spread 6.3 to 11.6 (median 7.7), 15 ultra-faint sizes; Graus Eq. 1 (V50 18 km/s, sigma 2.5) agrees on the width | gemini-batch1-check section 4; boss-papers 4.10 |
| GR5 Courteau scatter | 0.097 dex about the paper's own line (79 markers, Fig 2) | gemini-batch1-check |
| GR9a cluster size | 2.7 +- 0.35 pc, red 17 percent smaller than blue [Jor05] | gap plan 4.6 |
| GR4 and GR6 for the Milky Way | [SH20] Table 2, [AB09] Table 6, [Hay13] | boss-papers |
| BC2 star abundances | [Bre16] regressions | gap plan 4.7 |
| BC13 weights and limits | [IAEA14], read in this thread (neutron function excepted) | row BC13 |
| AD3 forcing laws | [Myh98] Table 3 | boss-papers 2.4 |
| AD4 for air, N2, O2, Ar, CO2 (and CO, CH4, N2O, SF6 at 532 nm) | [You81] [Bat84] [Buc95] [SU05]; Earth tau_R test value 0.09721 | boss-papers 2.1 |
| AD5 Mars surface dose | [Has14] | boss-papers 2.5 |
| AD7 water side and activities | [Tab97] [Zel91] | boss-papers 2.3 |
| BC3 Part II, BC4 crust and first-edition mantle, BC7 compaction | [SLB11], [RG14] with the 2003 Palme and O'Neill values, [BN19] | boss-papers 2.8 to 2.10 |
| Stellar tables: HB-to-RGB ratio and evolutionary flux | 1.45 (measured mean 1.47 +- 0.025) matches 1.4 [Sal04] (Decided (Boss, contradiction sweep card, 7:29 AM EDT 2026-10-21): 1.45); flux 1.7 +- 0.4e-11 [Bu97] (V-read in this thread) | stellar-evolution-tables-check |
| Giant mass-metallicity law | [Tho16] V-read: M_Z = (57.9 +- 7.03) M^(0.61 +- 0.08) | row BC5 |

## 5. Conflicts and filing errors, with the rule that decides each

**Errors in the project's own records** (a fix for the TODO thread; each was checked in Crossref or arXiv on 2026-10-17):

1. **"[R3734]" (astro-ph/9711072) is [R0618]**, single author, Evolutionary population synthesis models of primeval galaxies: a critical appraisal [Bu97]. Rule V (the paper read shows the author).
2. **Evolutionary flux 2.2e-11** is attributed to [R3735] [RB86], which was not opened. [R0618] [Bu97] (read here) gives B = 1.7 +- 0.4e-11 per solar luminosity per year, and the stellar-tables check found the same from [Bu89]. V gives 1.7; N gives 1.7 (the read paper is newer); M would keep 2.2 inside a wider band. Recommendation under the rules as Boss gave them: 1.7.
3. **"Yu and Ho 2018" for the bulge and pitch-angle result (GX2)** is probably Yu and Ho 2019 [YH19] (On the connection between spiral arm pitch angle and galaxy properties, ApJ 871, 194). [YH18] (ApJ 869, 29) is about the wavelength dependence of the pitch angle.
4. **"[R1836]"** can be either of two papers (MNRAS 468, 1850 and MNRAS 472, 2263 [Har17b]); the notes should name the one that holds the pitch and bulge result.
5. **The fO2 buffer.** [OE02] concerns activity coefficients of FeO, NiO, CoO and Mo oxides in silicate melts (by its title). The iron-wuestite calibration is [Hir21]; [Fro91] is the introduction to oxygen fugacity.
6. **Badhwar and O'Neill** exist as two papers, 1994 [BO94] (Long-term modulation of galactic cosmic radiation, Adv. Space Res. 14, 749) and 1996 [BO96] (Galactic cosmic radiation model and its applications, Adv. Space Res. 17, 7). The notes should name which one the code follows.
7. **Plotnykov.** The slope 0.94 (+1.02/-1.07) is [PV26] (the 2025 preprint [R3546]), not [PV20].

**Source against source** (V, N and M applied; the six "widen, keep both" answers stand and are not repeated):

| # | Conflict | V | N | M | Where it stands |
|---|---|---|---|---|---|
| 8 | Halo-fraction scatter: Harmsen 0.33 dex (six discs, 2017, read [Har17]) against Dragonfly 1.01 dex (eight galaxies, [Mer16], second-hand through the Knapen book) | 0.33 (verified over partly verified) | 0.33 if the Dragonfly paper is read and verified, since [Har17] is newer | 1.01 (widest) | M is now the widest range (Boss, 1:12 PM EDT), so 1.01 applies where both sources are in doubt; V and N give 0.33 |
| 9 | Thick-to-thin length ratio against Vc: [YD06] falls, [Com18] median 1.46 rises | no decision (both verified) | [Com18] | wider spread | Not decided; not one of the six |
| 10 | Giant heavy-element mass: [Tho16] against [Cha25] (1.4 to 1.6 times lower, 0.3 to 3 Jupiter masses) | no decision | [Cha25] | the union | Not decided |
| 11 | Moon ice fraction: [Sas10] against [Cil18] | no decision | [Cil18] | the 0.5 to 0.85 bracket (today's) | Not decided |
| 12 | Planet-star Fe/Mg slope: [PV20] against [PV26] | no decision | [PV26] | widest | Not decided |
| 13 | Dwarf stellar-to-halo scatter: [Fit17] 0.5 dex against [Nad20] under 0.19 dex | no decision | [Nad20] | uniform 0.19 to 0.5 (today's) | Already the M outcome |
| 14 | Water-rich share: [R3801] (at most 3 percent) against [R0715] (59 percent) (U) | open until read | 2026 paper | 0 to 0.87 (today's) | Both papers unverified here |
| 15 | Palme and O'Neill: 2003 edition values (read) against 2014 values quoted second-hand | 2003 | 2014 if opened and verified | wider | Open until the 2014 edition is read |
| 16 | Rossby saturation thresholds of different authors | Wright's, read by the lane; never mixed | | | Keep one author's pair of numbers |
| 17 | Rock fraction of Kuiper Belt objects: rising with size against constant 70 percent with porosity | no decision | dates not recorded | today's draw | Not decided |

**The twelve paper-against-project contradictions** of boss-papers section 3 (ice compaction reference pressure, [R3417] Th and U, mantle parameter set, "IPCC 1990" labels, generic Rayleigh error, Earth optical depth 0.0969 against 0.09721, Tabazadeh as an acid source, comet NH3 bound, constant 0.54 dispersion ratio, Wolff tables, Gaillard FMQ range, Mars column) are all decided by V: the paper wins and both values stay in the record. N and M do not apply because one side is a project value or a second-hand value, not a verified source. The 16 Gemini claims contradicted by a read source are also V: the verified figure wins.

## 6. Free sources and tools: what answered

Checked from this container at about 12:50 PM EDT on 2026-10-17 by requesting each site's page (HTTP 200 counts as reachable). Reachable does not mean a bulk download works.

**Reachable:** NIST Chemistry WebBook, CCCBDB, XCOM, ESTAR, PSTAR and chemical kinetics; PubChem; HITRAN, HITEMP, ExoMol; MIST, BaSTI, BHAC15 and ATMO 2020 model grids; CompOSE; Illustris and TNG data; MELTS; RELAB; EXPACS, CREME and OLTARIS; TUV; the Montreal White Dwarf Database; VizieR; NASA Exoplanet Archive; GEOROC and EarthChem; CoolProp; IAPWS; DACE; Colossus and HMFcalc; petitRADTRANS and virga; NuDat; the Crystallography Open Database; IGRF; MESA; Gaia archive; Hypatia; the IUPAC kinetics datasheets; the MPI-Mainz UV/VIS atlas; the Leiden photodissociation database; Crossref, arXiv, Open Library.

**Not reachable from this container (Boss's browser may work):** JPL kinetics, E-AIM, KIDA, UMIST, PARSEC, Kaye and Laby tables, NASA CEA, several GitHub repositories, the USGS spectral library page, Zenodo, the Meteoritical Bulletin Database, the [R2726] free copy (captcha), the Tabazadeh Wiley PDF. The IAEA nuclear-data API returned a Cloudflare challenge and the UNSCEAR 2000 Annex B address returned 404, so half-lives and UNSCEAR coefficients stay unverified.

## 7. References

Every reference cited above, in order of first citation. Mark: V-rec = record confirmed, text not opened; V-read = opened and read by this thread; PT = marked READ by a prior project thread; U = unverified (section 0).

- **[KT10]** [R2265] (Methane absorption coefficients for the jovian planets). *Mark: V-rec; text paid.*
- **[Kar94]** [R2258] (Spectrophotometry of the jovian planets and Titan at 300 to 1000 nm). *Mark: V-rec; text paid.*
- **[Har20]** [R1812] (HITEMP methane line list, 0-13,400 cm-1, i.e. longer than 746 nm, by its abstract). *Mark: V-rec; abstract read.*
- **[Lod03]** [R2726] (Solar system abundances and condensation temperatures of the elements). *Mark: V-rec; text not opened by this thread (free copy behind a captcha).*
- **[IAEA14]** IAEA 2014, Radiation Protection and Safety of Radiation Sources: International Basic Safety Standards, General Safety Requirements Part 3 (GSR Part 3), Vienna; glossary entries 'radiation weighting factor' and 'tissue weighting factor', and the occupational limit paragraphs. *Mark: V-read (the PDF opened; glossary values copied).*
- **[Bog78]** [R0443], J. Chem. Soc. Faraday Trans. 1, 74, 3008 (Rayleigh scattering depolarization ratio and molecular polarizability anisotropy for gases). *Mark: V-rec; text not opened (the files say it holds SO2 and H2S).*
- **[KWW12]** [R2364], Stellar Structure and Evolution, 2nd ed., Springer, [R2364] *Mark: V-rec.*
- **[Poi00]** [R3556], Introduction to the Physics of the Earth's Interior, 2nd ed., Cambridge University Press. *Mark: V-rec.*
- **[SD08]** [R4187], Physics of the Earth, 4th ed., Cambridge University Press, [R4187] *Mark: V-rec.*
- **[LP16]** [R2574], Phys. Rep. 621, 127 (The equation of state of hot, dense matter and neutron stars). *Mark: V-rec.*
- **[Cot14]** Cottaar, Heister, Geochem. Geophys. Geosyst. 15, 1164 (BurnMan: a lower mantle mineral physics toolkit). *Mark: V-rec.*
- **[BT08]** [R0398, R0399], Galactic Dynamics, 2nd ed., Princeton University Press, [R0398] *Mark: V-rec.*
- **[Ken98]** [R2321] (The global Schmidt law in star-forming galaxies). *Mark: V-rec; PT (index 1.4 +- 0.15).*
- **[KE12]** [R2320] (Star formation in the Milky Way and nearby galaxies). *Mark: V-rec.*
- **[Big08]** [R0392] (The star formation law in nearby galaxies on sub-kpc scales). *Mark: V-rec.*
- **[Pag09]** [R3405], Nucleosynthesis and Chemical Evolution of Galaxies, 2nd ed., Cambridge University Press, [R3405] *Mark: V-rec.*
- **[Nel19]** [R3258], The IllustrisTNG Simulations: Public Data Release. *Mark: V-rec.*
- **[ST99]** [R4060] (Large scale bias and the peak background split). *Mark: V-rec.*
- **[Tin08]** [R4367] (Toward a halo mass function for precision cosmology). *Mark: V-rec.*
- **[Gio08]** Giocoli, Tormen and van [R1574, R4492] (The population of dark matter subhaloes: mass functions and average mass loss rates). *Mark: V-rec; PT (Eq. 2 read).*
- **[Mos13]** [R3181]. *Mark: V-rec; PT.*
- **[MBW10]** Mo, van [R3106], Galaxy Formation and Evolution, Cambridge University Press, [R3106] *Mark: V-rec.*
- **[Die18]** [R1111] (COLOSSUS: a python toolkit for cosmology, large-scale structure and dark matter halos). *Mark: V-rec.*
- **[Mur13]** [R3221], Astronomy and Computing 3-4, 23, [R3221] (HMFcalc). *Mark: V-rec.*
- **[Bu97]** [R0618], astro-ph/9711072 (Evolutionary population synthesis models of primeval galaxies: a critical appraisal); the project files list it as '[R3734]'. *Mark: V-rec; V-read (B = 1.7 +- 0.4 e-11 per solar luminosity per year).*
- **[SC05]** [R3877], Evolution of Stars and Stellar Populations, Wiley. *Mark: V-rec (Open Library).*
- **[Bu89]** [R0617] (Evolutionary population synthesis in stellar systems I). *Mark: V-rec; PT.*
- **[Cho16]** [R0773] (MIST I). *Mark: V-rec.*
- **[You81]** [R4835], Appl. Opt. 20, 533 (Rayleigh scattering). *Mark: V-rec; PT.*
- **[Bat84]** [R0282], Planet. Space Sci. 32, 785 (Rayleigh scattering by air). *Mark: V-rec; PT.*
- **[Buc95]** [R0575], Appl. Opt. 34, 2765 (Rayleigh-scattering calculations for the terrestrial atmosphere). *Mark: V-rec; PT.*
- **[GY89]** [R1633], Atmospheric Radiation: Theoretical Basis, 2nd ed., Oxford University Press. *Mark: V-rec (Open Library).*
- **[LO91]** [R2529] (correlated k distribution method). *Mark: V-rec.*
- **[Sea10]** [R3990], Exoplanet Atmospheres: Physical Processes, Princeton University Press, [R3990] *Mark: V-rec.*
- **[Gor22]** [R1634], JQSRT 277, 107949 (The HITRAN2020 molecular spectroscopic database). *Mark: V-rec.*
- **[Ten20]** [R4323], JQSRT 255, 107228 (The 2020 release of the ExoMol database). *Mark: V-rec.*
- **[Pie10]** [R3519], Principles of Planetary Climate, Cambridge University Press, [R3519] *Mark: V-rec.*
- **[LF02]** [R2725] (Atmospheric chemistry in giant planets, brown dwarfs and low-mass dwarf stars I). *Mark: V-rec.*
- **[Bel14]** [R0329], Ind. Eng. Chem. Res. 53, 2498 (CoolProp). *Mark: V-rec.*
- **[Hir21]** [R1951], Geochim. Cosmochim. Acta 313, 74, [R1951] (Iron-wuestite revisited). *Mark: V-rec.*
- **[Fro91]** [R1464], Rev. Mineral. 25, ch. 1 (Introduction to oxygen fugacity and its petrologic importance). *Mark: V-rec.*
- **[CA79]** [R0913], Moon and Planets 21, 155 (Vaporization of comet nuclei: light curves and life times). *Mark: V-rec.*
- **[FS09]** [R1445], Planet. Space Sci. 57, 2053 (Sublimation of ices of astrophysical interest: a bibliographic review). *Mark: V-rec; text paid.*
- **[SS05]** [R4122], Adv. Space Res. 36, 2012 (A review of geomagnetic cutoff rigidities for earth-orbiting spacecraft). *Mark: V-rec.*
- **[Att86]** [R0191], Introduction to Radiological Physics and Radiation Dosimetry, Wiley. *Mark: V-rec.*
- **[Tur07]** [R4455], Atoms, Radiation, and Radiation Protection, 3rd ed., Wiley. *Mark: V-rec.*
- **[Sat15]** [R3924], PLoS ONE 10, e0144679 (PARMA: analytical model for estimating terrestrial cosmic ray fluxes). *Mark: V-rec.*
- **[BS05]** [R0520], Aeronomy of the Middle Atmosphere, 3rd ed., Springer. *Mark: V-rec.*
- **[Troe79]** [R4429], J. Phys. Chem. 83, 114 (Predictive possibilities of unimolecular rate theory). *Mark: V-rec.*
- **[Atk04]** [R0182], Atmos. Chem. Phys. 4, 1461 (IUPAC evaluated kinetic and photochemical data, Vol. I). *Mark: V-rec.*
- **[KR13]** [R2303], Earth Syst. Sci. Data 5, 365 (The MPI-Mainz UV/VIS spectral atlas). *Mark: V-rec.*
- **[Hea17]** [R1873] (Photodissociation and photoionisation of atoms and molecules of astrophysical interest). *Mark: V-rec.*
- **[Har17]** [R0326, R1816] (Diverse stellar haloes in nearby Milky Way mass disc galaxies). *Mark: V-rec; PT.*
- **[BLF20]** [R0338] (Figs 10 and 11). *Mark: V-rec; PT.*
- **[RG16]** [R3790] (Illustris stellar mass assembly). *Mark: V-rec; PT (abstract and Fig 3).*
- **[Nel15]** [R3257], The Illustris Simulation: Public Data Release. *Mark: V-rec.*
- **[Rem22]** [R3730] (Accreted or not accreted? The fraction of accreted mass in galaxies from simulations and observations). *Mark: V-rec.*
- **[DS14]** [R1169] (Parametrizing the stellar haloes of galaxies). *Mark: V-rec; PT.*
- **[YD06]** [R4824]. *Mark: V-rec; PT.*
- **[Com18]** [R0865]; VizieR J/A+A/610/A5. *Mark: V-rec; PT (catalogue fit computed).*
- **[Com21]** [R0860] (age of thick discs against stellar mass; Eq. 1 hot-component mass). *Mark: V-rec; PT.*
- **[SH20]** [R4043] (velocity dispersion of stars in the Milky Way; Table 2). *Mark: V-rec; PT (Table 2 read).*
- **[Mar13]** [R2917] (DiskMass Survey VI). *Mark: V-rec; PT.*
- **[AB09]** [R0202] (Table 6). *Mark: V-rec; PT.*
- **[Hay13]** [R1863]. *Mark: V-rec; PT.*
- **[Ben14]** [R0349] (title: Exploring the Milky Way stellar disk, 714 F and G dwarf stars); journal volume from memory. *Mark: V-rec (arXiv title); volume U.*
- **[Pin19a]** [R3531, R3532] (Fornax 3D: thick disk origin in FCC 170). *Mark: V-rec; PT.*
- **[Pin19b]** [R3531, R3532] (Fornax 3D: thick disks in a cluster environment). *Mark: V-rec; PT.*
- **[Sat24]** [R3530, R3928] (Relatively young thick disks in star-forming late-type galaxies). *Mark: V-rec; PT.*
- **[Cor15]** [R0894] (accretion history of dark matter halos I). *Mark: V-rec; PT.*
- **[FP12]** [R1410] (The baryonic halos of elliptical galaxies: radial distribution of globular clusters and diffuse hot gas). *Mark: V-rec.*
- **[HR14]** [R1811]. *Mark: V-rec.*
- **[Kar16]** [R2268] (SLUGGS: Leo II group). *Mark: V-rec; PT.*
- **[HR18]** [R2038] (sizes of globular cluster systems and their host dark matter haloes). *Mark: V-rec; PT.*
- **[For18]** [R1553], Proc. R. Soc. A 474, 20170616, [R1408] (Globular cluster formation and evolution in the context of cosmological galaxy assembly: open questions). *Mark: V-rec.*
- **[Geo10]** [R1539] (Globular cluster systems in nearby dwarf galaxies III). *Mark: V-rec.*
- **[Beh13]** [R0320]. *Mark: V-rec; PT.*
- **[Pfe18]** [R3508] (E-MOSAICS). *Mark: V-rec.*
- **[Kru14]** [R2505], Class. Quantum Grav. 31, 244006, [R4977] (Globular cluster formation in the context of galaxy formation and evolution). *Mark: V-rec.*
- **[Neu20]** [R3273] (Nuclear star clusters). *Mark: V-rec; PT.*
- **[TOS75]** [R4420] (The formation of the nuclei of galaxies I: M31). *Mark: V-rec.*
- **[Geo16]** [R2631], MNRAS (Masses and scaling relations for nuclear star clusters, and their co-existence with central black holes); volume and page not confirmed. *Mark: V-rec (title); volume U.*
- **[CC10]** [R0876] (Galactic stellar haloes in the CDM model). *Mark: V-rec.*
- **[GK14]** [R1523] (ELVIS; the low-mass slope 1.92). *Mark: V-rec; PT.*
- **[Nad20]** [R3232] (Milky Way satellite census II). *Mark: V-rec; PT.*
- **[Fit17]** [R1387] (FIRE in the field: simulating the threshold of galaxy formation). *Mark: V-rec.*
- **[GK17]** [R1524] (Not so lumpy after all). *Mark: V-rec.*
- **[Har17b]** [R1836] (Galaxy Zoo and SpArcFiRe: constraints on spiral arm formation mechanisms from spiral arm number); which [R1836] paper the project cites is not stated in the notes. *Mark: V-rec; which paper U.*
- **[YH19]** Yu and Ho 2019, ApJ 871, 194, [R4849] (On the connection between spiral arm pitch angle and galaxy properties). *Mark: V-rec.*
- **[YH18]** Yu and Ho 2018, ApJ 869, 29 (Dependence of the spiral arms pitch angle on wavelength as a test of the density wave theory). *Mark: V-rec.*
- **[DB18]** [R1168] (The masses and metallicities of stellar haloes reflect galactic merger histories); the journal reference is not confirmed. *Mark: V-rec (arXiv); journal U.*
- **[Kir13]** [R2368] (Universal stellar mass-stellar metallicity relation for dwarf galaxies). *Mark: V-rec.*
- **[Tre04]** [R4423] (Origin of the mass-metallicity relation, 53,000 SDSS galaxies). *Mark: V-rec.*
- **[Abb18]** [R0013], PRL 121, 161101 (GW170817: measurements of neutron star radii and equation of state). *Mark: V-rec.*
- **[Mil19]** [R3072] (PSR J0030+0451 mass and radius from NICER data). *Mark: V-rec.*
- **[Ril19]** [R3766] (A NICER view of PSR J0030+0451). *Mark: V-rec.*
- **[Mil21]** [R3075] (The radius of PSR J0740+6620 from NICER and XMM-Newton data). *Mark: V-rec.*
- **[Ril21]** [R3767] (A NICER view of the massive pulsar PSR J0740+6620). *Mark: V-rec.*
- **[Rea09]** [R3691], PRD 79, 124032 (Constraints on a phenomenologically parameterized neutron-star equation of state). *Mark: V-rec.*
- **[HPT00]** [R2055] (Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity). *Mark: V-rec.*
- **[Dot16]** [R1150] (MIST 0). *Mark: V-rec.*
- **[Sta83]** [R4189] (The birthline for low-mass stars). *Mark: V-rec.*
- **[Bar15]** [R0254] (BHAC15; tables at CDS). *Mark: V-rec; PT (paper has no tables).*
- **[Sie00]** Siess, A&A 358, 593 (pre-main-sequence tracks). *Mark: V-rec.*
- **[Wri11]** [R4780] (The stellar activity-rotation relationship and the evolution of stellar dynamos). *Mark: V-rec; PT.*
- **[Sku72]** [R4114]. *Mark: V-rec.*
- **[Joh21]** [R2183] (The active lives of stars). *Mark: V-rec; PT.*
- **[Jac12]** Jackson, MNRAS 422, 2024, [R2118] (coronal X-ray - age relation). *Mark: V-rec; PT.*
- **[Med22]** [R3017] (Galactic kinematics and observed flare rates of mid-to-late M dwarfs). *Mark: V-rec; PT.*
- **[Kir21]** [R2374] (Field substellar mass function from the 20-pc census); Table 13. *Mark: V-rec; PT.*
- **[BL93]** [R0607] (The science of brown dwarfs). *Mark: V-rec.*
- **[Phi20]** [R3512] (ATMO 2020 models). *Mark: V-rec.*
- **[Bar03]** [R0253] (evolutionary models for cool brown dwarfs and extrasolar giant planets). *Mark: V-rec.*
- **[Ser17]** [R4024] (The brightness of the red giant branch tip). *Mark: V-rec; PT.*
- **[Ree25]** [R3702] (A stellar evolutionary grid for binary population synthesis: from the main sequence to helium ignition); the core-mass-luminosity fit used in GEN.502. *Mark: V-rec; PT.*
- **[Sal04]** [R3879] (The initial helium abundance of the Galactic globular cluster system). *Mark: V-rec; PT.*
- **[LDZ94]** [R2619] (horizontal-branch stars in globular clusters II: the second parameter). *Mark: V-rec.*
- **[Dot10]** [R1151] (ACS Survey of Galactic Globular Clusters IX: horizontal branch morphology). *Mark: V-rec.*
- **[RB86]** [R3735], in Spectral Evolution of Galaxies, Reidel, p. 195 (Global properties of stellar populations and the spectral evolution of galaxies). *Mark: V-rec (chapter); text not opened.*
- **[BRL97]** [R0355] (chemical evolution of cool white dwarfs). *Mark: V-rec.*
- **[GBB19]** [R1533] (DA and DB white dwarfs from SDSS and Gaia). *Mark: V-rec.*
- **[Kep19]** [R2326] (White dwarf and subdwarf stars in SDSS DR14). *Mark: V-rec.*
- **[WSH19]** [R4755], Am. Mineral. 104, 844 (The condensation temperatures of the elements: a reappraisal). *Mark: V-rec.*
- **[Asp21]** [R0168] (The chemical make-up of the Sun: a 2020 vision); Table 2. *Mark: V-rec; PT.*
- **[Bre16]** [R0535, R4978] (Spectral properties of cool stars: extended abundance analysis of 1,617 planet-search stars); VizieR J/ApJS/225/32. *Mark: V-rec; PT (table downloaded and fitted).*
- **[SLB05]** [R4231], GJI 162, 610 (Thermodynamics of mantle minerals I). *Mark: V-rec.*
- **[SLB11]** [R4232], GJI 184, 1180 (Thermodynamics of mantle minerals II; Table A1). *Mark: V-rec; PT (Table A1 read).*
- **[PO14]** [R3417], in Treatise on Geochemistry 2nd ed. (Cosmochemical estimates of mantle composition); the 2003 first edition is the one read. *Mark: V-rec; 2014 values not opened.*
- **[MS95]** McDonough and Sun 1995, Chem. Geol. 120, 223 (The composition of the Earth). *Mark: V-rec.*
- **[RG14]** [R3842], in Treatise on Geochemistry 2nd ed. (Composition of the continental crust); Table 11. *Mark: V-rec; PT.*
- **[SCvH95]** [R3932] (An equation of state for low-mass stars and giant planets). *Mark: V-rec.*
- **[Tho16]** [R4352] (The mass-metallicity relation for giant planets). *Mark: V-rec; V-read (M_Z = 57.9 +- 7.03 M^(0.61 +- 0.08), 47 planets, scatter factor 1.82 +- 0.09).*
- **[For07]** [R1429] (Planetary radii across five orders of magnitude). *Mark: V-rec.*
- **[Cha25]** [R0710] (Revising the giant planet mass-metallicity relation). *Mark: V-rec.*
- **[Sas10]** [R3922] (Origin of the different architectures of the Jovian and Saturnian satellite systems). *Mark: V-rec; PT.*
- **[CW06]** [R0642] (A common mass scaling for satellite systems of gaseous planets). *Mark: V-rec.*
- **[CW02]** [R0641] (Formation of the Galilean satellites). *Mark: V-rec.*
- **[Cil18]** Cilibrasi, Szulagyi, MNRAS 480, 4355, [R0804] (Satellites form fast and late). *Mark: V-rec; PT.*
- **[BN19]** [R0390] (Using the density of Kuiper Belt objects to constrain their composition and formation history); Table 1. *Mark: V-rec; PT (Table 1 read); free copy PMC7058130 (U).*
- **[PV26]** [R3546]; preprint [R3546] (Evidence of 1:1 slope between rocky super-Earths and their host stars). *Mark: V-rec; PT (slope 0.94, band +1.02/-1.07).*
- **[PV20]** [R3545] (Chemical fingerprints of formation in rocky super-Earths' data). *Mark: V-rec.*
- **[Adi21]** [R0032] (A compositional link between rocky exoplanets and their host stars). *Mark: V-rec.*
- **[Fre09]** [R1457], Phys. Rev. B 79, 054107 (water at ultrahigh pressures). *Mark: V-rec.*
- **[Jou20]** [R2204], JGR Planets 125, e2019JE006176 (Gibbs representation of ice thermodynamics). *Mark: V-rec.*
- **[Hal20]** [R1754] (AQUA: a collection of H2O equations of state). *Mark: V-rec.*
- **[Hap12]** [R1806], Theory of Reflectance and Emittance Spectroscopy, 2nd ed., Cambridge University Press. *Mark: V-rec.*
- **[Kok17]** [R2429], USGS Data Series 1035 (USGS Spectral Library version 7), [R2429] *Mark: V-rec.*
- **[Ste99]** [R4219] (Life-sustaining planets in interstellar space?). *Mark: V-rec.*
- **[SW96]** [R4166] (equation of state for carbon dioxide). *Mark: V-rec.*
- **[ICRP07]** ICRP 2007, Publication 103, Ann. ICRP 37(2-4), Elsevier (the IAEA glossary cites it as the source of its wR and wT). *Mark: U (record not queried; not opened).*
- **[DC13]** [R1065] (taxonomic distribution of asteroids); Table 6. *Mark: V-rec; PT.*
- **[Wal11]** [R4603] (A low mass for Mars from Jupiter's early gas-driven migration). *Mark: V-rec.*
- **[RI17]** [R3681] (Origin of water in the inner Solar System). *Mark: V-rec.*
- **[Kro14]** [R2495], in Treatise on Geochemistry 2nd ed. (Classification of meteorites and their genetic relationships). *Mark: V-rec.*
- **[LB86]** [R0277], J. Petrology 27, 745 (A chemical classification of volcanic rocks based on the total alkali-silica diagram). *Mark: V-rec.*
- **[PA09]** [R3514], Principles of Igneous and Metamorphic Petrology, 2nd ed., Cambridge University Press. *Mark: V-rec.*
- **[Lic16]** [R3434], Isotopic enrichment of forming planetary systems from supernova pollution (journal reference not confirmed). *Mark: V-rec (arXiv title); journal U.*
- **[Die06]** [R1108, R1109] (Radioactive 26Al from massive stars in the Galaxy). *Mark: V-rec.*
- **[Rug09]** [R3844], PRL 103, 072502 (New measurement of the 60Fe half-life). *Mark: V-rec; number not opened.*
- **[Fis15]** [R1380], Geochim. Cosmochim. Acta 167, 177 (High pressure metal-silicate partitioning of Ni, Co, V, Cr, Si and O). *Mark: V-rec.*
- **[Hak18]** [R1752] (A new ab initio equation of state of hcp-Fe). *Mark: V-rec.*
- **[DA81]** [R1204], Phys. Earth Planet. Inter. 25, 297 (Preliminary reference Earth model). *Mark: V-rec.*
- **[Yur14]** [R4855] (ExoMol line lists IV: methane up to 1500 K). *Mark: V-rec.*
- **[HT74]** [R1799] (Light scattering in planetary atmospheres). *Mark: V-rec.*
- **[Wol09]** [R4743], E00D04 (Wavelength dependence of dust aerosol single scattering albedo on Mars). *Mark: V-rec; PT (figures only).*
- **[Khare84]** [R2335] (Optical constants of organic tholins). *Mark: V-rec.*
- **[Tom08]** [R4394], Planet. Space Sci. 56, 669 (A model of Titan's aerosols). *Mark: V-rec; text paid.*
- **[Myh98]** [R3229] (New estimates of radiative forcing due to well mixed greenhouse gases). *Mark: V-rec; PT (Table 3 read).*
- **[PG14]** [R3439] (A non-grey analytical model for irradiated atmospheres). *Mark: V-rec.*
- **[Gui10]** [R1715] (On the radiative equilibrium of irradiated planetary atmospheres). *Mark: V-rec.*
- **[Fre14]** [R1450] (Gaseous mean opacities for giant planet and ultracool dwarf atmospheres). *Mark: V-rec.*
- **[Kin23]** [R2360], Proc. R. Soc. A 104, 333. *Mark: V-rec; PT.*
- **[SU05]** [R4138], JQSRT 92, 293 (Direct measurement of the Rayleigh scattering cross section in various gases). *Mark: V-rec; PT (Table 2 read).*
- **[Has14]** [R1843], 1244797 (Mars' surface radiation environment measured with Curiosity's RAD). *Mark: V-rec; PT (0.21 mGy/day, 0.64 mSv/day).*
- **[BO96]** [R0225], Adv. Space Res. 17, 7, [R0225] (Galactic cosmic radiation model and its applications). *Mark: V-rec.*
- **[GA68]** [R1593] (Solar modulation of galactic cosmic rays). *Mark: V-rec.*
- **[Xap00]** [R4789], IEEE Trans. Nucl. Sci. 47, 486 (Probability model for cumulative solar proton event fluences). *Mark: V-rec.*
- **[SF09]** [R3943] (Chemistry of silicate atmospheres of evaporating super-Earths). *Mark: V-rec.*
- **[Kit16]** [R2380] (Atmosphere-interior exchange on hot rocky exoplanets). *Mark: V-rec.*
- **[Zel91]** [R4877] (Thermodynamic properties of the aqueous sulfuric acid system to 350 K). *Mark: V-rec; PT.*
- **[Tab97]** [R4289] (A new parameterization of H2SO4/H2O aerosol composition). *Mark: V-rec; PT (Eq. 1, Table 1 read).*
- **[Aye80]** [R0210] (On the vapor pressure of sulfuric acid). *Mark: V-rec; text paid.*
- **[KL90]** [R2516] (Binary nucleation of water-sulfuric acid system). *Mark: V-rec.*
- **[Gia60]** [R1550], J. Am. Chem. Soc. 82, 62 (aqueous sulfuric acid solutions and hydrates). *Mark: V-rec.*
- **[WC02]** [R4680] (E-AIM). *Mark: V-rec.*
- **[HM15]** [R2040], Planet. Space Sci. 106, 11 (Photoionization and photodissociation rates in solar and blackbody radiation fields). *Mark: V-rec.*
- **[JPL19]** Burkholder et al., JPL Publication 19-5 (Chemical kinetics and photochemical data for use in atmospheric studies, Evaluation No. 19); the project notes disagree on 2019 against 2020. *Mark: U.*
- **[Mil24]** [R3071] (The UMIST Database for Astrochemistry 2022). *Mark: V-rec.*
- **[Wak12]** [R4591] (KIDA). *Mark: V-rec.*
- **[Doy20]** [R1152] (title as returned by arXiv: Where are the extrasolar Mercuries?); Table 1 is the redox source in the tunable-gaps note. *Mark: V-rec; PT.*
- **[Doy23]** [R1153] (New chondritic bodies identified in eight oxygen-bearing white dwarfs); Table 4. *Mark: V-rec; PT.*
- **[Sos20]** [R4154], Sci. Adv. 6, eabd1387 (Redox state of Earth's magma ocean and its Venus-like early atmosphere). *Mark: V-rec; PT.*
- **[Gai22]** [R1488], Earth Planet. Sci. Lett. 577, 117255 (Redox controls during magma ocean degassing). *Mark: V-rec; PT.*
- **[GS14]** [R1491], Earth Planet. Sci. Lett. 403, 307 (volcanic degassing chemistry in a comparative planetology perspective). *Mark: V-rec.*
- **[Hir12]** [R1953], Earth Planet. Sci. Lett. 341-344, 48 (Magma ocean influence on early atmosphere mass and composition). *Mark: V-rec.*
- **[Ems21]** [R1268] (The New Generation Planetary Population Synthesis, NGPPS). *Mark: V-rec.*
- **[LP22]** [R2774] (Density, not radius, separates rocky and water-rich small planets orbiting M dwarf stars). *Mark: V-rec.*
- **[Hor18]** [R2009], Nat. Astron. 2, 303 (Haze production rates in super-Earth and mini-Neptune atmosphere experiments). *Mark: V-rec.*
- **[Tra06]** [R4415], PNAS 103, 18035 (Organic haze on Titan and the early Earth). *Mark: V-rec.*
- **[AM01]** [R0024] (Precipitating condensation clouds in substellar atmospheres). *Mark: V-rec.*
- **[Roo22]** [R3807] (virga sedimentation model). *Mark: V-rec.*
- **[ZC17]** [R4862] (The cosmic shoreline). *Mark: V-rec.*
- **[BT25]** [R0364] (The 3D cosmic shoreline). *Mark: V-rec; PT.*
- **[MG25]** [R3037] (An empirical determination of the cosmic shoreline). *Mark: V-rec; PT.*
- **[CK17]** [R0692], Atmospheric Evolution on Inhabited and Lifeless Worlds, Cambridge University Press. *Mark: V-rec.*
- **[DG14]** [R1137] (Abiotic ozone and oxygen in atmospheres similar to prebiotic Earth). *Mark: V-rec; PT.*
- **[Ben19]** [R0340] (Water vapor and clouds on the habitable-zone sub-Neptune exoplanet K2-18b). *Mark: V-rec.*
- **[HP22]** Harrington [R3535] (A survey of CO, CO2 and H2O in comets and Centaurs); Table 2. *Mark: V-rec; PT.*
- **[Gao20]** [R2612], Nat. Astron. 4, 951, [R1513] (Aerosol composition of hot giant exoplanets dominated by silicates and hydrocarbon hazes). *Mark: V-rec.*
- **[MC11]** [R3206] (The chemical composition of comets); Table 4. *Mark: V-rec; PT.*
- **[BN85]** [R0266] (Trapping and release of gases by water ice and implications for icy bodies). *Mark: V-rec.*
- **[BN87]** [R0265], Phys. Rev. B 35, 2427 (Amorphous water ice and its ability to trap gases). *Mark: V-rec.*
- **[BN88]** [R0267], Phys. Rev. B 38, 7749 (Trapping of gas mixtures by amorphous water ice). *Mark: V-rec.*
- **[NB96]** [R3316] (Enrichment of CO over N2 by their trapping in amorphous ice). *Mark: V-rec.*
- **[NB05]** [R3317] (A ~25 K temperature of formation for the submicron ice grains which formed comets). *Mark: V-rec.*
- **[Rub15]** [R3836] (Molecular nitrogen in comet 67P). *Mark: V-rec.*
- **[BB17]** [R0435], Phil. Trans. R. Soc. A 375, 20160252 (The composition of cometary ices). *Mark: V-rec.*
- **[CII]** Festou, Keller and Weaver (eds) 2004, Comets II, University of Arizona Press. *Mark: V-rec (Open Library).*
- **[Bie15]** [R0386] (Abundant molecular oxygen in the coma of comet 67P). *Mark: V-rec.*
- **[Bal15]** [R0249], Sci. Adv. 1, e1500377 (Detection of argon in the coma of comet 67P). *Mark: V-rec.*
- **[Rub18]** [R3834], Sci. Adv. 4, eaar6297 (Krypton isotopes and noble gas abundances in the coma of comet 67P). *Mark: V-rec.*
- **[BNO07]** [R0268] (Trapping of N2, CO and Ar in amorphous ice). *Mark: V-rec.*
- **[Rob24]** [R3782], cometary volatile abundance database (CC BY 4.0, as the tunable-gaps lane records); no bibliographic record found by Crossref in this thread. *Mark: U.*
- **[DR16]** Dello [R1063] (Emerging trends and a comet taxonomy based on the volatile chemistry measured in thirty comets). *Mark: V-rec.*
- **[Bar17]** [R0259] (Carbon-rich dust in comet 67P measured by COSIMA). *Mark: V-rec.*
- **[AH12]** [R0048] (Cometary volatiles and the origin of comets). *Mark: V-rec.*
- **[Jor05]** [R2193], astro-ph/0508219 (ACS Virgo Cluster Survey X). *Mark: V-rec; PT (abstract).*
- **[OE02]** [R3365], Chem. Geol. 186, 151 (the activity coefficients of FeO, NiO, CoO, MoO2 and MoO3 in silicate melts, by its title); not an IW-buffer calibration. *Mark: V-rec.*
- **[BO94]** [R0224], Adv. Space Res. 14(10), 749, [R0224] (Long-term modulation of galactic cosmic radiation and its model for space exploration). *Mark: V-rec.*
- **[Mer16]** [R3041] (The Dragonfly Nearby Galaxies Survey I). *Mark: V-rec; PT (second-hand, through the Knapen book).*

Count: 236 references: 3 read in this thread, 230 with the record confirmed, 3 unverified.

## 8. Limits of this report

- Book entries confirm the book exists (a Crossref or Open Library record); the chapter named is by topic and was not opened. Page and chapter locators are unverified.
- V-rec is a bibliographic check, not a reading. Only the IAEA Safety [R0618, R4352] were opened in this thread; every PT mark rests on an earlier thread's reading.
- The King-factor form with polarizability, the polytrope energy factors (3/7 and 3/4), and the 0.7 percent fusion fraction are textbook results I re-derived or recalled (the fusion fraction from the recalled atomic masses of H-1 and He-4; NIST atomic weights are reachable to check). Check each against a text before coding.
- The "quality" words are my judgement from the gap plan's closure checks; they are not new tests.
- I did not compare the project code with current main; the TODO items were read as filed.
- Several ids in earlier notes were not checked (anything not in section 7). Every id in this report was checked; ids I had guessed from memory and found wrong (about 20) were replaced, and the references that could not be confirmed carry U.
