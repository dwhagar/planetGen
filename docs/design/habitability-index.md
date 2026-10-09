# Planetary habitability index

How planetGen will score how habitable each planet and moon is, and what
the planet model needs before it can. This condenses Boss's research
documents in this folder:

- "Planetary Habitability and Speculative Xenobiology.md" (the master doc)
- "Planetary Habitability Index.md" and "Mathematical and Algorithmic
  Implementation of the Planetary Habitability Index.md" (PHI-4)
- "Atmospheric Toxicity.md"
- "Chemical Habitability.md"
- "Naturally Occuring Ionizing Radiation.md"
- "Speculative Xenobiology Examples.md" and "Speculative Xenobiology
  Extremes.md"

Boss (2026-10-03 05:38Z): "Create a habitability index based on the
pressure, temperature, composition, etc..." and "Refactor all planetary
classes to better align with our habitability index system and the
research for the different kinds of life, etc."

The research of 2026-10-09 that fixes the numbers is in two further notes:

- [atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md):
  when a planet keeps an atmosphere, the hot- and cold-zone classes (GEN.91),
  the class windows and life stage against time habitable (GEN.92).
- [activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md):
  stellar XUV and flares, planetary magnetic fields, surface dose, water,
  oceans and ocean chemistry (GEN.86, GEN.87, GEN.88).

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else
is a recommendation. Informs: GEN.83 to GEN.92.

## Decisions already taken

- **Score structure (GEN.84).** Boss (2026-10-07 17:11Z, confirmed 2026-10-09
  07:48Z): PHI-4's four domains and colour tiers for display, with the
  Xenobiology doc's three tiers (`PHI_bio`, `PHI_cpx`, `Phi_tech`) as the
  scores behind them. The recommendations below fill in the glue; they do not
  reopen that choice.
- Each planet also gets a human equipment profile: shirtsleeve, breathing
  mask, mask with scrubber, pressure suit, or full life support with
  radiation hardening (GEN.89).
- GEN.29 stays in phase 2 with the class refactor (same dates).

## 1. What the code has today

- **Per planet**: class and zone, radius, mass, density, gravity, orbit,
  rotation period, surface temperature (grey body times a per-class
  greenhouse factor), surface pressure, atmosphere density and scale
  height. `atmosphere` and `composition` are free text. Pressure scales
  with gravity linearly (`_atmosphere_retention_factor`, "a starting point,
  not derived"); there is no escape, temperature or age term.
- **Missing**: partial pressures of any gas, ocean or land fraction, pH,
  salinity, magnetic field, radiation dose, mantle redox, ozone, humidity.
- **Per star**: type and class, mass, radius, temperature, luminosity,
  age, habitable zone, heliosphere radius. No XUV or flare activity. The
  habitable zone is a fixed pair of fluxes (1.1 and 0.53) for every star, and
  a main-sequence star's luminosity is constant for its life.
- **Spin**: only moons are tidally locked; a planet's rotation is a uniform
  10 to 1,400 h draw. GEN.104 adds the spin vector.
- **Rogue planets** already have internal heat flux, a surface regime,
  ice-shell thickness and ocean depth (`rogue_surface.py`); nothing in a
  star system has these.
- **Nebulae**: a system's `surrounding_cloud` is set only when it is
  loaded from the database, so generation can't see it.
- **Life** is decided by class (`HABITABLE_PLANET_CLASSES`), and the life
  chemical by the star's spectral letter. Classes N and Q carry a life
  chemical and a stored timeline but are not in the habitable list, so the
  page hides the timeline.

## 2. The scores

The documents give two structures; Boss approved the layering (decision above).

- **PHI-4**: four domains, Pressure, Temperature, Chemistry and Radiation,
  each mapped through a logistic transfer function
  `S_k(x) = [1 + exp(-kappa_L (x - x_min))]^-1 [1 + exp(kappa_R (x - x_max))]^-1`
  (the Math doc's section 3; `kappa_L` and `kappa_R` are not given),
  combined as a weighted geometric mean times hard gates, with colour tiers
  (Blue, Green, Yellow, Red).
- **The Xenobiology doc**: three scores per planet:
  - `PHI_bio` (microbial life) from water activity and chaotropicity
    (`L_solv`), C/N/H and phosphorus (`L_chem`), available energy (`L_ener`:
    light, chemistry, radiolysis) and dose (`L_rad`);
  - `PHI_cpx` (complex life) = `PHI_bio` times a metazoan gate
    `T_metazoa = exp(-(pCO/pCO_tox)^2) exp(-(pCO2/pCO2_tox)^2) f(pO2)`;
  - `Phi_tech` (human operability) = `(M_press M_therm M_isru M_rad)^(1/4)`.

### 2.1 Recommended composite

1. **Inputs** from GEN.85 to GEN.88: total pressure, partial pressures,
   surface temperature (and wet-bulb), surface dose, water activity, energy
   flux, nutrients.
2. **Four domain scores for humans**, each 0 to 1 and each the minimum over
   its sub-limits. The tier comes from the worst sub-limit's band (2.3). The
   PHI-4 colour describes how hostile the unprotected environment is; the
   equipment tier (4) says what keeps a person alive in it.
3. **`PHI_bio` = Liebig blend** of the four terms:
   `PHI = L_min^0.65 * (L_solv L_chem L_ener L_rad)^(0.35/4)`. Plain `L_min`
   is the zero-constant alternative. The Math doc's weighted geometric mean has
   the same flaw as the Master doc's (one near-zero factor is forgiven), so
   the blend replaces it wherever a number is shown; the hard gates `Theta_j`
   (runaway, ice-sealed seafloor, lethal dose) stay.
4. **`PHI_cpx`** = `PHI_bio * T_metazoa` with `f(pO2)` a trapezoid: 0 below
   8 kPa, 1 from 16 to 40 kPa, falling to 0 at 160 kPa (the doc's 8 and 40,
   plus 160 for acute CNS toxicity).
5. **`Phi_tech`** keeps the Master doc's four-factor form, but `M_press` comes
   from the equipment tier (4) instead of its own pressure formula (which
   scores the 6.3 to 20 kPa band 1.0, wrongly, see C1). Default tier values
   1.0, 0.9, 0.8, 0.5, 0.15 for tiers 0 to 4 [D], to be calibrated against the
   doc's table (dense CO2 0.78 and stalled CO 0.88 at tier 2; Hycean 0.32 and
   Mars analog 0.62 at tier 3) once the other three factors exist. Drop `M_isru` from the first version (no
   in-situ resource data). Set `lambda` in `M_rad` so a 1,000 g/cm2 column
   halves the unshielded dose.
6. **Display**: the worst of the four human domain tiers gives the colour; the
   equipment tier comes from the same inputs.
7. An ESI-style "Earth-likeness" figure is not recommended: the literature
   (Cockell et al. 2016; phl.upr.edu) warns that it measures resemblance, and a
   second ranking would be misread.

### 2.2 Calibration of the blend

Re-evaluating the Master doc's own worked table (it gives four `L` values and
`PHI_bio` per archetype) [C]:

| World | Stated `PHI_bio` | Geometric mean | Min | Blend, weight 0.65 |
|---|---|---|---|---|
| Earth | 0.98 | 0.985 | 0.980 | 0.98 |
| Dense CO2, 3 bar | 0.82 | 0.934 | 0.880 | 0.90 |
| Hycean 20 bar H2 | 0.74 | 0.855 | 0.720 | 0.76 |
| Stalled CO | 0.79 | 0.914 | 0.860 | 0.88 |
| Ice-sealed ocean | 0.12 | 0.378 | 0.050 | 0.10 |
| Cryo-brine (Mars) | 0.18 | 0.353 | 0.150 | 0.20 |
| Desiccated dune | 0.05 | 0.337 | 0.020 | 0.05 |
| Europan ocean | 0.68 | 0.774 | 0.580 | 0.64 |

RMS error against the stated values: geometric mean 0.169, product 0.157, min
0.057, blend (weight 0.65) 0.046 (weights 0.55 to 0.75 all give 0.046 to
0.050). The Master text calls the mean "non-compensatory"; the geometric mean
is only weakly so. The blend lands within 0.1 of every row, so the table was
plausibly made with a Liebig-like rule. The `L` constants below still need
fixing, and the table must be recomputed once the `L` terms exist.

### 2.3 Display bands (reconciled)

Written numerically from PHI-4 with the corrections of section 3. Units: kPa
for pressure and partial pressure, mSv/yr for human dose, Gy/yr for life.

| Domain | Blue | Green | Yellow | Red |
|---|---|---|---|---|
| Pressure `P_tot` | 50 to 250 kPa | 20 to 50 (O2 mask) or 250 to 400 (diluent mask) | 6.3 to 20 (pressure suit or MCP) or 0.4 to 7 MPa (rebreather, then diving suit above 1 MPa) | under 6.3 kPa, or over 7 MPa (beyond the tested record; theoretical limit 10 to 12 MPa) |
| Temperature | 0 to 31 C dry, wet-bulb under 31 | -20 to 0 or 31 to 45 C, wet-bulb 31 to 35 | 45 to 120 C (liquid-cooled garment), -120 to -20 C (heated suit), wet-bulb over 35 | over 122 C (biology ceiling, sterilizes) or under -120 C |
| Chemistry, air | pO2 16 to 50; pCO2 < 0.5; pCO < 0.005; H2S < 0.001; SO2 < 0.0005 | pO2 8 to 16 or 50 to 160; pCO2 0.5 to 1 | pCO2 1 to 5; pCO 0.005 to 0.1; H2S or SO2 over chronic limit; pO2 under 8 (supplied) | pCO > 0.1; pCO2 > 5; pO2 under 6 or over 160 unprotected |
| Chemistry, water | Neutral pH, a_w > 0.9, no heavy metals or chaotropes | Mildly saline or carbonated | Chaotropic brine or soda ocean (pH 9 to 11.5); a_w toward 0.605 | Acid sulfate (pH 1 to 4.5, pH 0 to 1 the tail) or chaotropicity over 73.8 kJ/kg |
| Radiation | 20 mSv/yr or less (50 if Boss prefers PHI-4's figure) | 20 to 100 mSv/yr | 0.1 to 1 Sv/yr human (shielded work), or 1 to 10 Gy/yr (microbes fine) | over 10 Gy/yr |

PHI-4 left holes that the table fills: 90 to 122 C and below -50 C in
Temperature, 100 mSv/yr to 10 Gy/yr in Radiation, and a Chemistry Yellow
(CO up to 10 percent) that overlaps its Red (CO over 0.1 kPa, C16).

### 2.4 Constants that were undefined or inconsistent (defaults)

All [D] unless marked.

| Constant | Default | Note |
|---|---|---|
| `k_w` (water-activity sigmoid steepness) | 30 per unit a_w | 0.5 at 0.605, 0.95 at 0.705 |
| `chi_crit`, `sigma_chi` | 73.8, 10 kJ/kg | From the doc (Ball and Hallsworth 2015 [R]) |
| `K_P` | 1 micro-mol/L | From the doc |
| `X_ref` (C, N, H inventories) | 0.3 of Earth's | So Earth's `tanh` is 0.997 |
| `Phi_ref` (energy) | 5 W/m2 usable PAR-equivalent | About 5 percent of mean Earth surface PAR, chemosynthetic and radiolytic terms in the same units |
| `D_thresh`, `sigma_D` | 10, 100 Gy/yr | From the doc; Eigen's origin-of-replicators threshold, not an extremophile limit |
| `sigma_P` | 200 kPa | `M_press` 0.37 at 450 kPa, 0.05 at 850 kPa |
| `k_T`, `T_wb,crit` | 0.5 per K, 31 C | From the doc |
| `alpha`, `E_specific` (`M_isru`) | dropped | No in-situ resource data |
| `kappa_L`, `kappa_R` (PHI-4 transfer) | not needed | Domain scores are minimum-based; tiers carry the display |

## 3. Conflicts between the documents, and their resolution

Found by reading all of them (2026-10-09). "Master" is the Xenobiology master
doc, "Atm" Atmospheric Toxicity.md, "Chem" Chemical Habitability.md, "Math" the
Mathematical and Algorithmic doc, "PHI-4" Planetary Habitability Index.md.

| # | Quantity | Values found | Where | Resolution |
|---|---|---|---|---|
| C1 | Breathing-mask pressure band | 6.3 to 250 kPa (`M_press` = 1.0); 50 to 250 kPa; "10 to 20 kPa green, pure-O2 mask" | Master `M_press` (line 445) and conclusions (512); Atm "Breathing-Mask-Only Paradigm" and Sub-Baric Ice Worlds; PHI-4 Pressure | 50 to 250 shirtsleeve; 20 to 50 mask; below 20 a pressure suit (4). The 6.3 to 20 kPa pure-O2 claim fails the alveolar equation |
| C2 | Sub-baric pure-O2 rebreather | "10 to 20 kPa pO2 matches alveolar oxygenation" | Atm Sub-Baric Ice Worlds | At 10 to 20 kPa total, alveolar pO2 of pure O2 is -3 to +7 kPa; a mask needs 19 to 21 kPa or more (4) |
| C3 | CO2 toxicity | 0.93 kPa chronic; 1 to 2; 2.0 tox; 5.0 acute; "0.005 to 0.05 bar" (0.5 to 5 kPa); "0.005 to 0.01 bar" (0.5 to 1 kPa) | Master physiological table and `T_metazoa`; Atm CO2 row and HZCL; Chem thick-atmosphere paragraph; PHI-4 Chemistry | NASA's spacecraft limit is 3 mmHg (0.4 kPa) over one hour [S]; the docs' 0.93 exceeds it. 0.4 shirtsleeve ideal (limit 0.5), 1 scrubber line, 2 metazoan gate scale, 5 acute. Convert every "bar" figure to kPa |
| C4 | CO thresholds | 0.005 kPa chronic (50 ppm); 0.01 (100 ppm, `pCO_tox`); 0.01 to 0.1 "evolutionary barrier"; 0.1 acute; "1 to 10 percent, 0.01 to 0.1 bar" | Master, Atm table, Chem (lethal at 100 to 200 ppm), PHI-4 | 0.005 chronic (mask threshold), 0.01 metazoan gate scale, 0.1 acute. "0.01 to 0.1 bar" is a unit slip: the stalled-CO world has pCO of 1 to 10 kPa, far above lethal unprotected |
| C5 | O2 lower bound | 8.0 kPa (PHI-4 blue floor); 8.0 to 9.6; "< 6.0 hypoxic syncope"; `pO2 >= 10` (Math Category B) | PHI-4; Master; Atm; Math 1.2 | 16 shirtsleeve, 8 to 16 mask, under 8 supplied O2. 8 kPa is a short-exposure level, not a long-term range; La Rinconada (pO2 about 11 kPa) shows acclimatised permanence below 16 [R] |
| C6 | O2 upper bound | 40 (chronic, Master table); 40 to 53; 50 (`f(pO2)`); 53 (PHI-4 blue); 160 acute | Master; Atm; PHI-4 | 50 chronic (Lambertsen: about 0.5 bar tolerated indefinitely [S, forum summary]), 160 acute CNS |
| C7 | Wet-bulb limit | 31 C (Master `M_therm`, PHI-4); 35 C (Math Category B) | Master (line 447); PHI-4; Math 1.2 | 31 comfort limit (tier 0), 35 survivability limit (tier 3) |
| C8 | Maximum pressure for humans | 12 MPa (PHI-4 red); 10 to 12 MPa (Master text); 7.0 MPa (Master table, Hydra X 7.11 MPa) | PHI-4; Master Hyperbaric Regimes | 7 MPa tested record, 10 to 12 MPa theoretical work-of-breathing limit |
| C9 | Hyperbaric N2 | narcosis at 300 to 400 kPa; mask band to 250; PHI-4 green to 400, yellow 500 kPa to 10 MPa | Atm; Master; PHI-4 | Tiers in 4; Ar narcotic at half the N2 partial pressure |
| C10 | Radiation thresholds | `D_thresh` 10 Gy/yr; occupational 50 mSv/yr; "Mars analog `L_rad` = 0.65" | Master (line 437 and archetype table); PHI-4 | Two thresholds (life, human) are fine if stated; the Mars row is a unit slip: Mars is 0.077 Gy/yr, so `L_rad` is about 1. Tiers in 2.3 |
| C11 | Complex-life CO2 limit | HZCL from "0.005 to 0.05 bar" vs 0.93 to 2 kPa vs the dense-CO2 archetype's 2 to 10 bar | Atm | 0.5 to 5 kPa upper range; the dense-CO2 archetype is lethal to metazoans at any of these |
| C12 | `PHI_bio` reproduction | Formula does not reproduce 7 of its 8 table rows | Master table after line 471 | Liebig blend (2.1, 2.2) |
| C13 | Undefined constants | `k_w`, `X_ref`, `Phi_ref`, `sigma_P`, `alpha`, `E_specific`, `lambda`, `kappa_L`, `kappa_R` | Master; Math | Defaults in 2.4 |
| C14 | Desiccated super-Earth | "1 bar, simple masks, evaporative cooling" (tech 0.81) vs Venus-type 9 MPa, 737 K | Master archetype; Atm Post-Runaway Desiccated | Venus type default; 1-bar dune a rare variant (about 10 percent) ([atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md) 4.2) |
| C15 | Hycean vs ice-sealed | `PHI_bio` 0.74 vs 0.12 | Master archetype table | A deep-ocean Hycean has high-pressure ice at its floor; score it with the ice-sealed `L_chem` (same document, 4.3) |
| C16 | CO in PHI-4 Chemistry | Yellow allows CO up to 10 percent; Red starts above 0.1 kPa | PHI-4 Chemistry | Colour is the unprotected hazard: CO over 0.1 kPa is Red, and the equipment tier (mask with Hopcalite filter, tier 2) is what makes the stalled-CO world workable |

Related corrections to the radiation, ocean and atmosphere passages of these
documents are listed in section 6 of
[activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md)
and section 8 of
[atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md).

**The earlier list of problems, answered.** (1) `PHI_bio`'s geometric mean:
replaced (2.1). (2) Undefined constants: fixed (2.4). (3) The `L_rad` and Mars
disagreement: a mSv/Gy unit slip (C10). (4) Threshold differences: C1 to C9.
(5) Only M, K and G hosts covered: the activity defaults in the activity
document span 0.1 to 1.4 Msun and give A, F and B stars no flares; O and B
lifespans are too short for life (life-stage rules in the atmospheres
document, 7); white dwarfs, giants and binaries still have no research.
(6) Nebula effects on planets: section 7.

## 4. Equipment tiers

Derived from partial pressures with numeric thresholds. Reference
implementation: a pure function of dry ambient values (total and partial
pressures, dry-bulb and wet-bulb temperature, dose) [C]. "Alveolar" means
`PAO2 = x_O2 (P - 6.3) - 5.3/0.8` kPa at 37 C with a respiratory quotient of
0.8 (the Master doc's alveolar gas equation with `pH2O` = 6.3 kPa and
`pACO2` = 5.3 kPa).

The alveolar numbers behind the mask band [C]:

| Total P (kPa) | PAO2, pure O2 | PAO2, air | Altitude equivalent |
|---|---|---|---|
| 101.3 | 88.4 | 13.3 | sea level |
| 50 | 37.1 | 2.5 | about 5,700 m |
| 33.7 | 20.8 | -0.9 | Everest summit |
| 21 | 8.1 | -3.6 | |
| 18.8 | 5.9 | -4.0 | 40,000 ft |
| 15 | 2.1 | -4.8 | |
| 11.6 | -1.3 | -5.5 | 50,000 ft |
| 6.3 | -6.6 | -6.6 | Armstrong line |

Pure oxygen needs at least 18.9 kPa total for PAO2 = 6 kPa, 20.9 kPa for 8
kPa, and 26 kPa to match sea-level air (13.3 kPa) [C]. Aviation practice
agrees: 100 percent oxygen holds an "equivalent 10,000 ft" up to about
40,000 ft, positive-pressure breathing above that, a pressure suit above about
50,000 ft, ebullism at 62,000 to 63,000 ft [S: litfl, Wikipedia "Pressure
suit"]. The mask-only band therefore starts near 20 kPa, not 6.3.

| Tier | Name | Conditions (all must hold) |
|---|---|---|
| 0 | Shirtsleeve | 50 <= P_tot <= 250 kPa; 16 <= pO2 <= 50; pCO2 < 0.5; pCO < 0.005; pH2S < 0.001; pSO2 < 0.0005; T_dry -20 to 45 C and wet-bulb < 31 C (ordinary weather clothing); dose <= 50 mSv/yr |
| 1 | Mask | Any of: 20 <= P_tot < 50 (O2-enriched mask; pO2 of air too low); 8 <= pO2 < 16; pO2 50 to 160 (diluent mask); pCO2 0.5 to 1; 250 < P_tot <= 400; T 45 to 90 C, or -20 to -50 C with insulation; wet-bulb 31 to 35 C; dose 50 to 100 mSv/yr |
| 2 | Mask and scrubber | Any of: pCO2 >= 1; pCO >= 0.005 (Hopcalite); H2S or SO2 above the chronic limit; 400 < P_tot <= 1,000 kPa (heliox or trimix rebreather); pO2 < 8 with P_tot >= 20 and a hostile base gas |
| 3 | Pressure suit | P_tot < 20 kPa (6.3 to 20: gas-pressurised or mechanical-counterpressure suit plus helmet; below 6.3: full suit); P_tot > 1 MPa (atmospheric diving suit); pO2 > 160; T 90 to 120 C (liquid-cooled garment) or -50 to -120 C (heated suit); wet-bulb > 35 C; dose 0.1 to 1 Sv/yr (shielded habitat, limited surface time) |
| 4 | Full life support | P_tot > 7 MPa (beyond the tested human record, Hydra X at 7.11 MPa); T < -120 or > 120 C; corrosive acids (H2SO4, HF, Cl); dose > 1 Sv/yr (radiation hardening); no usable atmosphere and no suit-compatible thermal range |

The tier is the highest any condition demands. The dose rows follow the
bands of 2.3: a looser rule (tier 0 below 1 Sv/yr) would put Mars (0.24
Sv/yr) and the Moon (0.52) at shirtsleeve. The cold range of -50 to -120 C,
which PHI-4 leaves unassigned, is tier 3 [D].

Notes on sources and judgements:

- 16 kPa pO2 as the shirtsleeve floor: diving practice's lower O2 limit is about
  0.16 atm [S: secondary summary, weak]. The unaided human floor near 15 to 16
  kPa corresponds to about 2,500 m altitude.
- 50 kPa O2 as the long-exposure ceiling: Lambertsen about 0.5 bar [S, forum
  summary]; acute CNS toxicity at 160 kPa is the standard diving figure [R].
- pCO2: NASA's spacecraft limit is 3 mmHg (0.4 kPa) averaged over one hour,
  down from 3.8 to 7.5 mmHg [S: NASA OCHMO technical brief]. OSHA's 8-hour
  limit of 0.5 percent is 0.5 kPa at 1 bar [R]. 5 kPa is the acute-acidosis
  line.
- pCO 0.005 kPa is 50 ppm at 1 bar (OSHA PEL type value [R]); 0.1 kPa is about
  1,000 ppm, near the NIOSH IDLH of 1,200 ppm [R].
- Atmospheric diving suits run at 1 atm inside to around 300 m (3 MPa) in
  common products [R], so that tier starts at 1 MPa, not the 7 MPa record.
- The 122 C limit for life (PHI-4) is a biology ceiling for the Red band;
  suited humans work to about 120 C [R].

Worked results from the reference function [C]: Earth sea level and Class M
typical (82 kPa, 286 K) give tier 0; Everest summit air (33.7 kPa, -30 C)
gives 1; the stalled-CO world (CO 5 kPa), dense CO2 at 3 bar, and an 8-bar
N2/Ar world with 2.6 percent O2 give 2; Mars, Class L (2 kPa), Class P (7.3
kPa, 205 K), a 20-bar Hycean surface and the Moon (plus radiation) give 3;
Titan (-179 C), Venus and Class N give 4.

## 5. Inputs to build, in order

| # | Piece | Needs |
|---|---|---|
| 0 | Planet spin and tidal locking (GEN.104); atmosphere retention module ([atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md) 3) | mass, radius, insolation, star age |
| 1 | Mantle redox (reduced, intermediate, oxidized) from mass and differentiation | mass, radius |
| 2 | Atmosphere species: partial pressures of O2, CO2, CO, N2, Ar, H2, H2O, CH4, H2S, SO2, with the free text generated from them | 0, 1, pressure, class |
| 3 | Stellar activity: saturation phase, `L_XUV/L_bol`, flare rates by mass and age ([activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md) 2) | star mass, age, type |
| 4 | Magnetic field from mass, core fraction, age, rotation and locking (same document, 3) | mass, spin, 3 |
| 5 | Surface dose from column mass (P0/g), field, cosmic rays (higher inside a compressed heliosphere) and flares, with an ozone-loss flag (same document, 4) | 2, 3, 4, nebula containment |
| 6 | Hydrosphere: water fraction, ocean and land fraction, depth with the high-pressure-ice cap, ice shell; reuse `rogue_surface` in a shared module (same document, 5) | zone, temperature, class |
| 7 | Ocean chemistry: acid sulfate, neutral, soda, chloride brine or ice-sealed, with pH, water activity and phosphorus | 1, 2, 6 |
| 8 | The scores, tier and equipment profile | all of the above |

Values are stored in their own columns, not JSON, so they can be searched.

## 6. Classes and life

Known conflicts between today's classes and the research (checked 2026-10-09):

- **N** (Venus analog, about 737 K and 9.2 MPa) has a life chemical and an
  uncapped stored timeline although it is not in the habitable list; the
  thermal ceiling for life is 122 C (395 K). N also sits in the ecosphere at
  S about 1.05, below the runaway threshold (1.107 for 1 Earth mass); it
  should be a hot-zone, lifeless class. **E** (about 399 K) is above the
  ceiling too.
- **K** (about 605 Pa) and **L** (about 2 kPa) are below the Armstrong limit
  (6.3 kPa); L has land animals at a tiny O2 partial pressure.
- **P** (about 7 to 10 kPa, 204 K) can reach civilization; the research gives
  sub-baric snowballs minimal native life.
- **A** and **B** (13 to 17 kPa of SO2/CO2 or Na/He on 3.4 to 3.7 km/s bodies at
  S about 4) are about 50 times beyond the cosmic shoreline.
- **S** has no atmosphere in any zone, but a 2 to 10 Earth-mass rock in the
  ecosphere or cold zone should have one. **T** and **I** overlap. Rogue S
  reaches 17,600 km; a rocky 10 to 16 Earth-mass rogue should be an R.
- Classes whose atmosphere names O2 imply life made it, whatever their
  life stage.
- No class is an Earth-size world with air that life never processed.
  Planned class Z (lifeless temperate world) is that world; the research's
  "photochemically stalled CO" and "oxidized dense CO2" worlds are its
  likely forms.

Candidate new classes from the research: a dense-CO2 outer-edge world, a
hydrogen-rich (Hycean) super-Earth in the cold zone, a desiccated
post-runaway super-Earth in the hot zone (the "S and V in the other zones"
Boss asked for), an ice-sealed deep waterworld, a soda world and a
hyperbaric N2/Ar world. The recommended windows, letters and PR order are in
[atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md) 5:
the hot and cold counterparts are a re-scoped N and an extended X (plus a
Hycean flag on R), with Class S given air by rule, and no new letters.

Life, and how far it gets, should follow the scores rather than the class,
and its chemistry should follow the atmosphere, radiation and solvent
rather than the star's letter alone. The highest stage is the latest stage whose
time since habitable (star age minus a 0.05 Gyr onset, plus a pre-main-
sequence delay for M dwarfs) is reached on the planet's pace and whose score
gate passes; Earth's timeline fits `EVOLUTIONARY_TIMELINES["normal"]` (0.5,
1.5, 2.5, 4.0, 4.5 Gyr), while the `fast` and `slow` timelines have no
physical support (rule table in the same document, 7).

## 7. Nebulae

None of the documents cover nebulae. The only link they support is cosmic
ray flux, raised for a system inside a dense cloud by its compressed
heliosphere; the dose research limits that multiplier to about 2.5x and fades
it by 300 g/cm2 of column
([activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md) 4.3).
A separate feasibility study covers disk photoevaporation near O and B stars,
26Al and 60Fe heating from supernova enrichment, extinction of starlight and
cosmic rays near remnants, and turns the results into generation rules
(nebula-and-asteroid-field-classes.md).

## 8. Evidence notes

The equipment-tier and conflict numbers rest on [S] and [C] values (NASA CO2
limit, Armstrong altitudes, the alveolar equation) except the following [R],
to check when paper access is allowed (the research environment could only
read search-result text): NIOSH and OSHA CO2 and CO values; La Rinconada's
pO2; the atmospheric-diving-suit depth rating; the 160 kPa CNS limit;
the diving lower O2 limit (secondary source); Ball and Hallsworth's 73.8 kJ/kg.
The Phi_tech tier values (2.1) and the radiation rows of the tier table are
design defaults [D]. Source lists for the numbers in the two further notes are
in those notes.
