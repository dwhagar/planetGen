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

The research of 2026-10-09 behind the inputs is in two further notes:

- [atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md):
  when a planet keeps an atmosphere, the hot- and cold-zone classes (GEN.91),
  the class windows and life stage against time habitable (GEN.92).
- [activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md):
  stellar XUV and flares, planetary magnetic fields, surface dose, water,
  oceans and ocean chemistry (GEN.86, GEN.87, GEN.88).

Status: the scores (section 2) are settled and built (GEN.84); the rest is
research, 2026-10-09, and decisions marked "Boss" are his, everything else is a
recommendation. Informs: GEN.83 to GEN.92.

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
- **Spin**: planets and moons lock by the tidal-locking time (GEN.104, which
  also adds the spin vector); an unlocked rocky planet's day is log-uniform
  from 8 to 48 h (Kokubo and Genda 2010; Boss 2026-10-10), a gas giant's
  uniform from 8 to 20 h.
- **Rogue planets** already have internal heat flux, a surface regime,
  ice-shell thickness and ocean depth (`rogue_surface.py`); nothing in a
  star system has these.
- **Nebulae**: a system's `surrounding_cloud` is set only when it is
  loaded from the database, so generation can't see it.
- **Life** is decided by class (`HABITABLE_PLANET_CLASSES`), and the life
  chemical by the star's spectral letter. Classes N and Q carry a life
  chemical and a stored timeline but are not in the habitable list, so the
  page hides the timeline.

## 2. The scores (settled, GEN.84)

Boss chose the structure on 2026-10-07 (17:11Z): PHI-4's four domains and
colour tiers for display, with the Xenobiology doc's three scores as the
numbers behind them. `planetgen/physics/habitability.py` is the reference
maths: every constant below is named there with its source, and GEN.89
stores what it computes for every planet and moon.

**PHI-4 (display, for a human visitor).** Each domain's variables have
bands; a variable scores 1 in Blue, falls linearly to 0.75 at the Green
edge, 0.4 at the Yellow edge and 0 at the outer Red edge (in log10 for
pressures and doses). A domain scores its worst variable, and its tier is
Blue at 1, Green from 0.75, Yellow from 0.4, Red below. PHI-4 is the
geometric mean of the four domains (0 if any is 0, the doc's hard gates).

| Variable | Blue | Green | Yellow | Red |
|---|---|---|---|---|
| Pressure, kPa | 50 to 250 | 14.3 to 400 | 6.3 to 10,000 | beyond (0 a decade out) |
| Temperature, C | 0 to 31 | -20 to 45 | -50 to 122 | beyond (0 at -120, 200) |
| Wet bulb | | under 31 C | 31 C or more makes Blue or Green Yellow | |
| Inspired O2, F(P - 6.3) kPa | 8 to 50 | under 8 (a mask works) | 50 to 160 | over 160 |
| CO2, kPa | under 0.93 | to 2.0 | to 5.0 | over 5.0 |
| CO, H2S, SO2 | under the chronic limit | | to ten times the acute limit | beyond |
| Water pH | 6 to 8.5 | 5 to 9.5 | 1 to 11.5 | beyond |
| Water activity | 0.90 or more | 0.75 or more | 0.605 or more | below |
| Chaotropicity, kJ/kg | to 73.8 | | | over 73.8 |
| Dose, Sv/yr | under 0.05 | to 0.1 | to 10 | over 10 (0 at 1000) |

Where the numbers come from: 6.3 kPa is the Armstrong limit; 14.3 kPa is
the least pressure at which a pure-oxygen mask still gives the 8 kPa of
inspired O2 consciousness needs (6.3 + 8); 400 kPa is where nitrogen
narcosis starts; 31 C wet bulb is the uncompensable limit (Vecellio et al.
2022); 122 C is the limit of carbon-based life; gas limits are the
Atmospheric Toxicity table's chronic and acute values, and a filtering
respirator (NIOSH APF 10) makes Yellow out to ten times the acute limit;
CO2 can't be filtered, so its Yellow ends at the acute 5 kPa (scrubbers);
0.605 and 73.8 kJ/kg are the water-activity and chaotropicity limits for
cell division; 50 mSv/yr is the occupational limit and 10 Sv/yr the
Eigen threshold.

**Equipment**: shirtsleeve; breathing mask (low O2 or pressure below Blue);
mask with scrubber (a gas above its chronic limit); sealed suit (pressure
below Green, or any domain Red); full life support with radiation
hardening (Radiation Red).

**PHI_bio (microbial life)** = L_solv · L_chem · L_ener · L_rad, a product
(see 3.1).

- L_solv = solvent factor · f(aw) · exp(-0.5 ((chi - 73.8) / 10)^2) past
  73.8 kJ/kg; f(aw) = 1 / (1 + exp(-50 (aw - 0.605))). Solvent factor:
  water 1, liquid hydrocarbons 0.2, sulfuric acid 0.1, no stable liquid 0.
- L_chem = (prod over C, N, H of tanh(X_j / 0.25 Earth))^(1/3) · P / (P +
  0.1 umol/L).
- L_ener = log10(Phi / 1e-6) / log10(10 / 1e-6), clamped to [0, 1], Phi
  the best energy flux in W/m^2 (light, chemistry or radiolysis).
- L_rad = exp(-0.5 ((D - 10) / 100)^2) above 10 Gy/yr, D the dose where
  the solvent is (under the ice for an ice-sealed ocean).

**PHI_cpx (animal life)** = PHI_bio · exp(-(pCO / 0.01)^2) · exp(-(pCO2 /
2.0)^2) · f(pO2), with f = 1 for ambient O2 from 8 to 50 kPa and a half
Gaussian outside (sigma 2 kPa below, 20 kPa above).

**Phi_tech (human operability with equipment)** = the geometric mean of
four factors (one hard factor can be partly bought off):

- M_press: 0.2 below 6.3 kPa (a sealed base, as on the Moon), 1 to 250
  kPa, exp(-(P - 250) / 2500) above.
- M_therm: 1 / (1 + exp(0.5 (wet bulb - 31))) times a cold side 0.3 + 0.7
  / (1 + exp(-0.2 (T + 50))) (heating a habitat is cheap, so cold never
  makes a base impossible).
- M_isru by the water a base can draw on: liquid 1.0, ice 0.8, hydrated
  minerals or perchlorate brine 0.6, vapour only 0.4, none 0.1.
- M_rad = (0.05 / D)^0.2 above 50 mSv/yr.

**Surface dose**, until GEN.87 adds stellar particles and magnetic
fields: D = 0.45 exp(-X / 22) + 0.05 exp(-X / 200) Sv/yr from cosmic rays,
plus 2 mSv/yr from the ground, X = P0 / g in g/cm^2. This fits the Moon
(0.50), Mars (about 0.25 at 16 g/cm^2), a 0.05 bar world (about 0.08) and
Earth (2.4 mSv/yr). Inside a nebula the cosmic-ray part doubles (section
6).

### Worked examples

Inputs not listed take Earth's: 9.81 m/s^2, C, N and H at one Earth
inventory, 2.3 umol/L phosphate, 80 W/m^2, liquid water at aw 1 and chi
0, relative humidity 0.5, liquid water for a base. Doses are computed
unless given.

| World | Inputs | P / T / C / R | PHI-4 | PHI_bio | PHI_cpx | Phi_tech | Equipment |
|---|---|---|---|---|---|---|---|
| Earth | 101.3 kPa, 15 C, N2 78%, O2 21%, CO2 0.04%, pH 8.1, RH 0.7 | Blue / Blue / Blue / Blue | 1.00 | 0.96 | 0.96 | 1.00 | shirtsleeve |
| Dense CO2 | 300 kPa, 40 C, CO2 95%, N2 5%, pH 6.5, RH 0.8, 40 W/m^2 | Green / Yellow / Red / Blue | 0 | 0.96 | 0 | 0.48 | sealed suit |
| Hycean | 2 MPa, 60 C, H2 90%, He 10%, pH 7, RH 0.99, 20 W/m^2 | Yellow / Yellow / Green / Blue | 0.74 | 0.96 | 0 | 0.02 | sealed suit |
| Photochemical CO | 100 kPa, 10 C, N2 90%, CO2 5%, CO 5%, pH 7.5 | Blue / Blue / Red / Blue | 0.59 | 0.96 | 0 | 1.00 | sealed suit |
| Ice-sealed ocean | vacuum, -170 C, g 1.3, pH 8, 1e-3 W/m^2, surface 0.05 Sv/yr, ocean 2 mSv/yr, ice | Red / Red / Green / Blue | 0 | 0.41 | 0 | 0.47 | sealed suit |
| Mars-analog brine | 0.6 kPa, -60 C, CO2 95%, g 3.71, aw 0.5, chi 90, pH 8, C/N/H 0.05/0.01/0.05, 0.5 umol/L P, 40 W/m^2, ice | Red / Red / Red / Yellow | 0 | 0.00 | 0 | 0.46 | sealed suit |
| Dune world | 90 kPa, 45 C, N2 78%, O2 21%, RH 0.05, aw 0.7, pH 8, C/N/H 0.5/0.5/0.05, vapour | Blue / Green / Yellow / Blue | 0.83 | 0.54 | 0.54 | 0.80 | shirtsleeve |
| Europan ocean | vacuum, -160 C, g 1.31, pH 9, 1e-4 W/m^2, surface 2000 Sv/yr, ocean 2 mSv/yr, ice | Red / Red / Green / Red | 0 | 0.27 | 0 | 0.28 | full life support |

## 3. How the documents' problems were settled

1. PHI_bio is the product of its four likelihoods, not their geometric
   mean: each is a requirement, and the mean forgave a missing solvent
   (the doc's dune world computed to 0.34 against its stated 0.05).
   Phi_tech keeps the doc's geometric mean, since equipment can make up
   for one hard factor in part.
2. The undefined constants are fixed above: k_w 50, X_ref a quarter of
   Earth's inventory, K_P 0.1 umol/L (the doc's 1.0 scored Earth's ocean
   0.70), L_ener log-scaled from 1e-6 to 10 W/m^2 (no single Phi_ref in
   tanh fits both sunlight and radiolysis), sigma_P 2500 kPa, M_isru a
   lookup in place of alpha and E_specific, M_rad's exponent 0.2 in place
   of lambda, solvent factors 1 / 0.2 / 0.1 / 0, and the surface dose
   formula.
3. L_rad keeps 10 Gy/yr; Mars's 0.25 Sv/yr then scores 1, not the table's
   0.65. Mars is held back by its solvent and chemistry instead.
4. One threshold each: CO chronic 0.005 kPa, acute 0.1, and T_metazoa's
   0.01 scale; CO2 chronic 0.93, Green 2.0, acute 5.0; the mask band is
   14.3 to 400 kPa (Green), with 6.3 to 14.3 kPa needing a counterpressure
   suit (Yellow).
5. Other stars (GEN.86, GEN.87): F stars act like G with a shorter
   saturated phase (about 50 Myr); A, B and O stars have no flare term but
   a high UV flux; white dwarfs and giants use their own luminosity and
   age; binaries sum both stars' fluxes at the planet.
6. Nebulae: section 6.

## 4. Inputs to build, in order

| # | Piece | Needs |
|---|---|---|
| 0 | Planet spin and tidal locking (GEN.104); atmosphere retention module ([atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md) 3) | mass, radius, insolation, star age |
| 1 | Mantle redox (reduced, intermediate, oxidized) from mass and differentiation | mass, radius |
| 2 | Atmosphere species: partial pressures of O2, CO2, CO, N2, Ar, H2, H2O, CH4, H2S, SO2, with the free text generated from them | 1, pressure, class |
| 3 | Stellar activity: saturation phase, L_XUV/L_bol, flare and particle-event rates by mass and age ([activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md) 2) | star mass, age, type |
| 4 | Magnetic field from mass, rotation and tidal locking | mass, rotation (the spin vector) |
| 5 | Surface dose from column mass (P0/g), field, cosmic rays (higher inside a compressed heliosphere) and flares, with an ozone-loss flag | 2, 3, 4, nebula containment |
| 6 | Hydrosphere: water fraction, ocean and land fraction, depth or ice shell | zone, temperature, class; reuse `rogue_surface` in a shared module (same document, 5) |
| 7 | Ocean chemistry: acid sulfate, neutral, soda, chloride brine or ice-sealed, with pH, water activity and phosphorus | 1, 2, 6 |
| 8 | The scores, tier and equipment profile | all of the above |

Values are stored in their own columns, not JSON, so they can be searched.

Rows 1 and 2 are built (GEN.85, schema v70, `planetgen/physics/atmosphere.py`).
The redox mean is `IW+3.5 + 2 log10(M / M_earth)` with a sigma of 1 [D],
reduced at or below IW and oxidized from IW+2.5. Each class's mix is a
template read from its old atmosphere text and its solar-system analogue,
scattered log-normally. A reduced mantle under anoxic air moves 60 percent of
CO2, SO2 and H2O to CO and CH4, H2S and H2 (15 percent for intermediate). No
gas exceeds its vapour pressure at the surface temperature (the
`pN2_cap` of the retention document's 3.3, applied to every gas), and the
surface pressure drops by what freezes out. Helium, ammonia and sodium are a
class's unstored balance. `atmosphere.hydrogen_kpa` gives the H inventory for
the water-loss clock. An H2-dominated mix (`x_H2` over 0.5) is never called
"reducing" in the text, so the sub-Neptune and Hycean classes can reuse it
once GEN.90 adds them.

## 5. Classes and life

Known conflicts between today's classes and the research:

- **N** (Venus analog, about 737 K and 9.2 MPa) carries life with no stage
  cap; the thermal ceiling for life is 122 C (395 K). **E** (about 399 K)
  is above it too.
- **K** (about 611 Pa) and **L** (about 2 kPa) are below the Armstrong
  limit (6.3 kPa); L has land animals at a tiny O2 partial pressure.
- **P** (about 9.5 kPa, 204 K) can reach civilization; the research gives
  sub-baric snowballs minimal native life.
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
hyperbaric N2/Ar world.

Life, and how far it gets, should follow the scores rather than the class,
and its chemistry should follow the atmosphere, radiation and solvent
rather than the star's letter alone.

## 6. Nebulae

None of the documents cover nebulae. The only link they support is cosmic
ray flux, raised for a system inside a dense cloud by its compressed
heliosphere. A separate feasibility study covers disk photoevaporation
near O and B stars, 26Al and 60Fe heating from supernova enrichment,
extinction of starlight and cosmic rays near remnants, and turns the
results into generation rules (nebula-and-asteroid-field-classes.md).

## 7. Research record (2026-10-09)

Material from the research pass that the settled maths in section 2 does not
use. It is kept as an alternative with measurements, to be weighed when
GEN.89 is calibrated; nothing here is built.

### 7.1 A Liebig blend for PHI_bio, measured against the Master doc's table

Section 3 item 1 settles on the product of the four likelihoods. Re-evaluating
the Master doc's own worked table (it gives four `L` values and `PHI_bio` per
archetype) [C] shows another rule fits the stated values better:
`PHI = L_min^0.65 * (L_solv L_chem L_ener L_rad)^(0.35/4)`.

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
0.050).

The Master text calls the mean "non-compensatory"; the geometric mean is only
weakly so. The blend lands within 0.1 of every row, so the table was plausibly
made with a Liebig-like rule. The product built in section 2 is stricter (it
scores the Master doc's Europan ocean 0.27 against a stated 0.68). The two are
not a like-for-like test, because this table uses the Master doc's own `L`
values and section 2 fixed different constants; the blend is the option if the
product proves too harsh for ocean worlds.

### 7.2 Alveolar oxygen behind the mask band

Alveolar PAO2 = x_O2 (P - 6.3) - 5.3/0.8 kPa at 37 C (respiratory quotient
0.8) [C]:

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

Section 2 starts the mask band at 14.3 kPa (6.3 plus the 8 kPa of inspired
O2 consciousness needs). The alveolar equation, which also subtracts the CO2
the lungs add, puts a pure-oxygen mask's useful floor nearer 19 to 21 kPa. The
settled 14.3 kPa stands until GEN.89's tests show a case where the difference
matters.

### 7.3 Five equipment tiers

An alternative to the equipment list in section 2: one numbered tier taken from
the worst condition, a pure function of dry ambient values [C].

| Tier | Name | Conditions (all must hold) |
|---|---|---|
| 0 | Shirtsleeve | 50 <= P_tot <= 250 kPa; 16 <= pO2 <= 50; pCO2 < 0.5; pCO < 0.005; pH2S < 0.001; pSO2 < 0.0005; T_dry -20 to 45 C and wet-bulb < 31 C (ordinary weather clothing); dose <= 50 mSv/yr |
| 1 | Mask | Any of: 20 <= P_tot < 50 (O2-enriched mask; pO2 of air too low); 8 <= pO2 < 16; pO2 50 to 160 (diluent mask); pCO2 0.5 to 1; 250 < P_tot <= 400; T 45 to 90 C, or -20 to -50 C with insulation; wet-bulb 31 to 35 C; dose 50 to 100 mSv/yr |
| 2 | Mask and scrubber | Any of: pCO2 >= 1; pCO >= 0.005 (Hopcalite); H2S or SO2 above the chronic limit; 400 < P_tot <= 1,000 kPa (heliox or trimix rebreather); pO2 < 8 with P_tot >= 20 and a hostile base gas |
| 3 | Pressure suit | P_tot < 20 kPa (6.3 to 20: gas-pressurised or mechanical-counterpressure suit plus helmet; below 6.3: full suit); P_tot > 1 MPa (atmospheric diving suit); pO2 > 160; T 90 to 120 C (liquid-cooled garment) or -50 to -120 C (heated suit); wet-bulb > 35 C; dose 0.1 to 1 Sv/yr (shielded habitat, limited surface time) |
| 4 | Full life support | P_tot > 7 MPa (beyond the tested human record, Hydra X at 7.11 MPa); T < -120 or > 120 C; corrosive acids (H2SO4, HF, Cl); dose > 1 Sv/yr (radiation hardening); no usable atmosphere and no suit-compatible thermal range |

The tier is the highest any condition demands. The cold range of -50 to
-120 C, which PHI-4 leaves unassigned, is tier 3 here. Differences from
section 2: the dose rows here are stricter than a rule that passes tier 0
below 1 Sv/yr (which would put Mars at 0.24 Sv/yr and the Moon at 0.52
at shirtsleeve), and the pressure-suit band here is under 20 kPa where
section 2 uses 14.3.

Sources behind the thresholds: aviation practice (100 percent oxygen holds an
equivalent 10,000 ft to about 40,000 ft, a pressure suit above about 50,000
ft, ebullism at 62,000 to 63,000 ft) [S: litfl, Wikipedia "Pressure suit"];
NASA's spacecraft CO2 limit of 3 mmHg (0.4 kPa) over one hour [S: NASA OCHMO
technical brief]; atmospheric diving suits at 1 atm inside to about 300 m
[R]; the human hyperbaric record, Hydra X at 7.11 MPa.

### 7.4 Evidence notes

The numbers above rest on [S] and [C] values except the following [R], to
check when paper access is allowed (the research environment could only read
search-result text): NIOSH and OSHA CO2 and CO values; La Rinconada's pO2; the
atmospheric-diving-suit depth rating; the 160 kPa CNS limit; the diving lower
O2 limit (secondary source); Ball and Hallsworth's 73.8 kJ/kg. Source lists for
the numbers in the two further notes are in those notes.
