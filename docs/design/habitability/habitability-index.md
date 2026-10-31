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
- **Spin**: a rocky planet is born spinning at the rate and tilt the giant impacts that built it leave (Kokubo and Genda 2010, GEN.731:
  a median day of 2 to 3 hours, an isotropic axis), a moon log-uniform from 8 to 48 h, a gas giant uniform from 8 to 20 h; the tide of
  what it orbits then slows the spin for the age of the star (`tidal_spin`), to a lock (one turn per orbit), to a spin-orbit resonance
  such as Mercury's 3:2, or not at all.
- **"Habitable" is the PHI-6 rating, not a class** (Boss, 2026-10-23, GEN.643): a
  world is habitable when its overall PHI-6 colour is Blue or Green
  (`habitability.is_habitable`, the one function the API's `habitable`
  flag, the system page chips, the per-system count `hab_count` and the
  description text all call). The API serves the class flag
  separately as `habitable_class` (the "Habitable class" chip). A request for a
  habitable world (`+habitable_world`) raises the odds of a body that meets its conditions (`generation/habitable_request.py`, GEN.1130) and steers by no class list or score, so the world it gives may not score Blue
  or Green; `StarSystem.habitable_class_count` is a report column only, and
  `habitability.habitable_request_result` says whether the finished world met the
  request ("aim"), settled for a Yellow with Suitability at least Yellow
  ("fallback"), or missed.
- **Rogue planets** already have internal heat flux, a surface regime,
  ice-shell thickness and ocean depth (`rogue_surface.py`); nothing in a
  star system has these. A terrestrial rogue is scored from its surface too (`rogue_habitability.py`, [rogue-planet-surface.md](../systems/rogue-planet-surface.md)).
- **Nebulae**: a system's `surrounding_cloud` is set only when it is
  loaded from the database, so generation can't see it.
- **Life** is decided by class (`HABITABLE_PLANET_CLASSES`), and the life
  chemical by the star's spectral letter. Classes N and Q carry a life
  chemical and a stored timeline but are not in the habitable list, so the
  page hides the timeline.

## 2. The scores (settled, GEN.84)

Boss chose the structure on 2026-10-07 (17:11Z): PHI-4's domains and
colour tiers for display, with the Xenobiology doc's three scores as the
numbers behind them. GEN.207 (2026-10-15) made the four domains six: Chemistry
split into Air and Water, and a long-term Suitability was added. `planetgen/physics/habitability.py` is the reference
maths: every constant below is named there with its source, and GEN.89
stores what it computes for every planet and moon.

**PHI-6 (display, for a human visitor; called PHI-4 until GEN.207 made the four domains six).** Each domain's variables have
bands; a variable scores 1 in Blue, falls linearly to 0.75 at the Green
edge, 0.4 at the Yellow edge and 0 at the outer Red edge (in log10 for
pressures and doses). A domain scores its worst variable, and its tier is
Blue at 1, Green from 0.75, Yellow from 0.4, Red below. The six domains are
Pressure, Temperature, Air Chemistry (is the air breathable: oxygen, CO2, CO,
H2S, SO2), Water Chemistry (is the water drinkable: pH, water activity,
chaotropicity; a world with no liquid water scores 0, Red), Radiation and
Suitability. PHI-6 (still stored as `phi4`) is their geometric mean held under
the colour gates below.

| Variable | Blue | Green | Yellow | Red |
|---|---|---|---|---|
| Pressure, kPa | 50 to 250 | 14.3 to 400 | 6.3 to 10,000 | beyond (0 a decade out) |
| Temperature, C | 0 to 31 | -20 to 45 | -50 to 122 | beyond (0 at -120, 200) |
| Wet bulb | | under 31 C | 31 C or more makes Blue or Green Yellow | |
| Inspired O2, F(P - 6.3) kPa | 8 to 50 | under 8 (a mask works) | 50 to 160 | over 160 |
| CO2, kPa | under 0.93 | to 2.0 | to 5.0 | over 5.0 |
| CO, H2S, SO2 | under the chronic limit | | to ten times the acute limit | beyond |
| Water pH (Water) | 6 to 8.5 | 5 to 9.5 | 1 to 11.5 | beyond |
| Water activity | 0.90 or more | 0.75 or more | 0.605 or more | below |
| Chaotropicity, kJ/kg | to 73.8 | | | over 73.8 |
| Dose, Sv/yr | under 0.02 | to 0.1 | to 1 | over 1 (0 at 1000) |

Where the numbers come from: 6.3 kPa is the Armstrong limit; 14.3 kPa is
the least pressure at which a pure-oxygen mask still gives the 8 kPa of
inspired O2 consciousness needs (6.3 + 8); 400 kPa is where nitrogen
narcosis starts; 31 C wet bulb is the uncompensable limit (Vecellio et al.
2022); 122 C is the limit of carbon-based life; gas limits are the
Atmospheric Toxicity table's chronic and acute values, and a filtering
respirator (NIOSH APF 10) makes Yellow out to ten times the acute limit;
CO2 can't be filtered, so its Yellow ends at the acute 5 kPa (scrubbers);
0.605 and 73.8 kJ/kg are the water-activity and chaotropicity limits for
cell division; the dose bands are ICRP's (GEN.148, Boss's defaults of
2026-10-10): Blue to the 20 mSv/yr five-year average, Green to 100 mSv/yr,
Yellow to 1 Sv/yr (the threshold of radiation sickness; shelter or shielding
needed) and Red beyond, which closes the gap the PHI-4 doc left between
100 mSv/yr and 10 Gy/yr. 50 mSv in any one year is ICRP's single-year limit
and where `M_rad` starts charging for shielding. The Eigen threshold of 10
Gy/yr belongs to `L_rad` (microbes), not to a person's tiers.

**Ultraviolet in the Radiation domain (GEN.445).** The Radiation score is the lower of the dose's and the ground ultraviolet's. The UV is the DNA-damaging (Setlow) dose at the ground as a multiple of Earth's (`World.uv_index`, from the stored `uv_surface_index`, which `radiation.generate_dose` derives from the two-stream spectrum under the world's own air, ozone, haze and cloud; radiation-model.md section 7). Blue below 3 times Earth's (where sun-cover clothing starts), Green to 30 (a full-cover suit), Yellow to 300 (the unfiltered light of a world with oxygen and no ozone layer, 280 to 356 times Earth's: shelter from the sun), Red beyond, 0 at 3000 (`protection_limits` rows `uv_cover_index`, `uv_suit_index`, `uv_lethal_index`, `uv_outer_index`; the last two are engineering edges from the research table, not published limits). A world with no UV index scores on its dose alone. A body with no stored dose reads the cosmic-ray tables of the radiation model (`World.bare_cosmic_ray_dose_sv_yr`); the earlier second curve is gone.

**Habitable request fallback and food (Boss, 6:54 PM EDT 2026-10-23).** A request for a habitable world aims for Blue or Green overall; if it must, it accepts Yellow as technically habitable so long as the sustainability factor (read as Suitability) is at least Yellow. Food comes from the Flora and Fauna scores: Flora is always good food, both are better, Fauna includes all animal life (GEN.653).

**The four life scores (Boss, 6:14 PM EDT 2026-10-23).** Microbial, Flora, Fauna and Sentient, each 0 (none) to 9 (abundant), from the evolutionary stage and the surface conditions, shown as four small icons with a digit next to the PHI-6 icons. They replace the microbial, complex-life and human-operability scores (phi_bio, phi_cpx, phi_tech; GEN.649); section 12 has the rules, and `/help/habitability` says them to the reader (DOC.77).

**Habitable defined (Boss, 6:02 PM EDT 2026-10-23).** A world is habitable when its overall score is Blue or Green; Yellow and Red are not. One predicate (GEN.643) serves the sector counts, the markers, the lists, the search and the descriptions. Boss also wants the PHI-6 tied into the planet classification; the proposal is GEN.644. Decided (Boss, 6:09 PM EDT): classes are conditions-based and assigned after the body is complete, from the atmosphere composition and pressure, the surface temperature, the terrestrial and the gas composition (GEN.645 to GEN.648).

**The gates (GEN.207; Boss, 2026-10-11 and 2026-10-14).** Importance, most
first: Radiation, Temperature, Pressure, Air, Water, Suitability. The colour
gates come first, then the geometric mean inside them
(`habitability.overall_score`):

- Radiation and Temperature cap the overall colour at their own: Yellow caps
  at Yellow, Red at Red.
- Pressure: Red makes the overall Red and Yellow caps it at Yellow, unless a
  full suit is needed anyway (Temperature, Radiation or Air Red), when it
  changes nothing about the overall (it is still listed among the reasons for
  the suit).
- Air in Yellow or Red caps the overall at Yellow; Red air alone cannot make
  it Red. Water likewise.
- Suitability Red drops the overall to Yellow, never to Red.
- Blue overall needs all six Blue.

The ones that cannot make the overall Red (Air, Water, Suitability, and
Pressure when a suit is needed anyway) enter the mean no lower than 0.4 (the
Yellow edge), so a dead one does not pull the mean to 0.

**Suitability (GEN.207)** is long-term survival on the surface: the lowest of
three rungs, water within reach, food and shelter, each read off
Exoplanetary Habitability.md Table 2, then held under a ceiling the other
five set (`habitability.suitability_ceiling`): Red if Radiation, Temperature,
Pressure or Water is Red (a world with no liquid water has Red Water and so
Red Suitability); Yellow if any of the five is Yellow or Air is Red; Green only
if the five are Green or Blue; Blue only if all five are Blue.

| Rung | Blue (1) | Green (0.9) | Yellow | Red |
|---|---|---|---|---|
| Water within reach | liquid at the surface and Blue Water | liquid at the surface | ice, or an ocean under ice (0.6) | water in minerals (0.25), vapour only (0.2), none (0) |
| Food (native life) | a full biosphere of forests and animals (1) | the first plants and animals (0.85) | microbial mats and plankton (0.5) | none (0.2); life only in an ocean under ice counts as none |
| Shelter | the surface dose is Blue | Green | Yellow (caves or burial) | Red (deep lava tubes or buried habitats) |

The food stage is the latest milestone of the planet's evolutionary timeline
(abiogenesis to complex cells are microbial, multicellularity is the first
plants and animals, a technological civilization is a full biosphere), or
microbial for life with no timeline (a moon, a world outside the habitable
zone) and none when no life chemistry was drawn. The numbers are planetGen's,
one per tier, to be tuned with sample worlds. Life is drawn after the score, so
generation scores again once it is known (`habitability_world.rescore`).

**Equipment**: ideal (called "shirtsleeve" until Boss renamed it, 2026-10-10 08:08Z); breathing mask (low O2 or pressure below Blue);
mask with scrubber (a gas above its chronic limit); environment suit (pressure
below Green, any domain Red, a temperature past the Green heat or cold edge, or
oxygen above the chronic limit, which is where the description asks for one); full life support with radiation
hardening (Radiation Red). The kit is about the air, the cold or heat and the
radiation; Water Chemistry and Suitability are about drinking and eating and do
not change it.

**PHI_bio (microbial life)** = the Liebig blend of L_solv, L_chem, L_ener
and L_rad: `L_min^0.65 · (L_solv · L_chem · L_ener · L_rad)^(0.35/4)` (GEN.148;
section 7.1 has the measurements against the Master doc's table; section 3
item 1 records the product it replaced).

- L_solv = solvent factor · f(aw) · exp(-0.5 ((chi - 73.8) / 10)^2) past
  73.8 kJ/kg; f(aw) = 1 / (1 + exp(-50 (aw - 0.605))). Solvent factor:
  water 1, liquid hydrocarbons 0.2, sulfuric acid 0.1, no stable liquid 0.
- L_chem = (prod over C, N, H of tanh(X_j / 0.25 Earth))^(1/3) · P / (P +
  0.1 umol/L).
- L_ener = log10(Phi / 1e-6) / log10(10 / 1e-6), clamped to [0, 1], Phi
  the best energy flux in W/m^2 (light, chemistry or radiolysis).
- L_rad = exp(-0.5 ((D - 10) / 100)^2) above 10 Gy/yr, D the absorbed dose
  where the solvent is (under the ice for an ice-sealed ocean). D is in
  **Gy/yr**; the person's bands are in **Sv/yr**, and the two are different
  quantities (Mars: 0.077 Gy/yr, 0.236 Sv/yr, so about 3.06 Sv per Gy for a
  cosmic-ray dose; `habitability.SV_PER_GY`). The maths keeps them as two
  numbers by hand because astropy 6.0.1 has no `Gy` or `Sv` unit (and astropy
  is no longer a runtime dependency).

**Microbial dose tolerance** (D10, the dose for a 90 percent kill; for any
biology table and for `LETHAL_HABITAT_DOSE_GY_YR`): E. coli about 0.3 kGy,
bacterial spores 2.1 kGy, *Deinococcus radiodurans* 5 to 10 kGy (10.4 kGy on
average; use 5 to 10), a tardigrade's LD50 1.1 to 6.2 kGy
([activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md) 4.6).
`L_rad`'s 10 Gy/yr is Eigen's threshold for the origin of replicating
molecules, a thousand times under these. A compact object's planet is rated
at the Deinococcus figure, 10 kGy/yr, so it keeps no life.

**PHI_cpx (animal life)** = PHI_bio · exp(-(pCO / 0.01)^2) · exp(-(pCO2 /
2.0)^2) · f(pO2), with f = 1 for ambient O2 from 8 to 50 kPa and a half
Gaussian outside (sigma 2 kPa below, 20 kPa above).

**Phi_tech (human operability with equipment)** was the Xenobiology document's third number (the geometric mean of a
pressure, a thermal, a water-supply and a radiation-shielding factor). Boss removed it on 2026-10-23 (GEN.649): the
equipment tier and the six domains already say what a person needs, and the number was shown beside two life scores it
had nothing to do with. Its four factors, the constants behind them (sigma_P, the cold floor, the water-source
lookup, the shielding exponent) and the column are gone.

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

| World | Inputs | P / T / Air / Water / R / Suit | PHI-6 | PHI_bio | PHI_cpx | Equipment |
|---|---|---|---|---|---|---|
| Earth | 101.3 kPa, 15 C, N2 78%, O2 21%, CO2 0.04%, pH 8.1, RH 0.7, a full biosphere | Blue x6 | 1.00 | 0.97 | 0.97 | ideal |
| Dense CO2 | 300 kPa, 40 C, CO2 95%, N2 5%, pH 6.5, RH 0.8, 40 W/m^2, no life | Green / Yellow / Red / Blue / Blue / Red | 0.69 | 0.97 | 0 | environment suit |
| Hycean | 2 MPa, 60 C, H2 90%, He 10%, pH 7, RH 0.99, 20 W/m^2, microbes | Yellow / Yellow / Green / Blue / Blue / Yellow | 0.73 | 0.97 | 0 | environment suit |
| Photochemical CO | 100 kPa, 10 C, N2 90%, CO2 5%, CO 5%, pH 7.5, no life | Blue / Blue / Red / Blue / Blue / Red | 0.74 | 0.97 | 0 | environment suit |
| Ice-sealed ocean | vacuum, -170 C, g 1.3, pH 8, 1e-3 W/m^2, surface 0.05 Sv/yr, ocean 2 mSv/yr, ice, microbes under it | Red / Red / Green / Blue / Green / Red | 0 | 0.53 | 0 | environment suit |
| Mars-analog brine | 0.6 kPa, -60 C, CO2 95%, g 3.71, aw 0.5, chi 90, pH 8, C/N/H 0.05/0.01/0.05, 0.5 umol/L P, 40 W/m^2, ice, no life | Red / Red / Green / Red / Yellow / Red | 0.3999 | 0.01 | 0 | environment suit |
| Dune world | 90 kPa, 45 C, N2 78%, O2 21%, RH 0.05, aw 0.7, pH 8, C/N/H 0.5/0.5/0.05, vapour, first plants | Blue / Green / Blue / Yellow / Blue / Red | 0.7499 | 0.64 | 0.64 | ideal |
| Europan ocean | vacuum, -160 C, g 1.31, pH 9, 1e-4 W/m^2, surface 2000 Sv/yr, ocean 2 mSv/yr, ice, microbes under it | Red / Red / Green / Green / Red / Red | 0 | 0.40 | 0 | full life support |

Recomputed for GEN.207 (six domains and the gates): the Water column is the old
Chemistry's water inputs, the Air column its gas inputs, and Suitability is new;
Red air or Red Suitability alone now gives a Yellow overall instead of 0 (the dense CO2 world
0 to 0.69, the photochemical CO world 0.59 to 0.74); the Mars-analog brine
has Red Pressure and Temperature (0.3999, the top of Red). Recomputed for GEN.148 (the Liebig blend, the 20 mSv/yr Blue edge): the
ice-sealed ocean's surface dose of 0.05 Sv/yr is now Green (it was Blue);
PHI_bio rose where one likelihood was short (the Europan ocean 0.27 to 0.40,
the dune world 0.54 to 0.64, the ice-sealed ocean 0.41 to 0.53) and
Earth 0.96 to 0.97. The Master doc's own stated values for the same
archetypes are in 7.1; a test (`test_habitability.py`) pins this table.

## 3. How the documents' problems were settled

1. PHI_bio was first built as the product of its four likelihoods, not
   their geometric mean: each is a requirement, and the mean forgave a
   missing solvent (the doc's dune world computed to 0.34 against its stated
   0.05). GEN.148 replaced the product with the Liebig blend (7.1), which
   keeps that gate (a term at 0 gives 0) and lands within 0.1 of the Master
   doc's stated value for every archetype.
2. The undefined constants are fixed above: k_w 50, X_ref a quarter of
   Earth's inventory, K_P 0.1 umol/L (the doc's 1.0 scored Earth's ocean
   0.70), L_ener log-scaled from 1e-6 to 10 W/m^2 (no single Phi_ref in
   tanh fits both sunlight and radiolysis), solvent factors 1 / 0.2 / 0.1 / 0, and the surface dose
   formula.
3. L_rad keeps 10 Gy/yr, read in Gy (not Sv); Mars's 0.077 Gy/yr (0.236
   Sv/yr) then scores 1, not the table's 0.65. Mars is held back by its
   solvent and chemistry instead.
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
class's unstored balance. The H inventory for
the water-loss clock is read from the stored partial pressures (H2, H2O, 2 CH4, H2S) where that clock is built. An H2-dominated mix (`x_H2` over 0.5) is never called
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

**Built first (GEN.206, 2026-10-14).** The descriptions now follow the scores, in
`generation/life.py` (`life_ceiling`, `apply_life_data`):

- `phi_bio` below `tuning.HABITABILITY_LIFE_SCORE_MIN` (0.005, what shows as 0.00):
  the world carries no life chemistry, speed, spectrum or timeline. The class texts
  claim no life (GEN.221 took "bacterial life", "simple life" and "with vegetation"
  out of F, G and L), so nothing in them needs replacing.
- `phi_cpx` below the same line with microbial life allowed: the timeline stops at
  "photosynthesis" (`get_evolutionary_timeline(..., max_stage=)`, the lower of this
  and the class ceiling `PLANET_CLASS_MAX_LIFE_STAGE`).
- A body with no scores (a gas giant, or one saved without them) keeps the old class
  rules. Of 31 bodies of a habitable class in 120 generated systems, 14 score no
  microbial life, 6 allow microbial life only and 11 are unlimited, so fewer worlds
  carry life than before (a class used to be enough).

The chemistry itself (which pigment, which solvent) still follows the star's letter;
that, and a stage that follows the scores rather than the star's age, are left.

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

Section 3 item 1 first settled on the product of the four likelihoods (GEN.84);
GEN.148 built this blend in its place (`LIEBIG_MIN_WEIGHT` 0.65). Re-evaluating
the Master doc's own worked table (it gives four `L` values and `PHI_bio` per
archetype) [C] shows another rule fits the stated values better:
`PHI = L_min^0.65 * (L_solv L_chem L_ener L_rad)^(0.35/4)`.

| World | Stated `PHI_bio` | Geometric mean | Min | Blend, weight 0.65 |
|---|---|---|---|---|
| Earth | 0.98 | 0.985 | 0.98 |
| Dense CO2, 3 bar | 0.82 | 0.934 | 0.880 | 0.90 |
| Hycean 20 bar H2 | 0.74 | 0.855 | 0.720 | 0.76 |
| Stalled CO | 0.79 | 0.914 | 0.860 | 0.88 |
| Ice-sealed ocean | 0.12 | 0.378 | 0.10 |
| Cryo-brine (Mars) | 0.18 | 0.353 | 0.150 | 0.20 |
| Desiccated dune | 0.05 | 0.337 | 0.020 | 0.05 |
| Europan ocean | 0.68 | 0.774 | 0.64 |


RMS error against the stated values: geometric mean 0.169, product 0.157, min
0.057, blend (weight 0.65) 0.046 (weights 0.55 to 0.75 all give 0.046 to
0.050).

The Master text calls the mean "non-compensatory"; the geometric mean is only
weakly so. The blend lands within 0.1 of every row, so the table was plausibly
made with a Liebig-like rule. The product built in section 2 is stricter (it
scores the Master doc's Europan ocean 0.27 against a stated 0.68). The two are
not a like-for-like test, because this table uses the Master doc's own `L`
values and section 2 fixed different constants. The blend is built; the
product it replaced scored ocean worlds too harshly.

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

### 7.3 Five equipment tiers (superseded)

Superseded by the survival equipment line of GEN.655 (`physics/equipment.py`, docs/design/survival-equipment-fields.md):
the one numbered tier is gone, and the table below is kept as the research it was built from, not as a rule. One numbered
tier taken from the worst condition, a pure function of dry ambient values [C].

| Tier | Name | Conditions (all must hold) |
|---|---|---|
| 0 | Ideal | 50 <= P_tot <= 250 kPa; 16 <= pO2 <= 50; pCO2 < 0.5; pCO < 0.005; pH2S < 0.001; pSO2 < 0.0005; T_dry -20 to 45 C and wet-bulb < 31 C (ordinary weather clothing); dose <= 50 mSv/yr |
| 1 | Mask | Any of: 20 <= P_tot < 50 (O2-enriched mask; pO2 of air too low); 8 <= pO2 < 16; pO2 50 to 160 (diluent mask); pCO2 0.5 to 1; 250 < P_tot <= 400; T 45 to 90 C, or -20 to -50 C with insulation; wet-bulb 31 to 35 C; dose 50 to 100 mSv/yr |
| 2 | Mask and scrubber | Any of: pCO2 >= 1; pCO >= 0.005 (Hopcalite); H2S or SO2 above the chronic limit; 400 < P_tot <= 1,000 kPa (heliox or trimix rebreather); pO2 < 8 with P_tot >= 20 and a hostile base gas |
| 3 | Pressure suit | P_tot < 20 kPa (6.3 to 20: gas-pressurised or mechanical-counterpressure suit plus helmet; below 6.3: full suit); P_tot > 1 MPa (atmospheric diving suit); pO2 > 160; T 90 to 120 C (liquid-cooled garment) or -50 to -120 C (heated suit); wet-bulb > 35 C; dose 0.1 to 1 Sv/yr (shielded habitat, limited surface time) |
| 4 | Full life support | P_tot > 7 MPa (beyond the tested human record, Hydra X at 7.11 MPa); T < -120 or > 120 C; corrosive acids (H2SO4, HF, Cl); dose > 1 Sv/yr (radiation hardening); no usable atmosphere and no suit-compatible thermal range |

The tier is the highest any condition demands. The cold range of -50 to
-120 C, which PHI-4 leaves unassigned, is tier 3 here. Differences from
section 2: the dose rows here are stricter than a rule that passes tier 0
below 1 Sv/yr (which would put Mars at 0.24 Sv/yr and the Moon at 0.52
at ideal), and the pressure-suit band here is under 20 kPa where
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

## 8. The score of every planet and moon (GEN.89, built)

`planetgen/physics/habitability_world.py` turns a rocky body's stored values
into the `World` that section 2's maths reads and stores 19 columns on
`planets` and `moons` (schema v76; a gas giant has none). It runs last in a
body's generation and again when `refresh_containment` changes the dose, from
the row alone. Every rule below is a planetGen default, not a source's number.

| Input | Rule |
|---|---|
| Pressure, temperature, gravity | As stored (`atmospheric_pressure_pa`, `surface_temperature_k`, `gravity_g`). |
| Gases | Each of O2, CO2, CO, H2S and SO2 as its partial pressure over the total; the balance gas only counts toward the total. |
| Relative humidity | `p_h2o_kpa` over water's saturation pressure at the surface temperature (Magnus formula, clamped to -45 to 60 C), limited to 5 to 99%. |
| Solvent | Water for a surface ocean, a Hycean ocean or an ocean under ice; liquid hydrocarbon when CH4 is at 98% of its vapour pressure and the surface is at most 190 K (a Titan); else none, which zeroes `l_solv`. Sulfuric acid is not generated yet. |
| pH, water activity | The ocean's stored `ocean_ph` and `water_activity`. A chloride brine is chaotropic: 200 kJ/kg per unit of water activity below 1 (70 to 140 kJ/kg across its 0.65 to 0.30), past the 73.8 kJ/kg limit for the saltiest. |
| Inventories (C, N, H) | Hydrogen is the water, `water_mass_fraction` x mass over Earth's 1.4e21 kg of ocean. Carbon and nitrogen are at least as much (volatiles arrive with water) or the air's own carbon (CO2 + CO + CH4, against Earth's 7.5e19 kg) or nitrogen (N2, against 4.0e18 kg) when more. |
| Phosphate | The ocean's `phosphorus` class: high 1000, limited 2.3 (Earth's), starved 0.01 umol/L; a methane sea takes 2.3 (the term has no non-aqueous model); no liquid, 0. |
| Energy | A surface liquid: the greater of 6% of the starlight (1361 W/m2 gives Earth's 80 of usable light) and 1% of the body's heat flow (Earth's 0.09 W/m2 gives 9e-4, inside the radiolytic 1e-4 to 1e-3). An ocean under ice: the heat flow's 1% only, and the heat flow of an icy moon is its own plus its tide (`tidal_flux_w_m2`, kept since v137; Europa's 4.6e-3 W/m2 is about fifty times the radiogenic 0.09 of an Earth-like core scaled to its size). Stored as `energy_flux_w_m2`. |
| Dose | `surface_dose_msv_yr` in Sv/yr; an ocean under ice sees only `dose_ground_msv_yr`. |
| Water source for a base | Liquid at a surface ocean, ice for an icy surface or ice lid, vapour for a steam world, hydrated minerals for a dry body with 1e-5 or more water by mass, else vapour in the air, else none. |

**Equipment.** The survival equipment line (GEN.655, `physics/equipment.py`,
docs/design/survival-equipment-fields.md) lists what a human visitor needs to
stand on the surface, in four slots (body protection, air filter or scrubber,
air supply, radiation protection), each item required or recommended, and only
what is needed: a world that needs nothing shows no line. It replaces the old
`equipment_tier` (0 to 4), which is gone. The gear follows the physics wherever
that is stricter than the score's bands; the PHI-6 scores are not changed by it.
Its stored levels are the `eq_*` columns, and the name and the short code
(`B4p F100 S2 [R2]`) are derived from them. Each domain stores its score (`phi4_<domain>`) and its colour
(`tier_<domain>`: 0 Blue to 3 Red); `hab_note` lists the domains short of
Blue.

**Compact hosts.** A planet of a pulsar, a non-pulsing neutron star or a black
hole is rated as if its surface took 1000 Sv/yr, the outer Red edge of the
radiation domain: Red radiation, a hardened habitat required, microbial and complex life
0, and a `hab_note` naming the host. The rule is by the host's class, so an
X-ray binary (a normal star with a compact companion; none is generated yet)
will need its companion's class tested where `generate` and the dose refresh
find the host. The dose refresh leaves such a system's scores alone.

**Open question for Boss (default stands).** GEN.89's note asks whether to
store the two flare numbers or derive the rest; GEN.86 already stores them on
the star, and everything else the note lists is stored per planet.

## 9. One table of limits and one paragraph of description (GEN.208, built)

`planetgen/physics/protection_limits.py` holds every limit of the human visitor in
one table (`LIMIT_ROWS`): the pressure, oxygen, temperature, wet-bulb and dose edges
that the score's bands (`habitability.PRESSURE_BAND`, `TEMPERATURE_BAND`,
`DOSE_BAND`, `ARMSTRONG_KPA`, `INSPIRED_O2_*`, `WET_BULB_CRIT_C`, `COLD_LIMIT_C`)
are built from, and that the gear logic in `planetgen/physics/protection.py` reads.
A later per-species limit replaces a row in this table and nothing else; the score,
the gear and the description cannot disagree.

Where the table the code had (`built`) and the research notes (`notes`) gave two
values, Boss's rule (2026-10-14) is that the more permissive one wins. The six pairs
decided, the value in use first:

| Limit | In use | Other | Reading |
|---|---|---|---|
| Mask minimum pressure | 14.3 kPa | 20 kPa | A pure-oxygen mask still gives 8 kPa inspired O2 down to it. |
| No-gear minimum pressure | 50 kPa | 70 kPa | Lower edge of gear-free air. |
| Hyperbaric suit maximum | 10,000 kPa | 5,000 kPa | Higher edge of a hyperbaric suit. |
| Inspired O2 minimum | 8 kPa | 19.5 kPa | The consciousness limit, not the 19.5 kPa working limit. |
| Cold, insulated suit to | -50 C | -40 C | Colder edge of an insulated suit. |
| Wet-bulb limit | 31 C | 30 C | The higher uncompensable wet-bulb value. |

The one new row, a heat-shielded suit past 140 C, comes from the notes alone. The
Armstrong limit (6.3 kPa against 6.27) is the same number to the notes' rounding.

`protection.rating` gives the description of a rocky body as data (GEN.216, UX.104): the overall value and name, the gear
a visitor needs, and one row per category with its colour, score and one short sentence. The describer
(`describe/build.py`) puts it in the body's Document and the renderers write it. It opens "PHI-6 habitability is 0.00 (Red),
and a visitor needs a pressure suit + dust-tight seals." (the required items of the equipment line, then "and should add" the recommended ones; a world that needs nothing opens "PHI-6 habitability is 1.00 (Blue)." with no gear clause) The six categories follow, each as "Name: Colour (value)." and one short sentence with
the number and the gear its category asks for ("The air pressure is 0 kPa. Required: hardened pressure suit.", UX.153, UX.163). Air with no atmosphere reads
"There is no atmosphere (0 kPa)." (GEN.215), Radiation names its band, and Suitability says which of food, water and shelter
are readily available. Gear is named in short words (`SHORT_GEAR`); the first line is the equipment line itself (`protection.gear_clause`), so the
headline and the stored levels cannot disagree. A planet of a pulsar, a neutron star or a black hole says its dose is rated lethal.

## 10. What a species can live in (POP.12, built without storage)

A species' environment limits are section 9's table with some rows replaced
(`planetgen/physics/species_limits.py`): `overrides`, `{limit name: value}`,
empty for the human limits. `habitability.bands(overrides)` builds the pressure,
temperature, dose and oxygen bands from them, and `habitability.phi4`,
`core_scores`, `equipment.compute` and the domain functions take those bands (the human
`HUMAN` by default), so the six domains are scored by the same code with the
species' edges. Nothing else about the score changes: the Water and Suitability
rungs, the toxic gases and PHI_bio are the human visitor's.

**Generated from the homeworld** (Boss's default of 2026-10-15: every species
passes its own homeworld). `species_limits.generated(world)` leaves a band alone
when the human limits already hold the homeworld inside Blue with a margin
(`tuning.SPECIES_LIMIT_MARGIN`). Past it, the side the homeworld lies beyond moves
as a whole, so the Green and Yellow widths stay the human ones: pressure edges
scale by the factor that puts the homeworld 1.5 times inside the Blue edge (low
side or high side), temperature edges move by the degrees that put the cycle's
extreme 10 degrees inside, the wet-bulb limit sits 2 degrees above the
homeworld's, all the dose edges scale so the homeworld's dose is a third of the
Blue edge, and the inspired oxygen limits follow the homeworld's air (the
minimum a half of it, the chronic limit a quarter above it, the acute limit in
the human ratio). A species with an Earth-like homeworld has no overrides.

**What could live where.** `species_limits.assess(world, overrides, tech_level)`
scores the six domains. A species could live there **unaided** when no domain is
Red. It could live there **with equipment** when Water is not Red and its tech
level reaches the highest of the required items its limits ask for (`equipment.tech_needed`;
`tuning.SPECIES_EQUIPMENT_TECH` gives a level per item); a species with no tech
level has no kit. A base brings its own food and shelter, so Suitability does not
count for the equipped list. The lists are `db/species_fit.py`:
`species_for_body` and `species_for_system` (the species that could live on a
body), and `bodies_for_species` (the worlds of a system or sector one species
could live on), through `GET /api/planets|moons/<id>/habitable-species`,
`/api/systems/<id>/habitable-species` and `/api/species/<id>/habitable-worlds/...`;
the system page shows the first as "Species that could live here" and the species
page shows its limits.

**Storage.** `species_fit.limits_for_species` is the one function that says what a
species' limits are. It generates them from the homeworld now; DB.30 adds the
columns, and `species_fit.stored_overrides` then returns the stored values, which
win over the generated ones. A species whose homeworld has no habitability score
has the human limits. POP.13 (the species editor) will write the stored values and
should run `species_limits.validate` on them.

### 10.1 A species is only placed where it can live (POP.14)

`species_fit.check_placement(conn, species_id, kind, body_id, habitat=False)` is the one check every way of putting
a species on a planet or moon calls before it writes. It passes when the species could live there unaided (no Red
domain against its own limits). With `habitat` true, for a placement that is a habitat or a vessel and brings its own
kit, it also passes when the species could live there with the equipment its tech level can make (section 10).
Otherwise it raises `PlacementError`, whose message names the failing conditions in plain words with no numbers
("The Terrans could not live on Kepler II: it is too hot; the radiation is far beyond what it can survive.") and
whose `failing` lists the Red domains. A gas giant fails (no surface). A body with no habitability score cannot be
judged and passes.

A species on its own homeworld is a native there: its Suitability is not judged (`assess(..., native=True)`), because
the water and food of its world are its own. Any other placement is judged on all six domains.

Today a species has one placement, its homeworld, made by the population pass (`population.scan_life_worlds`). It
calls `species_fit.check_homeworld` first and makes no species on a world whose own limits cannot hold it (its water
is undrinkable); the planet is passed over like any other without a civilization. A later settlement (POP.11) and the
species editor (POP.13) call `check_placement` before they write, and the editor runs it on save and shows the failing
domain. Facilities are not placements of a species: a facility's affiliation says who runs it, not who lives there.

`python3 -m planetgen.cli.species_check` is the consistency report: it lists every species placed outside its limits
(exit 1) or says all are inside (exit 0). While the limits are generated from the homeworld they follow it, so the
report finds nothing; it matters once the limits are stored (DB.30, POP.13) and one is edited or its world changes.
A species it lists is fixed by regenerating its limits or the placement; there is no tolerance flag.


## 11. Where the rating is shown, and the zone (DOC.75)

The in-app text for all of this is `/help/habitability` (under the Help menu, UX.180; every chip links to its section, UX.181; a body's symbols and life scores link to the Habitability and Life groups of its own page, UX.197) (`web/lib/habitability_ref.py`, `templates/habitability.html`); its numbers are read from the code named in sections 2 and 9, so it cannot disagree with the score. The seven rating symbols (overall first, then the six domains; UX.101) are drawn on the Planets and Moons lists (`/lists/planet`, `/lists/moon`, a Habitability column that is not sortable; UX.100), under the title of a planet or moon page, in the rows of a system page and in the Search panels. A body with no score shows none; a gas giant has none. The four life scores (0 to 9, the cluster of life shapes beside the symbols) are described on the same page and in GEN.649. The sector colouring (MAP.143), the habitable-world markers (MAP.219, MAP.220) and the radiation terms (GEN.445) are added to the page as each lands.

The habitable zone is a guide, not a rating. Its inner edge is the runaway greenhouse flux and its outer edge the maximum greenhouse flux of Kopparapu et al. (2014, Table 1), each a polynomial in the star's temperature minus 5780 K, valid from 2600 to 7200 K (a star outside takes the nearest end's edges); for the Sun they are 0.95 and 1.68 AU (`physics/orbits.calculate_habitable_zone`, GEN.146). The planet-mass dependence of the inner edge and the early runaway history are GEN.639.

## 12. The four life scores (GEN.649, built)

Boss, 2026-10-23: the line "Microbial 0.00 . Complex life 0.00 . Human operability 0.52" becomes four whole numbers, Microbial,
Flora, Fauna and Sentient, each 0 (none) to 9 (abundant), "based on the evolutionary path and surface conditions", and the
page shows them as small icons with a number beside the PHI-6 symbols. `physics/life_scores.py` makes them; they are stored on
every scored planet and moon (`life_microbial`, `life_flora`, `life_fauna`, `life_sentient`, schema v150, NULL for a gas giant
or a body with no score) and replace the stored `phi_bio`, `phi_cpx` and `phi_tech`. The two likelihoods remain as in-memory
values (they decide, in `life.life_ceiling`, what the life draw may claim); they are not stored.

**The path** is the latest milestone the body's timeline reached (`life_stage`: abiogenesis, photosynthesis, complex cells,
multicellularity, technological civilization). Life with a chemistry but no timeline (a moon, or a world outside the habitable
zone) is microbial and nothing more, as `habitability_world.native_life` already said. A body with no life chemistry, or a
timeline that has reached no milestone, scores 0 on all four.

**The surface** is the likelihood: `ceil(9 x PHI_bio)` for microbial life, `ceil(9 x PHI_cpx)` for complex life, never below 1
once the likelihood clears `tuning.HABITABILITY_LIFE_SCORE_MIN` (0.005), 0 under it.

| Score | Needs | Value |
|---|---|---|
| Microbial | any milestone | from PHI_bio |
| Flora | photosynthesis | from PHI_bio but at most 3 until multicellularity (mats and algae); from multicellularity the larger of that and PHI_cpx's step |
| Fauna | multicellularity | from PHI_cpx, never more than Flora (animals eat plants) |
| Sentient | technological civilization | Fauna's step times the share of 10 Myr (`LIFE_SCORE_SENTIENT_FULL_GYR`) the civilization has had since its milestone, at least 1 |

The thresholds (3 for simple flora, 10 Myr, the ceil step) are DESIGN values, in `tuning.py`, to be changed by a sweep or by Boss.
What the sweep of 150 systems (2,556 scored bodies) showed on 2026-10-23: 44 bodies carry a life chemistry, 38 of those have a
timeline (21 reached a technological civilization, 15 photosynthesis, 2 abiogenesis), 41 score above 0 Microbial; Sentient
equals Fauna on nearly every civilization world because the generator places a civilization's milestone more than 10 Myr before
"now" on almost all of them.

**Who has a path (Boss, 2026-10-23).** Every planet and moon, in any zone, that has a life chemistry gets an evolutionary
path ("moons always get an evolutionary path if appropriate for them, using the same classifications as planets"); the class's
own life chemistries and stage ceiling, and the habitability index's ceiling (`life.life_ceiling`), are what make it
appropriate. Before this only an ecosphere planet that was not a moon had one.

**Food comes from Flora and Fauna (Boss, 2026-10-23).** The Suitability food rung is no longer read from the stage of the
timeline: `life_scores.food_stage` turns the three scores into the old rungs. Flora is always good food (any Flora above 0 is
Green, 0.85), Flora with Fauna is better (Blue, 1.0), Fauna is all animal life and counts as food on its own (Green); microbes
alone stay Yellow (0.5) and no life Red (0.2); microbes that live only under an ice lid (or in a deep aquifer) are out of reach as before, but flora and fauna are food wherever the world scores them (Boss, 2026-10-24: a world with Flora 5, Fauna 5 and a civilization had no food because its ocean is under ice, GEN.765). The Suitability colour is
still held under the five core domains' ceiling.
