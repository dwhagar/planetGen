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

## 1. What the code has today

- **Per planet**: class and zone, radius, mass, density, gravity, orbit,
  rotation period, surface temperature (grey body times a per-class
  greenhouse factor), surface pressure, atmosphere density and scale
  height. `atmosphere` and `composition` are free text.
- **Missing**: partial pressures of any gas, ocean or land fraction, pH,
  salinity, magnetic field, radiation dose, mantle redox, ozone, humidity.
- **Per star**: type and class, mass, radius, temperature, luminosity,
  age, habitable zone, heliosphere radius. No XUV or flare activity.
- **Rogue planets** already have internal heat flux, a surface regime,
  ice-shell thickness and ocean depth (`rogueSurface.py`); nothing in a
  star system has these.
- **Nebulae**: a system's `surrounding_cloud` is set only when it is
  loaded from the database, so generation can't see it.
- **Life** is decided by class (`HABITABLE_PLANET_CLASSES`), and the life
  chemical by the star's spectral letter.

## 2. The scores

The documents give two structures, and one has to be chosen (a decision
for Boss):

- **PHI-4**: four domains, Pressure, Temperature, Chemistry and Radiation,
  each mapped through a logistic transfer function, combined as a weighted
  geometric mean times hard gates, with colour tiers (Blue, Green, Yellow,
  Red).
- **The Xenobiology doc**: three scores per planet:
  - PHI_bio (microbial life) = (L_solv · L_chem · L_ener · L_rad)^(1/4),
    from water activity and chaotropicity, C/N/H and phosphorus,
    available energy (light, chemistry, radiolysis) and dose;
  - PHI_cpx (complex life) = PHI_bio times a metazoan gate on CO, CO2 and
    O2 partial pressures;
  - Phi_tech (human operability) from pressure, wet-bulb temperature,
    local resources and shielded dose.

Default until Boss decides: PHI-4's domains and tiers for display, with
the three Xenobiology scores as the numbers behind them. Each planet also
gets a human equipment profile: shirtsleeve, breathing mask, mask with
scrubber, pressure suit, or full life support with radiation hardening.

## 3. Problems in the documents to settle first

1. PHI_bio's geometric mean doesn't reproduce the doc's own table (for
   example the dune world computes to 0.34 against a stated 0.05). A
   geometric mean forgives one near-zero factor; a product, a minimum
   gate or a weighted mean is needed.
2. Constants are undefined: k_w, X_ref, Phi_ref, sigma_P, alpha,
   E_specific, lambda, the scaling for non-water solvents, and how surface
   dose is computed.
3. L_rad's threshold (10 Gy/yr) disagrees with the table's Mars value.
4. Thresholds differ between docs: CO (0.005, 0.01 or 0.1 kPa), CO2
   (0.93, 2.0 or 0.5 to 5 kPa), the breathing-mask pressure band (6.3 to
   250 kPa or 50 to 250 kPa), and archetype pressure ranges.
5. Only M, K and G hosts are covered; the generator also makes F, A, B, O
   stars, white dwarfs, giants and binaries.
6. No document covers nebula effects on planets (section 6).

## 4. Inputs to build, in order

| # | Piece | Needs |
|---|---|---|
| 1 | Mantle redox (reduced, intermediate, oxidized) from mass and differentiation | mass, radius |
| 2 | Atmosphere species: partial pressures of O2, CO2, CO, N2, Ar, H2, H2O, CH4, H2S, SO2, with the free text generated from them | 1, pressure, class |
| 3 | Stellar activity: saturation phase, L_XUV/L_bol, flare and particle-event rates by mass and age | star mass, age, type |
| 4 | Magnetic field from mass, rotation and tidal locking | mass, rotation (the spin vector) |
| 5 | Surface dose from column mass (P0/g), field, cosmic rays (higher inside a compressed heliosphere) and flares, with an ozone-loss flag | 2, 3, 4, nebula containment |
| 6 | Hydrosphere: water fraction, ocean and land fraction, depth or ice shell | zone, temperature, class; reuse `rogueSurface` |
| 7 | Ocean chemistry: acid sulfate, neutral, soda, chloride brine or ice-sealed, with pH, water activity and phosphorus | 1, 2, 6 |
| 8 | The scores, tier and equipment profile | all of the above |

Values are stored in their own columns, not JSON, so they can be searched.

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
