# The wiring register: what is built or downloaded and not read yet (GEN.754)

Status: hand-kept, checked by `scripts/wiring_guard.py` and `src/tests/test_wiring_guard.py`. Every module, setup dataset, dataset reader, data file and stand-in registry row (GEN.725) that nothing in the app reaches is named here with the open item that will wire it, or with the words `by design` and the reason. The test fails when something unread is missing from this page, when a line names something that is read now (delete the line in the same change that wires it), and when a named item is not open in `docs/TODO.md`. No shim: an unread piece is wired, filed or deleted.

What counts as read: the app is every module outside `physics`, `generation`, `galaxy`, `datasets` and `tuning` (the CLI, API, web, admin, queue, population, backup and describe packages). A library module is read when an import chain from the app reaches it, a lazy import inside a function counts, and so does a module named by a string (`importlib.import_module("planetgen.x")`, the registry's `"planetgen.x:function"`). A dataset reader is read when a reachable module outside `datasets` uses it. A data file is read when a reachable module names it. A registry row is read when one of its `module.CONSTANT` names is used by a reachable module.

`python scripts/wiring_guard.py` prints what the guard finds; `--check` is the test.

| kind | name | item | reason |
|---|---|---|---|
| dataset | harris_clusters | GEN.297 | the Harris catalog of 157 globular clusters; the calibration of the cluster generator reads it |
| dataset | noaa_proton_events | GEN.805 | the 265 NOAA solar proton events; the solar particle event model fits its rate and slope to them |
| module | planetgen.physics.alpha_tracks | GEN.744 | the equivalent [Fe/H] of an alpha-enhanced star; the track lookup of a high-alpha star will read it |
| module | planetgen.physics.age_metallicity | GEN.744 | the disc's age-metallicity relation by branch; the star's [Fe/H] draw will follow its age |
| row | 6.34:AL1 | GEN.744 | `alpha_tracks.ALPHA_SHARE`, not read until the track lookup is wired |
| row | 6.34:AR1 | GEN.744 | `age_metallicity.THIN_BIRTH_FE_H`, not read until the [Fe/H] draw follows the age |
| row | 6.34:AR2 | GEN.744 | `age_metallicity.THICK_ANCHOR_FE_H`, not read until the [Fe/H] draw follows the age |
| row | 6.34:AR3 | GEN.744 | `age_metallicity.SIGMA`, not read until the [Fe/H] draw follows the age |
| module | planetgen.physics.line_opacity | by design | the band opacities of the radiative-convective model, which the greenhouse table is fitted to and the checks hold it to; the generator reads the table, never the model |
| module | planetgen.physics.radiative_convective | by design | the model the greenhouse table (`greenhouse_table.py`) is fitted to by `scripts/fit_greenhouse_table.py`, too slow for the fixed point; run by the fit and the tests, not by generation |
| module | planetgen.physics.radiative_kernel | by design | the numba sweep of `radiative_convective`; the web app must not load numba, so only the fit and the tests reach it |
| module | planetgen.physics.clearance | by design | the independent overlap verifier (TEST.477); the tests and the redraw audit read it, the generator never does, so the proof cannot share a bug with the code it checks |
| module | planetgen.physics.photolysis | GEN.982 | the photolysis lifetime table for the exosphere and the coma losses; no caller yet |
| data | photolysis_lifetimes.json | GEN.982 | read by `physics/photolysis.py` only |
| module | planetgen.physics.mixture_trace | GEN.524 | the trace gases of an air drawn after the life stage; no caller yet |
| module | planetgen.physics.design_flags | GEN.524 | gates the trace rows `mixture_trace` draws |
| data | design_flags.json | GEN.524 | read by `physics/design_flags.py` only |
| module | planetgen.physics.equilibrium_edges | GEN.994 | the edges of the hot-Jupiter chemistry regimes, sourced and wired with the hot-giant thermochemistry |
| module | planetgen.physics.ground_albedo | GEN.684 | the colour of the ground by band, wired by the ground albedo item |
| module | planetgen.physics.giant_dynamo | GEN.809 | the dynamo region and field of a giant from its interior layers (lock-free); the wiring step has `magnetism.generate_field` take a giant's moment from it in place of the three anchors, once the structure of GEN.733 has placed an ice giant's layers (GEN.967 wires `dynamo` for the rest) |
| module | planetgen.physics.outside_hours | GEN.768 | the hours a person can be outside in a day, stored per planet by the item |
| module | planetgen.physics.surface_liquids | GEN.1035 | the liquids a surface can hold; the stored surface_liquid, its coverage and the ocean tide read it once the wiring turn (train 9) lands |
| row | 6.69:LQ6 | GEN.1035 | `surface_liquids.SOLUBILITY_FACTOR`, read only by the ideal sea until the wiring turn (train 9) |
| data | surface_liquids.json | GEN.1035 | read by `physics/surface_liquids.py` only |
| module | planetgen.physics.radiation_belts | GEN.983 | trapped radiation belts of a magnetised body, read when the stored magnetosphere reaches dose |
| module | planetgen.galaxy.influence | GEN.342 | the influence set of an object, read by the acceleration kernel of the orbital update |
| module | planetgen.physics.closest_approach | GEN.197 | the pair closest-approach library, read by the orbital event writer |
| module | planetgen.physics.collision_outcome | GEN.197 | what a collision does, read by the orbital event writer |
| module | planetgen.physics.body_positions | MAP.253 | where every body of a system is at a time, read by the scene and the maps once they read the parent chain |
| module | planetgen.physics.light_travel | VIEW.1 | light-travel positions for the view from a planet |
| module | planetgen.physics.radius_valley | by design | the simulated Kepler survey of the radius valley, a check run by the tuning tests; the valley itself comes from the envelope-loss rows of the generator (GEN.722 retired) |
| module | planetgen.physics.iron_branch_tuning | by design | the simulated sample the iron branch table is tuned to; run by the tuning tests, not by generation |
| module | planetgen.generation.redraw_kinds | by design | the names and one-line meanings of the redraw counters (TEST.476); `draw.redrawn` is called by name from the generator and the census reads this registry, so no generation code imports it |
| module | planetgen.generation.class_sample | by design | a sample that holds every class, a helper of the sheet audit tests |
| module | planetgen.generation.plausibility | by design | the anomaly finder engine; run by the plausibility tests and their command line wrapper in `src/tests` |
| module | planetgen.generation.phenomena_plausibility | by design | the exotic-phenomena counterpart of the anomaly finder, run by its tests |
| row | 6.7:AM18 | GEN.524 | `mixture_trace.ABUNDANCES`, not drawn until the trace gases are wired |
| row | 6.8:OH8 | GEN.768 | `protection_limits.GARMENTS` insulation, read by `outside_hours` |
| row | 6.8:OH9 | GEN.768 | `protection_limits.GARMENTS` cooling and heating power, read by `outside_hours` |
| row | 6.8:OH11 | GEN.768 | `outside_hours.LOCKED_BELT_LIGHT_HOURS`, read by `outside_hours.for_body`, which nothing calls until the stored columns are wired |
| row | 6.10:RV2 | by design | `radius_valley.CORE_MASS_SIGMA`, a width of the simulated survey, read by the tuning tests |
| row | 6.10:RV3 | by design | `radius_valley.ENVELOPE_SIGMA`, a width of the simulated survey, read by the tuning tests |
| row | 6.15:RB1 | GEN.983 | `radiation_belts.WIND_EDGE_FRACTION`, read when the belts are |
| row | 6.15:RB2 | GEN.983 | `radiation_belts.ROTATION_EDGE_FACTOR`, read when the belts are |
| row | 6.15:RB3 | GEN.983 | `radiation_belts.WIND_PEAK_FRACTION`, read when the belts are |
| row | 6.15:RB4 | GEN.983 | `radiation_belts.INNER_PEAK_FRACTION`, read when the belts are |
| row | 6.15:RB5 | GEN.983 | `radiation_belts.ATMOSPHERE_EDGE_SCALE_HEIGHTS`, read when the belts are |
| row | 6.15:RB6 | GEN.983 | `radiation_belts.EARTH_OUTER_PEAK_GY_YR`, read when the belts are |
| row | 6.15:RB7 | GEN.983 | `radiation_belts.OUTER_RISE_EXPONENT`, read when the belts are |
| module | planetgen.physics.tidal_rheology | GEN.955 | the rheology of a rocky or icy body under a tide (stage 1a, lock-free); stage 1b wires it into tidal heating, interior heat and circularisation |
| module | planetgen.physics.spin_torques | GEN.957 | the star, moon and thermal-tide torques on a planet's spin (stage 2a, lock-free); stage 2b wires the sum into planet spin |
| module | planetgen.physics.moon_orbits | GEN.959 | moon orbit drift, spiral-in and the Roche fall (stage 3a, lock-free); stage 3b wires the history into the planet |
| module | planetgen.physics.impact_spin | GEN.979 | Kokubo and Genda equation 4 with the gamma law of the last impacts (GEN.978, lock-free); GEN.979 wires it into the rocky planet's birth spin |
| row | 6.23:TR1 | GEN.955 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.23:TR2 | GEN.955 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.23:TR4 | GEN.955 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.23:TR5 | GEN.955 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.23:TR8 | GEN.955 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.23:TR10 | GEN.955 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.24:TQ1 | GEN.957 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.24:TQ2 | GEN.957 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.24:TQ3 | GEN.957 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.24:TQ4 | GEN.957 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.25:TM1 | GEN.959 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.25:TM2 | GEN.959 | a constant of the tides chain's lock-free module, read when its wiring step lands |
| row | 6.26:GS2 | GEN.979 | a constant of the birth-spin module, read when its wiring step lands |
| row | 6.26:GS3 | GEN.979 | a constant of the birth-spin module, read when its wiring step lands |
| row | 6.26:GS4 | GEN.979 | a constant of the birth-spin module, read when its wiring step lands |
| row | 6.26:GS5 | GEN.979 | a constant of the birth-spin module, read when its wiring step lands |
| row | 6.26:GS6 | GEN.979 | a constant of the birth-spin module, read when its wiring step lands |
| module | planetgen.physics.planet_occurrence | by design | the reader of the observed occurrence bands; the architecture and the census read it (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
| module | planetgen.physics.kepler_multiplicity | by design | the Kepler transit multiplicity harness, read by the tuning tests until the wiring (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
| data | planet_occurrence.json | by design | read by `physics/planet_occurrence.py` only (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
| module | planetgen.physics.substellar_ceiling | GEN.494 | the Baraffe 2015 hydrogen-burning limit temperature by age (lock-free); the wiring step replaces the 2800 K clip of `brown_dwarf` and the hot edge of the atmosphere row g1 with it, riding the Baraffe brown dwarf tracks of GEN.948 |
| row | 6.37:SB1 | GEN.494 | `substellar_ceiling.LIMIT_MASS_SOLAR`, read when the wiring step lands |
| row | 6.37:SB2 | GEN.494 | `substellar_ceiling.HOTTEST_K`, read when the wiring step lands |
| module | planetgen.physics.birth_radioisotopes | GEN.476 | aluminium-26 and iron-60 at a star's birth from the gamma-ray masses, the region's depletion time and the granularity of the sources (lock-free); the wiring step has `radiogenic.draw_content` take the 26Al ratio from it, and row o6 its enriched share |
| row | 6.38:BR1 | GEN.476 | `birth_radioisotopes.MEAN_LIFE_YR`, read when the wiring step lands |
| row | 6.38:BR2 | GEN.476 | `birth_radioisotopes.GALAXY_AL26_MASS_SOLAR, FE60_AL26_FLUX_RATIO`, read when the wiring step lands |
| row | 6.38:BR3 | GEN.476 | `birth_radioisotopes.GALAXY_STAR_FORMATION_SOLAR_MASSES_PER_YR`, read when the wiring step lands |
| row | 6.38:BR4 | GEN.476 | `birth_radioisotopes.RECURRENCE_RANGE_YR`, read when the wiring step lands |
| row | 6.38:BR5 | GEN.476 | `birth_radioisotopes.EVENT_YIELD_SIGMA_LN`, read when the wiring step lands |
| module | planetgen.physics.saturation_time | GEN.493 | the X-ray saturation time by stellar mass and initial rotation percentile from Johnstone et al. 2021 (lock-free); the wiring step has row p6 and the activity table's `t_sat` read it with the star's own mass |
| row | 6.39:RS1 | GEN.493 | `saturation_time.ANCHORS_MYR`, read when the wiring step lands |
| module | planetgen.physics.taxon_kernel | GEN.473 | the taxa of a belt from the disk's ice lines and the early heating by aluminium-26 (lock-free); the wiring step replaces `small_body_composition.kernel` with it and hands it the belt's temperature and the star's 26Al ratio |
| row | 6.41:TK1 | GEN.473 | `taxon_kernel.WATER_LINE_K, AMMONIA_LINE_K`, read when the wiring step lands |
| row | 6.41:TK2 | GEN.473 | `taxon_kernel.PRIMITIVE_BELOW_K, WIDTH_LN_T`, read when the wiring step lands |
| row | 6.41:TK3 | GEN.473 | `taxon_kernel.IMPLANT_FRACTION`, read when the wiring step lands |
| row | 6.41:TK4 | GEN.473 | `taxon_kernel.MELT_ENERGY_J_PER_KG, ALUMINIUM_MASS_FRACTION`, read when the wiring step lands |
| row | 6.41:TK5 | GEN.473 | `taxon_kernel.ACCRETION_SPREAD_MYR`, read when the wiring step lands |
| row | 6.41:TK6 | GEN.473 | `taxon_kernel.DRY_SPLIT, IGNEOUS_SPLIT, HYDRATED_SPLIT`, read when the wiring step lands |
| module | planetgen.physics.belt_mass | GEN.477 | the mass left in a planetesimal belt after collisional grinding (Wyatt et al. 2007; lock-free); the wiring step gives `AsteroidBelt` its mass and its density word from it |
| row | 6.43:BM1 | GEN.477 | `belt_mass.GRAIN_DENSITY_KG_M3, STRENGTH_FLOOR_J_PER_KG`, read when the wiring step lands |
| row | 6.43:BM2 | GEN.477 | `belt_mass.DENSE_FACTOR`, read when the wiring step lands |
| module | planetgen.physics.meteorite_groups | GEN.475 | the chondrite groups as regions of a composition space with a nearest-group label (lock-free); the wiring step labels belts and fields from the derived bulk and replaces `taxon_kernel` TK6's fixed within-family splits |
| row | 6.44:MG1 | GEN.475 | `meteorite_groups.CHONDRITE_GROUPS`, read when the wiring step lands |
| row | 6.44:MG2 | GEN.475 | `meteorite_groups.AXIS_SCALES, NEAREST_LIMIT`, read when the wiring step lands |
| module | planetgen.physics.rock_types | GEN.474 | the CIPW norm, the staged melt of a rock and the ground type of its oxides (lock-free); the wiring step has the surface draw take the crust from the stored mantle through it, replacing the fixed assemblages of `surface_composition.GROUND_TYPES` |
| row | 6.45:RT1 | GEN.474 | `rock_types.MELT_STAGES`, read when the wiring step lands |
| row | 6.45:RT2 | GEN.474 | `rock_types.ANORTHOSITE_PLAGIOCLASE, ULTRABASIC_SILICA, KOMATIITE_MAGNESIA, BASIC_SILICA, ACID_SILICA`, read when the wiring step lands |
| row | 6.45:RT3 | GEN.474 | `rock_types.INCOMPATIBLE_D`, read when the wiring step lands |
| module | planetgen.physics.comet_ice_budget | GEN.595 | the CO2/H2O and CO/H2O of a comet's ice from the disc zone it formed in (lock-free); the wiring step has `comet_population.draw_volatiles` take its CO and CO2 medians from the zone and the interstellar bracket from the disc's carbon and oxygen |
| row | 6.46:CI1 | GEN.595 | `comet_ice_budget.ICE_LINES_K, DISK_T0_K, DISK_INDEX`, read when the wiring step lands |
| row | 6.46:CI2 | GEN.595 | `comet_ice_budget.LUMINOSITY_EXPONENT`, read when the wiring step lands |
| row | 6.46:CI3 | GEN.595 | `comet_ice_budget.NEBULA_PER_1E4_H`, read when the wiring step lands |
| row | 6.46:CI4 | GEN.595 | `comet_ice_budget.TRAPPED_FRACTION`, read when the wiring step lands |
| module | planetgen.physics.dry_carbon | GEN.1005 | the carbon a dry body keeps in its mantle from its building blocks and the core's share; the flip step has `mixture_rocky.INVENTORY['co2_dry_wt_percent']` read its bracket |
| row | 6.58:DC1 | GEN.1005 | `dry_carbon.NC_UNDIFFERENTIATED_C_PPM, NC_DIFFERENTIATED_C_PPM`, read when the flip step lands |
| row | 6.58:DC2 | GEN.1005 | `dry_carbon.DIFFERENTIATED_SHARE`, read when the flip step lands |
| row | 6.58:DC3 | GEN.1005 | `dry_carbon.MANTLE_SHARE, CO2_PER_C`, read when the flip step lands |
| data | potassium_depletion.json | GEN.492 | the K/U of six bodies (Dauphas et al. 2022 Table 1) with CI chondrites; `radiogenic.draw_content` will take the potassium depletion and CI K/U from it |
| module | planetgen.physics.potassium_depletion | GEN.492 | the potassium depletion drawn by body class; replaces the log-uniform `tuning.RADIOGENIC_K_DEPLETION_RANGE` once `radiogenic` calls it with the body's mass |
| row | 6.61:KD3 | GEN.492 | `potassium_depletion.CLASS_MASS_EARTH`, not read until `radiogenic` calls the module |
| row | 6.61:KD4 | GEN.492 | `potassium_depletion.MAX_DEPLETION`, not read until `radiogenic` calls the module |
| data | co2_far_wing.json | GEN.1019 | the carbon dioxide far-wing rules (Perrin and Hartmann chi, widths, reach); `line_opacity` will read them when the wing replaces IR3 |
| module | planetgen.physics.co2_far_wing | GEN.1019 | the far wing of carbon dioxide; replaces the flat `CO2_FAR_WING_M2_KG` of IR3 once `line_opacity` calls it |
| module | planetgen.physics.clearance | GEN.1016 | the independent overlap checker (audit of a stored system, TEST.477); read by tests/test_clearance.py now and by the stability gate of GEN.1016 once it lands |
| data | sulfuric_cloud.json | GEN.1019 | the sulfuric acid cloud rules (base temperature, published depth and radius brackets, cover); `aerosols.draw_cloud` will read them in the wiring step |
| module | planetgen.physics.sulfuric_cloud | GEN.1019 | the acid cloud of an air with sulfur; `aerosols.draw_cloud` will add it to the condensate cloud once wired |
| data | acid_cloud_optics.json | by design | the acid droplets' infrared optics per band (made by `scripts/build_acid_cloud_optics.py`); `radiative_convective` will read it through `cloud_infrared` in the wiring step (read by the radiative-convective model, which the generator reads through its stored table, like line_opacity; GEN.850 built) |
| module | planetgen.physics.cloud_infrared | by design | the acid cloud's infrared depth in the layers of a column; called by the column solver once the cloud is wired (read by the radiative-convective model, which the generator reads through its stored table, like line_opacity; GEN.850 built) |
| module | planetgen.physics.mie | by design | Mie spheres, used by the optics builder script only; the generator reads the table it makes, so this stays unread by the app (delete when the builder moves out of the repo) |
| module | planetgen.physics.pair_history | GEN.1012 | the history rule of a pair of members; `multistar.classify_pair` will read it once the member record (GEN.1009) provides each member's largest mass and radius |
| row | 6.63:CF3 | by design | `cloud_infrared.TOP_OVER_BASE_PRESSURE`, read when the column solver calls the cloud (read by the radiative-convective model, which the generator reads through its stored table, like line_opacity; GEN.850 built) |
| module | planetgen.physics.clearance_ledger | GEN.1014 | the spans where a planet, moon or belt may sit; replaces `_orbit_floor_au`, the per-star `a_crit_au` reads and the spacing and `min_orbit_distance` pushes of `space_orbits` once the system build draws from it (train 10, after GEN.1009 and GEN.1011) |
| data | quadruple_stability_3p1.json | GEN.1016 | the 3+1 quadruple stability rules (validity ranges, the stability definition); `multistar` limits will read them through `quadruple_stability` in the wiring step |
| data | quadruple_stability_3p1_weights.json | GEN.1016 | the 8,301 numbers of the authors' network (made by `scripts/build_quadruple_stability.py`); read by `quadruple_stability` |
| data | quadruple_stability_2p2.json | GEN.1016 | the 2+2 quadruple stability rules (validity ranges, the stability definition); read by `quadruple_stability` |
| data | quadruple_stability_2p2_weights.json | GEN.1016 | the 8,301 numbers of the authors' 2+2 network (made by `scripts/build_quadruple_stability.py`); read by `quadruple_stability` |
| module | planetgen.physics.quadruple_stability | GEN.1016 | the stability verdict of a 3+1 quadruple and (GEN.1031) of a 2+2 quadruple; the stability margin and quadruple fit (P10) will call it |
| module | planetgen.generation.member_record | GEN.1009 | the frozen member record (one evolution of a member); the tree, the build and the stored tree will read it once GEN.1009's wiring merges on its train |
| module | planetgen.physics.ocean_tides | GEN.1034 | the open-ocean equilibrium tide and its spring and neap ranges; the fill pass will call it with the stored spin and moons and write `tide_range_spring_m` and `tide_range_neap_m` (train 7)|
| row | 6.68:OT1 | GEN.1034 | a constant of the ocean tide module, read when its wiring step lands |
| row | 6.68:OT2 | GEN.1034 | a constant of the ocean tide module, read when its wiring step lands |
| row | 6.68:OT6 | GEN.1034 | a constant of the ocean tide module, read when its wiring step lands |
| module | planetgen.physics.small_body_rings | GEN.1057 | the resonance ring of a small body (R4) and the gate against rings inside 0.1 AU (R9), pure functions (lock-free); the wiring step calls them from the ring pass |
| row | 6.71:RR1 | GEN.1057 | `small_body_rings.RESONANCE_RING_RATE, RESONANCE_ORDER`, read when the wiring step lands |
| row | 6.71:RR2 | GEN.1057 | `small_body_rings.INNER_RING_LIMIT_AU`, read when the wiring step lands |
| row | 6.71:RR3 | GEN.1057 | `small_body_rings.HILL_FRACTION_MAX`, read when the wiring step lands |
| row | 6.71:RR4 | GEN.1057 | `small_body_rings.RESONANCE_RING_MASS_KG, RESONANCE_WIDTH_FRACTION`, read when the wiring step lands |
| module | planetgen.physics.ring_origin | GEN.1057 | the ring as a result of a fallen moon or a comet inside the Roche limit; the conditional ring draw (GEN.1057) will call it and replace the `rings.draw_ring` chances |
| row | 6.70:RO1 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO2 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO3 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO4 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO5 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO6 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO7 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO8 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| module | planetgen.physics.oddities | GEN.1047 | the oddity table and its rolls (candidates, probability, roll, limit_failures); the planet and moon builds call it once GEN.1047 wires the table, with the effects the modules of GEN.1040 and its siblings define |
| data | oddities.json | GEN.1047 | the oddity table (empty until the modules add their rows); read by physics/oddities.py, which nothing calls yet |
| row | 6.72:OD1 | GEN.1047 | `tuning.WEIRDNESS_DIAL`, read by `oddities.probability` once the table is wired; the dial also goes in the settings page |
| row | 6.72:OD2 | GEN.1047 | `tuning.WEIRDNESS_RATE_SCALE`, read by `oddities.probability` once the table is wired |
| module | planetgen.physics.disintegrating_planet | GEN.1050 | the wind, lifetime and gate of a disintegrating rocky planet (W14); the oddity row reads it when its wiring step lands |
| row | 6.73:DP1 | GEN.1050 | the olivine vapour constants, not read until the oddity row is wired |
| row | 6.73:DP2 | GEN.1050 | the vapour mean mass and index, not read until the oddity row is wired |
| row | 6.73:DP3 | GEN.1050 | the full model's factor and the solid angle, not read until the oddity row is wired |
| row | 6.73:DP4 | GEN.1050 | the duty cycle, not read until the oddity row is wired |
| row | 6.73:DP5 | GEN.1050 | the gate, not read until the oddity row is wired |
| module | planetgen.physics.atmosphere_retention | GEN.1047 | the shoreline ratio, verdict, pressure cap and Jeans gate of the design note; the pressure stage reads the logistic of rows d until the envelope-survival rule (GEN.963) reads the module and the pressure stage takes its cap from it |
| row | 6.74:SH1 | GEN.1047 | `atmosphere_retention.SHORELINE_K`, not read until the pressure stage takes its cap from the module |
| row | 6.74:SH2 | GEN.1047 | `atmosphere_retention.VERDICT_EDGES`, not read until the pressure stage takes its cap from the module |
| row | 6.74:SH3 | GEN.1047 | `atmosphere_retention.PRESSURE_CAPS_PA`, not read until the pressure stage takes its cap from the module |
| row | 6.74:SH4 | GEN.1047 | `atmosphere_retention.JEANS_KEEP`, not read until the pressure stage takes its cap from the module |
| module | planetgen.physics.nightside_collapse | GEN.1051 | whether a locked planet's air freezes onto its night side and how much stays; the oddity wiring (GEN.1051) will call it after the climate solver |
| row | 6.75:NC1 | GEN.1051 | a constant of the nightside collapse module, read when its wiring step lands |
| row | 6.75:NC3 | GEN.1051 | a constant of the nightside collapse module, read when its wiring step lands |
| row | 6.75:NC6 | GEN.1051 | a constant of the nightside collapse module, read when its wiring step lands |
| module | planetgen.physics.ring_origin | GEN.1057 | the ring as a result of a fallen moon or a comet inside the Roche limit; the conditional ring draw (GEN.1057) will call it and replace the `rings.draw_ring` chances |
| row | 6.70:RO1 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| module | planetgen.physics.venus_spin | GEN.1048 | the spin equilibria of a thick-air planet (thermal tide against the solid tide, Cunha, Correia and Laskar 2015); the weird-worlds wiring (W2) will give a not-locked thick-air planet its slow or retrograde spin |
| row | 6.90:VS1 | GEN.1048 | a constant of the Venus-type spin module, read when its wiring step lands |
| row | 6.90:VS2 | GEN.1048 | a constant of the Venus-type spin module, read when its wiring step lands |
| row | 6.90:VS4 | GEN.1048 | a constant of the Venus-type spin module, read when its wiring step lands |
| row | 6.90:VS5 | GEN.1048 | a constant of the Venus-type spin module, read when its wiring step lands |
| row | 6.90:VS6 | GEN.1048 | a constant of the Venus-type spin module, read when its wiring step lands |
| row | 6.90:VS7 | GEN.1048 | a constant of the Venus-type spin module, read when its wiring step lands |
| row | 6.90:VS8 | GEN.1048 | a constant of the Venus-type spin module, read when its wiring step lands |
| module | planetgen.physics.seasonal_swing | GEN.1048 | the temperature through the year at a latitude from the tilt and the orbit; the weird-worlds wiring (W1) will add the seasonal swing to the regional temperatures |
| row | 6.91:SS1 | GEN.1048 | a constant of the seasonal swing module, read when its wiring step lands |
| row | 6.91:SS2 | GEN.1048 | a constant of the seasonal swing module, read when its wiring step lands |
| row | 6.91:SS3 | GEN.1048 | a constant of the seasonal swing module, read when its wiring step lands |
| module | planetgen.physics.resonant_chains | GEN.1049 | the resonant chains of a system's periods and the eccentricity they force; the weird-worlds wiring (W5) will give a chain planet its forced orbit and tidal heating |
| row | 6.92:RC1 | GEN.1049 | a constant of the resonant chains module, read when its wiring step lands |
| row | 6.92:RC2 | GEN.1049 | a constant of the resonant chains module, read when its wiring step lands |
| row | 6.92:RC3 | GEN.1049 | a constant of the resonant chains module, read when its wiring step lands |
| row | 6.92:RC4 | GEN.1049 | a constant of the resonant chains module, read when its wiring step lands |
| module | planetgen.physics.belt_density | UX.240 | the objects per square AU above a diameter, the mass per square AU and the size-table bins of a belt; the belt page data builder reads it |
| row | 6.99:BD1 | UX.240 | a constant of the belt density module, read when its wiring step lands |
| row | 6.99:BD2 | UX.240 | a constant of the belt density module, read when its wiring step lands |
| row | 6.99:BD3 | UX.240 | a constant of the belt density module, read when its wiring step lands |
| row | 6.99:BD4 | UX.240 | a constant of the belt density module, read when its wiring step lands |
| row | 6.99:BD5 | UX.240 | a constant of the belt density module, read when its wiring step lands |
| module | planetgen.physics.band_winds | GEN.1110 | the mean surface wind by day and by night in each latitude band; the seven-band climate (GEN.1108) will feed it band temperatures, and GEN.1110's storage step reads it |
| row | 6.98:BW1 | GEN.1110 | a constant of the band winds module, read when its wiring step lands |
| row | 6.98:BW2 | GEN.1110 | a constant of the band winds module, read when its wiring step lands |
| row | 6.98:BW3 | GEN.1110 | a constant of the band winds module, read when its wiring step lands |
| row | 6.98:BW4 | GEN.1110 | a constant of the band winds module, read when its wiring step lands |
| row | 6.98:BW5 | GEN.1110 | a constant of the band winds module, read when its wiring step lands |
| row | 6.98:BW6 | GEN.1110 | a constant of the band winds module, read when its wiring step lands |
| module | planetgen.physics.rogue_drift | GEN.1111 | the closed-form drift of an ejected planet from its star; the virtual-donor placement of rogue planets reads it |
| row | 6.89:GD1 | GEN.1111 | a constant of the rogue drift module, read when its wiring step lands |
| row | 6.89:GD2 | GEN.1111 | a constant of the rogue drift module, read when its wiring step lands |
| row | 6.89:GD3 | GEN.1111 | a constant of the rogue drift module, read when its wiring step lands |
| row | 6.89:GD4 | GEN.1111 | a constant of the rogue drift module, read when its wiring step lands |
| module | planetgen.physics.rogue_age | GEN.1111 | a rogue's age, ejection delay and time since ejection; the rogue rows read them for the heat law, the drift and the history sentence |
| row | 6.88:RG1 | GEN.1111 | a constant of the rogue age module, read when its wiring step lands |
| row | 6.88:RG2 | GEN.1111 | a constant of the rogue age module, read when its wiring step lands |
| module | planetgen.physics.rogue_kick | GEN.1111 | the ejection speed of a rogue and its velocity through the Galaxy; the rogue rows read them |
| row | 6.87:RK1 | GEN.1111 | a constant of the rogue kick module, read when its wiring step lands |
| row | 6.87:RK2 | GEN.1111 | a constant of the rogue kick module, read when its wiring step lands |
| module | planetgen.physics.rogue_mass_function | GEN.1111 | the masses of the bodies a system ejects and their shares; the rogue rows draw their masses from it |
| row | 6.110:RM1 | GEN.1111 | a constant of the rogue mass function module, read when its wiring step lands |
| row | 6.110:RM2 | GEN.1111 | a constant of the rogue mass function module, read when its wiring step lands |
| module | planetgen.physics.rogue_volatiles | GEN.1111 | the envelope, water, giant class and frozen air of an ejected planet; the rogue rows read them |
| row | 6.111:EV1 | GEN.1111 | a constant of the rogue volatiles module, read when its wiring step lands |
| row | 6.111:EV2 | GEN.1111 | a constant of the rogue volatiles module, read when its wiring step lands |
| module | planetgen.physics.rogue_moons | GEN.1111 | whether a moon leaves with its planet and the tide it raises there; the rogue rows read them |
| row | 6.112:MO1 | GEN.1111 | a constant of the rogue moons module, read when its wiring step lands |
| row | 6.112:MO3 | GEN.1111 | a constant of the rogue moons module, read when its wiring step lands |
| module | planetgen.physics.cluster_unbinding | GEN.1111 | the share of planets a birth cluster has unbound; the rogue rows add it to the ejections |
| row | 6.70:RO2 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO3 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO4 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO5 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO6 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO7 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| row | 6.70:RO8 | GEN.1057 | a constant of the ring origin module, read when its wiring step lands |
| module | planetgen.physics.super_puff | GEN.1048 | the super-puff radius, fade and gate (W15); the bulk step reads it when the wiring step lands |
| row | 6.76:PF1 | GEN.1048 | a constant of the super-puff module, read when its wiring step lands |
| row | 6.76:PF2 | GEN.1048 | a constant of the super-puff module, read when its wiring step lands |
| row | 6.76:PF3 | GEN.1048 | a constant of the super-puff module, read when its wiring step lands |
| row | 6.76:PF4 | GEN.1048 | a constant of the super-puff module, read when its wiring step lands |
| module | planetgen.physics.helium_fractionation | GEN.1047 | the helium share of an envelope that loses hydrogen first (W16); the oddity row reads it when the oddity wiring lands |
| row | 6.77:HF1 | GEN.1047 | a constant of the helium fractionation module, read when its wiring step lands |
| row | 6.77:HF2 | GEN.1047 | a constant of the helium fractionation module, read when its wiring step lands |
| row | 6.77:HF3 | GEN.1047 | a constant of the helium fractionation module, read when its wiring step lands |
| module | planetgen.physics.giant_tilt | GEN.1058 | the tilt a giant merger gives an ice giant and the ring and moon plane that follow the equator (G1); the spin draw reads it when the wiring lands |
| row | 6.78:GT1 | GEN.1058 | a constant of the giant tilt module, read when its wiring step lands |
| row | 6.78:GT2 | GEN.1058 | a constant of the giant tilt module, read when its wiring step lands |
| row | 6.78:GT3 | GEN.1058 | a constant of the giant tilt module, read when its wiring step lands |
| module | planetgen.physics.mantle_stripping | GEN.1048 | the iron-rich core of a giant impact (Marcus and others 2010 scaling); the weirdness wiring (W10) gives a rocky planet that had a recorded merger a core fraction from the impacts the architecture recorded |
| row | 6.93:MT4 | GEN.1048 | `mantle_stripping.Q_STAR_SCALE` of the mantle stripping module, read when its wiring step lands |
| row | 6.93:MT6 | GEN.1048 | `mantle_stripping.MAX_CORE_MASS_FRACTION` of the mantle stripping module, read when its wiring step lands |
| row | 6.93:MT7 | GEN.1048 | `mantle_stripping.IMPACT_ANGLE_DEG` of the mantle stripping module, read when its wiring step lands |
| module | planetgen.physics.hycean_candidate | GEN.1047 | the Hycean box as a flag on the classes T, I and V; the oddity line names the planets it selects when the oddity wiring lands |
| row | 6.79:HC1 | GEN.1047 | the four ranges of the Hycean box, read when its wiring step lands |
| module | planetgen.physics.irregular_zone | GEN.1076 | the stable limits, Kozai band and inner edge of the irregular zone; the irregular orbit draw (GEN.1078) and the captured-moon wiring will call it |
| row | 6.94:IZ1 | GEN.1076 | a constant of the irregular zone module, read when its wiring step lands |
| row | 6.94:IZ2 | GEN.1076 | a constant of the irregular zone module, read when its wiring step lands |
| row | 6.94:IZ3 | GEN.1076 | a constant of the irregular zone module, read when its wiring step lands |
| module | planetgen.physics.rogue_collision | GEN.1103 | the rate, disruption speed, largest remnant, fragment speeds and shock of a rogue-rogue collision; the field columns and host rule (GEN.1103, GEN.1104) call it |
| module | planetgen.physics.irregular_orbits | GEN.1083 | the a / r_H, eccentricity, inclination and family draw of irregular moons; the captured-moon wiring will call it |
| row | 6.96:IO1 | GEN.1083 | a constant of the irregular orbits module, read when its wiring step lands |
| row | 6.96:IO2 | GEN.1083 | a constant of the irregular orbits module, read when its wiring step lands |
| row | 6.96:IO3 | GEN.1083 | a constant of the irregular orbits module, read when its wiring step lands |
| row | 6.96:IO4 | GEN.1083 | a constant of the irregular orbits module, read when its wiring step lands |
| module | planetgen.physics.scatter | GEN.1092 | the shares of an outer disc's mass in each scattered population, the rule of the scatter phase; the wiring step feeds it the architecture's giants |
| row | 6.83:SX1 | GEN.1092 | a constant of the scatter module, read when its wiring step lands |
| row | 6.83:SX2 | GEN.1092 | a constant of the scatter module, read when its wiring step lands |
| row | 6.83:SX3 | GEN.1092 | a constant of the scatter module, read when its wiring step lands |
| row | 6.83:SX4 | GEN.1092 | a constant of the scatter module, read when its wiring step lands |
| row | 6.83:SX5 | GEN.1092 | a constant of the scatter module, read when its wiring step lands |
| module | planetgen.physics.belt_ledger | GEN.1069 | the stored bodies of a belt by the mass cut, the size laws and the largest body settled with the collisional mass (lock-free, GEN.1061); the belt store and the belt generator call it, the outer disc (GEN.1071) and the scatter step (GEN.1092) after them |
| row | 6.95:BR17 | GEN.1069 | `belt_ledger.L1_MASS_KG`, read when the wiring step lands |
| row | 6.95:BR18 | GEN.1069 | `belt_ledger.LEDGER_COVERAGE, LEDGER_MAX_BODIES`, read when the wiring step lands |
| row | 6.95:BR19 | GEN.1069 | `belt_ledger.ICY_LAW`, read when the wiring step lands |
| row | 6.95:BR20 | GEN.1069 | `belt_ledger.ROCKY_LAW`, read when the wiring step lands |
| row | 6.95:BR21 | GEN.1069 | `belt_ledger.LEDGER_MAX_SHARE, SOLVE_TOLERANCE, SOLVE_MIN_KM`, read when the wiring step lands |
| module | planetgen.physics.body_status | GEN.1069 | the status word (planet, dwarf planet, minor body, moon) and its numbers for each stored body (lock-free, GEN.1063); the belt store and the member pages call it |
| row | 6.95:BR22 | GEN.1069 | `body_status.ZONE_HILL_RADII`, read when the wiring step lands |
| row | 6.95:BR23 | GEN.1069 | `body_status.LAMBDA_K`, read when the wiring step lands |
| row | 6.95:BR24 | GEN.1069 | `body_status.ROUND_ROCK_SHIFT_DECADES, ROUND_THRESHOLD`, read when the wiring step lands |
| module | planetgen.physics.belt_origin | GEN.1092 | the cause of each belt (gap, tail, snow line, outer disc, destroyed planet, rubble), what the giants cleared and a moving giant swept, and the sentence for the belt page; GEN.1068 gives it the final chain and GEN.1092 stores its result |
| data | belt_origins.json | GEN.1092 | the causes of a belt and the sentence each prints; read by physics/belt_origin.py, which nothing calls yet |
| row | 6.97:BG1 | GEN.1092 | `belt_origin.GAP_SPACING_K`, read when the belt origin wiring lands |
| row | 6.97:BG2 | GEN.1092 | `belt_origin.MIN_BELT_RATIO`, read when the belt origin wiring lands |
| row | 6.97:BG3 | GEN.1092 | `belt_origin.MIN_BELT_MASS_EARTH`, read when the belt origin wiring lands |
| row | 6.97:BG4 | GEN.1092 | `belt_origin.SNOW_BUMP_WIDTH`, read when the belt origin wiring lands |
| row | 6.97:BG5 | GEN.1092 | `belt_origin.CHAOTIC_ZONE_COEFFICIENT`, read when the belt origin wiring lands |
| row | 6.97:BG7 | GEN.1092 | `belt_origin.RESONANCE_HALF_WIDTH`, read when the belt origin wiring lands |
| row | 6.97:BG8 | GEN.1092 | `belt_origin.STIR_ZONE_WIDTHS`, read when the belt origin wiring lands |
| row | 6.97:BG9 | GEN.1092 | `belt_origin.INSTABILITY_SURVIVAL_RANGE`, read when the belt origin wiring lands |
| row | 6.97:BG10 | GEN.1092 | `belt_origin.DESTROYED_RETAINED_SHARE`, read when the belt origin wiring lands |
| row | 6.97:BG11 | GEN.1092 | `belt_origin.DEBRIS_MIN_WIDTH`, read when the belt origin wiring lands |
| row | 6.97:BG12 | GEN.1092 | `belt_origin.COLLISIONAL_DEBRIS_PER_GYR`, read when the belt origin wiring lands |
| row | 6.97:BG13 | GEN.1092 | `belt_origin.COLLISIONAL_DEBRIS_MASS_EARTH`, read when the belt origin wiring lands |
| module | planetgen.galaxy.rogue_density | GEN.1111 | the smoothed star density a rogue is drawn from (lock-free, GEN.1115); the virtual-donor sector draw takes its Poisson mean from it and the low-zoom map reads `old_density` |
| row | 6.102:RSM1 | GEN.1111 | `rogue_density.RADIAL_SMOOTH_PC`, read when the wiring step lands |
| row | 6.102:RSM2 | GEN.1111 | `rogue_density.VERTICAL_SMOOTH_PC`, read when the wiring step lands |
| row | 6.102:RSM3 | GEN.1111 | `rogue_density.YOUNG_AGE_MYR`, read when the wiring step lands |
| row | 6.102:RSM4 | GEN.1111 | `rogue_density.RING_STEPS`, read when the wiring step lands |
| module | planetgen.physics.dwarf_population | GEN.1092 | the count, round mass, orbit-clearing status and Class D of a belt's stored bodies (lock-free, GEN.1090); the wiring step lists the dwarf planets as members of their belt |
| row | 6.95:BR6 | GEN.1092 | `dwarf_population.SOTER_MU_PLANET`, read when the wiring step lands |
| row | 6.95:BR7 | GEN.1092 | `dwarf_population.MARGOT_COEFFICIENT_EARTH, MARGOT_C`, read when the wiring step lands |
| row | 6.95:BR9 | GEN.1092 | `dwarf_population.ROUND_MASS_ICE_KG, ROUND_MASS_ROCK_KG`, read when the wiring step lands |
| row | 6.104:DPM1 | GEN.1092 | `dwarf_population.ICE_BODY_FRACTION`, read when the wiring step lands |
| module | planetgen.physics.field_dispersal | GEN.1103 | the radius, age, retirement, freeze-out, shear and bodies of an interstellar asteroid field (lock-free, GEN.1100); the AsteroidField columns (GEN.1103) read the radius and the retirement, the bring-to-epoch step (GEN.1107) retires fields |
| row | 6.105:IF15 | GEN.1103 | `field_dispersal.RETIREMENT_RADIUS_M`, read when the field columns land |
| row | 6.105:IF16 | GEN.1105 | `field_dispersal.OORT_A_KMS_KPC, SHEAR_VALID_YEARS`, read when the field members land |
| row | 6.105:IF17 | GEN.1105 | `field_dispersal.DEFAULT_DENSITY_KG_M3`, read when the field members land |
| module | planetgen.physics.field_rates | GEN.1104 | the physical birth rates of interstellar asteroid fields and the host rule of a supernova remnant (lock-free, GEN.1101); the supernova host rule (GEN.1104) draws `hosted_belt`, the columns (GEN.1103) take the standing numbers |
| row | 6.106:IF19 | GEN.1104 | `field_rates.GALAXY_CORE_COLLAPSE_RATE_PER_YEAR`, read when the host rule lands |
| row | 6.106:IF20 | GEN.1104 | `field_rates.TRUNCATION_FRACTION`, read when the cluster swarm hook lands |
| row | 6.106:IF21 | GEN.1104 | `field_rates.CLUSTER_BELT_SHARE`, read when the cluster swarm hook lands |
| module | planetgen.physics.rogue_topup | GEN.1120 | the stage 3 top-up of the region run: how many rogues, their masses, places, velocities and apparent origins, and the remainder for the asteroid fields |
| row | 6.86:TU1 | GEN.1120 | a constant of the rogue top-up planner, read when the region run wires it |
| row | 6.86:TU2 | GEN.1120 | a constant of the rogue top-up planner, read when the region run wires it |
| row | 6.86:TU3 | GEN.1120 | a constant of the rogue top-up planner, read when the region run wires it |
| row | 6.86:TU4 | GEN.1120 | a constant of the rogue top-up planner, read when the region run wires it |
| row | 6.86:TU5 | GEN.1120 | a constant of the rogue top-up planner, read when the region run wires it |
