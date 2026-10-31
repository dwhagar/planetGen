# The wiring register: what is built or downloaded and not read yet (GEN.754)

Status: hand-kept, checked by `scripts/wiring_guard.py` and `src/tests/test_wiring_guard.py`. Every module, setup dataset, dataset reader, data file and stand-in registry row (GEN.725) that nothing in the app reaches is named here with the open item that will wire it, or with the words `by design` and the reason. The test fails when something unread is missing from this page, when a line names something that is read now (delete the line in the same change that wires it), and when a named item is not open in `docs/TODO.md`. No shim: an unread piece is wired, filed or deleted.

What counts as read: the app is every module outside `physics`, `generation`, `galaxy`, `datasets` and `tuning` (the CLI, API, web, admin, queue, population, backup and describe packages). A library module is read when an import chain from the app reaches it, a lazy import inside a function counts, and so does a module named by a string (`importlib.import_module("planetgen.x")`, the registry's `"planetgen.x:function"`). A dataset reader is read when a reachable module outside `datasets` uses it. A data file is read when a reachable module names it. A registry row is read when one of its `module.CONSTANT` names is used by a reachable module.

`python scripts/wiring_guard.py` prints what the guard finds; `--check` is the test.

| kind | name | item | reason |
|---|---|---|---|
| dataset | harris_clusters | GEN.297 | the Harris catalog of 157 globular clusters; the calibration of the cluster generator reads it |
| dataset | noaa_proton_events | GEN.805 | the 265 NOAA solar proton events; the solar particle event model fits its rate and slope to them |
| module | planetgen.physics.alpha_tracks | GEN.744 | the equivalent [Fe/H] of an alpha-enhanced star; the track lookup of a high-alpha star will read it |
| row | 6.34:AL1 | GEN.744 | `alpha_tracks.ALPHA_SHARE`, not read until the track lookup is wired |
| module | planetgen.galaxy.thick_ages | GEN.1025 | the old population's age law with the thick disc's 8 to 13 Gy; part 2 wires it into the galaxy's age laws |
| row | 6.146:TD1 | GEN.1025 | `tuning.THICK_DISC_AGE_RANGE_GY`, read by `thick_ages` only until part 2 |
| module | planetgen.physics.chemically_peculiar | GEN.1027 | the Ap and Am kinds, slow rotation and field; part 2 draws them in the star build |
| row | 6.147:PEC1 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
| row | 6.147:PEC2 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
| row | 6.147:PEC3 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
| row | 6.147:PEC4 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
| row | 6.147:PEC5 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
| row | 6.147:PEC6 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
| row | 6.147:PEC7 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
| row | 6.147:PEC8 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
| row | 6.147:PEC9 | GEN.1027 | `chemically_peculiar` constant, not read until part 2 |
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
| module | planetgen.physics.latitude_bands | GEN.1108 | the climate of a body in seven latitude bands (hours of day and night, polar day and night, day and night temperatures, state), a pure function; the storage part replaces the four stored bands of `regional_temperature` with it |
| module | planetgen.physics.radiation_belts | GEN.983 | trapped radiation belts of a magnetised body, read when the stored magnetosphere reaches dose |
| module | planetgen.galaxy.influence | GEN.342 | the influence set of an object, read by the acceleration kernel of the orbital update |
| module | planetgen.physics.closest_approach | GEN.197 | the pair closest-approach library, read by the orbital event writer |
| module | planetgen.physics.body_positions | MAP.253 | where every body of a system is at a time, read by the scene and the maps once they read the parent chain |
| module | planetgen.physics.light_travel | VIEW.1 | light-travel positions for the view from a planet |
| module | planetgen.physics.radius_valley | by design | the simulated Kepler survey of the radius valley, a check run by the tuning tests; the valley itself comes from the envelope-loss rows of the generator (GEN.722 retired) |
| module | planetgen.physics.iron_branch_tuning | by design | the simulated sample the iron branch table is tuned to; run by the tuning tests, not by generation |
| module | planetgen.generation.redraw_kinds | by design | the names and one-line meanings of the redraw counters (TEST.476); `draw.redrawn` is called by name from the generator and the census reads this registry, so no generation code imports it |
| module | planetgen.galaxy.designations | GEN.1138 | the pure rules of designation history (retired IDs reserved, cascade, resolution), read by the writers and the lookup when the table lands |
| module | planetgen.generation.class_sample | by design | a sample that holds every class, a helper of the sheet audit tests |
| module | planetgen.generation.plausibility | by design | the anomaly finder engine; run by the plausibility tests and their command line wrapper in `src/tests` |
| module | planetgen.generation.phenomena_plausibility | by design | the exotic-phenomena counterpart of the anomaly finder, run by its tests |
| row | 6.7:AM18 | GEN.524 | `mixture_trace.ABUNDANCES`, not drawn until the trace gases are wired |
| row | 6.10:RV2 | by design | `radius_valley.CORE_MASS_SIGMA`, a width of the simulated survey, read by the tuning tests |
| row | 6.10:RV3 | by design | `radius_valley.ENVELOPE_SIGMA`, a width of the simulated survey, read by the tuning tests |
| row | 6.15:RB1 | GEN.983 | `radiation_belts.WIND_EDGE_FRACTION`, read when the belts are |
| row | 6.15:RB2 | GEN.983 | `radiation_belts.ROTATION_EDGE_FACTOR`, read when the belts are |
| row | 6.15:RB3 | GEN.983 | `radiation_belts.WIND_PEAK_FRACTION`, read when the belts are |
| row | 6.15:RB4 | GEN.983 | `radiation_belts.INNER_PEAK_FRACTION`, read when the belts are |
| row | 6.15:RB5 | GEN.983 | `radiation_belts.ATMOSPHERE_EDGE_SCALE_HEIGHTS`, read when the belts are |
| row | 6.15:RB6 | GEN.983 | `radiation_belts.EARTH_OUTER_PEAK_GY_YR`, read when the belts are |
| row | 6.15:RB7 | GEN.983 | `radiation_belts.OUTER_RISE_EXPONENT`, read when the belts are |
| module | planetgen.physics.planet_occurrence | by design | the reader of the observed occurrence bands; the architecture and the census read it (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
| module | planetgen.physics.kepler_multiplicity | by design | the Kepler transit multiplicity harness, read by the tuning tests until the wiring (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
| data | planet_occurrence.json | by design | read by `physics/planet_occurrence.py` only (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
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
| row | 6.63:CF3 | by design | `cloud_infrared.TOP_OVER_BASE_PRESSURE`, read when the column solver calls the cloud (read by the radiative-convective model, which the generator reads through its stored table, like line_opacity; GEN.850 built) |
| module | planetgen.physics.clearance_ledger | GEN.1014 | the spans where a planet, moon or belt may sit; replaces `_orbit_floor_au`, the per-star `a_crit_au` reads and the spacing and `min_orbit_distance` pushes of `space_orbits` once the system build draws from it (train 10, after GEN.1009 and GEN.1011) |
| data | quadruple_stability_3p1.json | GEN.1016 | the 3+1 quadruple stability rules (validity ranges, the stability definition); `multistar` limits will read them through `quadruple_stability` in the wiring step |
| data | quadruple_stability_3p1_weights.json | GEN.1016 | the 8,301 numbers of the authors' network (made by `scripts/build_quadruple_stability.py`); read by `quadruple_stability` |
| data | quadruple_stability_2p2.json | GEN.1016 | the 2+2 quadruple stability rules (validity ranges, the stability definition); read by `quadruple_stability` |
| data | quadruple_stability_2p2_weights.json | GEN.1016 | the 8,301 numbers of the authors' 2+2 network (made by `scripts/build_quadruple_stability.py`); read by `quadruple_stability` |
| module | planetgen.physics.quadruple_stability | GEN.1016 | the stability verdict of a 3+1 quadruple and (GEN.1031) of a 2+2 quadruple; the stability margin and quadruple fit (P10) will call it |
| module | planetgen.physics.hierarchy_stability | GEN.1016 | the stability of a tree of stars level by level (triple levels, the 3+1 quadruple with its margin, 2+2 by nested triples); read by tests/test_hierarchy_stability.py now and by the `multistar` limits once the wiring merges on its train |
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
| module | planetgen.physics.belt_density | UX.241 | the objects per square AU above a diameter, the mass per square AU and the size-table bins of a belt; the belt page data builder reads it |
| row | 6.99:BD1 | UX.241 | a constant of the belt density module, read when its wiring step lands |
| row | 6.99:BD2 | UX.241 | a constant of the belt density module, read when its wiring step lands |
| row | 6.99:BD3 | UX.241 | a constant of the belt density module, read when its wiring step lands |
| row | 6.99:BD4 | UX.241 | a constant of the belt density module, read when its wiring step lands |
| row | 6.99:BD5 | UX.241 | a constant of the belt density module, read when its wiring step lands |
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
| row | 6.87:RK3 | GEN.1111 | a constant of the rogue kick module, read when its wiring step lands |
| row | 6.87:RK4 | GEN.1111 | a constant of the rogue kick module, read when its wiring step lands |
| row | 6.87:RK5 | GEN.1111 | a constant of the rogue kick module, read when its wiring step lands |
| row | 6.87:RK6 | GEN.1111 | a constant of the rogue kick module, read when its wiring step lands |
| row | 6.87:RK7 | GEN.1111 | a constant of the rogue kick module, read when its wiring step lands |
| row | 6.87:RK8 | GEN.1111 | a constant of the rogue kick module, read when its wiring step lands |
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
| module | planetgen.physics.rogue_radiogenic | GEN.1111 | the radiogenic content and heat of a rogue from the star it left; the rogue rows read them in place of the host-less draw |
| row | 6.114:RY1 | GEN.1111 | a constant of the rogue radiogenic module, read when its wiring step lands |
| module | planetgen.physics.super_puff | GEN.1048 | the super-puff radius, fade and gate (W15); the bulk step reads it when the wiring step lands |
| row | 6.76:PF1 | GEN.1048 | a constant of the super-puff module, read when its wiring step lands |
| row | 6.76:PF2 | GEN.1048 | a constant of the super-puff module, read when its wiring step lands |
| row | 6.76:PF3 | GEN.1048 | a constant of the super-puff module, read when its wiring step lands |
| row | 6.76:PF4 | GEN.1048 | a constant of the super-puff module, read when its wiring step lands |
| module | planetgen.physics.helium_fractionation | GEN.1047 | the helium share of an envelope that loses hydrogen first (W16); the oddity row reads it when the oddity wiring lands |
| row | 6.77:HF1 | GEN.1047 | a constant of the helium fractionation module, read when its wiring step lands |
| row | 6.77:HF2 | GEN.1047 | a constant of the helium fractionation module, read when its wiring step lands |
| row | 6.77:HF3 | GEN.1047 | a constant of the helium fractionation module, read when its wiring step lands |
| module | planetgen.physics.mantle_stripping | GEN.1048 | the iron-rich core of a giant impact (Marcus and others 2010 scaling); the weirdness wiring (W10) gives a rocky planet that had a recorded merger a core fraction from the impacts the architecture recorded |
| row | 6.93:MT4 | GEN.1048 | `mantle_stripping.Q_STAR_SCALE` of the mantle stripping module, read when its wiring step lands |
| row | 6.93:MT6 | GEN.1048 | `mantle_stripping.MAX_CORE_MASS_FRACTION` of the mantle stripping module, read when its wiring step lands |
| row | 6.93:MT7 | GEN.1048 | `mantle_stripping.IMPACT_ANGLE_DEG` of the mantle stripping module, read when its wiring step lands |
| module | planetgen.physics.hycean_candidate | GEN.1047 | the Hycean box as a flag on the classes T, I and V; the oddity line names the planets it selects when the oddity wiring lands |
| row | 6.79:HC1 | GEN.1047 | the four ranges of the Hycean box, read when its wiring step lands |
| module | planetgen.physics.irregular_zone | GEN.1083 | the stable limits, Kozai band and inner edge of the irregular zone; the irregular orbit draw (GEN.1078) and the captured-moon wiring will call it |
| row | 6.94:IZ1 | GEN.1083 | a constant of the irregular zone module, read when its wiring step lands |
| row | 6.94:IZ2 | GEN.1083 | a constant of the irregular zone module, read when its wiring step lands |
| row | 6.94:IZ3 | GEN.1083 | a constant of the irregular zone module, read when its wiring step lands |
| module | planetgen.physics.rogue_collision | GEN.1103 | the rate, disruption speed, largest remnant, fragment speeds and shock of a rogue-rogue collision; the field columns and host rule (GEN.1103, GEN.1104) call it |
| module | planetgen.physics.irregular_orbits | GEN.1083 | the a / r_H, eccentricity, inclination and family draw of irregular moons; the captured-moon wiring will call it |
| row | 6.96:IO1 | GEN.1083 | a constant of the irregular orbits module, read when its wiring step lands |
| row | 6.96:IO2 | GEN.1083 | a constant of the irregular orbits module, read when its wiring step lands |
| row | 6.96:IO3 | GEN.1083 | a constant of the irregular orbits module, read when its wiring step lands |
| row | 6.96:IO4 | GEN.1083 | a constant of the irregular orbits module, read when its wiring step lands |
| module | planetgen.physics.belt_ledger | GEN.1069 | the stored bodies of a belt by the mass cut, the size laws and the largest body settled with the collisional mass (lock-free, GEN.1061); the belt store and the belt generator call it, the outer disc (GEN.1071) and the scatter step (GEN.1092) after them |
| row | 6.95:BR17 | GEN.1069 | `belt_ledger.L1_MASS_KG`, read when the wiring step lands |
| row | 6.95:BR18 | GEN.1069 | `belt_ledger.LEDGER_COVERAGE, LEDGER_MAX_BODIES`, read when the wiring step lands |
| module | planetgen.physics.belt_member | GEN.473 | the taxon, differentiation, grains, density and radius of one stored belt member from the kernel at its own distance (lock-free, GEN.1067); the belt wiring of the taxa calls it for each member of the ledger |
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
| module | planetgen.galaxy.resident_mass | GEN.1120 | the stage 0 budget of a sector once its residents are off, and the stage 3 shortfall; the generation phases read it |
| module | planetgen.physics.field_composition | GEN.1103 | the layers, escaping shards, remnant core fraction, collision outcome and makeup of a collision field (lock-free, GEN.1102); the field generation (GEN.1103) calls `collision_field` and `field_makeup` in place of `draw_field("collisional")` |
| row | 6.107:IF23 | GEN.1103 | `field_composition.CORE_LEAK_SHARE`, read when the field generation lands |
| row | 6.107:IF24 | GEN.1103 | `field_composition.MANTLE_ROWS, CORE_ROW`, read when the field generation lands |
| module | planetgen.physics.topup_field | GEN.1120 | the light top-up asteroid fields of the remainder the Rogues level leaves unmatched (lock-free, GEN.1126); the region run's stage 3 calls `plan` with `rogue_topup.TopUp.remainder_kg`, and the field items (GEN.1103 to GEN.1107) store and retire what it returns |
| row | 6.108:TF2 | GEN.1120 | `topup_field.LARGEST_BODY_MASS_KG, BULK_DENSITY_KG_M3`, read when the region run wires stage 3 |
| row | 6.108:TF3 | GEN.1120 | `topup_field.SPEED_KM_S_RANGE, MIN_RADIUS_LY`, read when the region run wires stage 3 |
| row | 6.108:TF4 | GEN.1120 | `topup_field.MIN_FIELD_MASS_KG, MAX_FIELD_MASS_KG, MAX_FIELDS`, read when the region run wires stage 3 |
| module | planetgen.physics.resonance_gaps | GEN.1069 | the resonance zones, chaotic zones, survivors, Yarkovsky cut and stirring of a belt from its giants (lock-free, GEN.1062); the belt generator reads the map and the clearance ledger (GEN.1014) can reuse the zones |
| row | 6.116:RGP1 | GEN.1069 | `resonance_gaps.CHAOTIC_ZONE_COEFFICIENT, EXTENDED_ZONE_COEFFICIENT, CRITICAL_ECCENTRICITY_COEFFICIENT, MAX_ECCENTRICITY_COEFFICIENT, LOW_ECCENTRICITY_WIDTH, HIGH_ECCENTRICITY_WIDTH`, read when the wiring step lands |
| row | 6.116:RGP2 | GEN.1069 | `resonance_gaps.ORDER_WIDTH_EXPONENT`, read when the wiring step lands |
| row | 6.116:RGP3 | GEN.1069 | `resonance_gaps.RESONANCE_CLEARING_MYR`, read when the wiring step lands |
| row | 6.116:RGP4 | GEN.1069 | `resonance_gaps.MAX_ORDER, MAX_P`, read when the wiring step lands |
| row | 6.116:RGP5 | GEN.1069 | `resonance_gaps.FORCED_ECCENTRICITY_COEFFICIENT`, read when the wiring step lands |
| row | 6.116:RGP6 | GEN.1069 | `resonance_gaps.YARKOVSKY_COEFFICIENT_AU_KM_PER_MYR, YARKOVSKY_REFERENCE_AU`, read when the wiring step lands |
| row | 6.116:RGP7 | GEN.1069 | `resonance_gaps.ZONE_ECCENTRICITY_FLOOR`, read when the wiring step lands |
| module | planetgen.physics.belt_history | GEN.1069 | the events that act on a belt (giant migration or instability, the giant branch, mass loss, a companion) and the pollution rate of the rock it scatters (lock-free, GEN.1064); the belt generator runs the events a system's history holds and the white dwarf pollution row reads the rate |
| row | 6.117:BHS1 | GEN.1069 | `belt_history.GRAND_TACK_SURVIVAL, INSTABILITY_SURVIVAL`, read when the wiring step lands |
| row | 6.117:BHS2 | GEN.1069 | `belt_history.SURFACE_DENSITY_SLOPE`, read when the wiring step lands |
| row | 6.117:BHS3 | GEN.1069 | `belt_history.MASS_LOSS_SCATTER_SHARE`, read when the wiring step lands |
| row | 6.117:BHS4 | GEN.1069 | `belt_history.DELIVERY_MYR`, read when the wiring step lands |
| row | 6.117:BHS5 | GEN.1069 | `belt_history.MIN_POLLUTION_RATE_G_S`, read when the wiring step lands |
| module | planetgen.physics.irregular_population | GEN.1083 | the capture windows, ground-down count and stored irregulars of a planet (lock-free, GEN.1077); the captured-moon wiring calls `make_pool`, `draw_irregulars` and `draw_exchange` |
| row | 6.118:IP1 | GEN.1083 | `irregular_population.REFERENCE_DISC_MASS_EARTH`, read when the captured-moon wiring lands |
| row | 6.118:IP2 | GEN.1083 | `irregular_population.CAPTURED_DENSITY_KG_M3`, read when the captured-moon wiring lands |
| row | 6.118:IP3 | GEN.1083 | `irregular_population.ENCOUNTER_YIELD`, read when the captured-moon wiring lands |
| row | 6.118:IP4 | GEN.1083 | `irregular_population.REFERENCE_HILL_AU`, read when the captured-moon wiring lands |
| row | 6.118:IP5 | GEN.1083 | `irregular_population.GAS_WINDOW_STRENGTH`, read when the captured-moon wiring lands |
| row | 6.118:IP6 | GEN.1083 | `irregular_population.GAS_SIZE_WINDOW_KM`, read when the captured-moon wiring lands |
| row | 6.118:IP7 | GEN.1083 | `irregular_population.EXCHANGE_CHANCE`, read when the captured-moon wiring lands |
| row | 6.118:IP8 | GEN.1083 | `irregular_population.EXCHANGE_MIN_RADIUS_KM`, read when the captured-moon wiring lands |
| row | 6.118:IP9 | GEN.1083 | `irregular_population.TEMPORARY_MIN_PLANET_RADII`, read when the captured-moon wiring lands |
| row | 6.118:IP10 | GEN.1083 | `irregular_population.DISRUPTION_SIZE_RATIO`, read when the captured-moon wiring lands |
| row | 6.118:IP11 | GEN.1083 | `irregular_population.MODEL_MIN_RADIUS_KM`, read when the captured-moon wiring lands |
| row | 6.118:IP12 | GEN.1083 | `irregular_population.STORE_FLOOR_KG`, read when the captured-moon wiring lands |
| module | planetgen.physics.captured_composition | GEN.1083 | the taxon, ice and density of a captured body from its source and a family's lumps (lock-free, GEN.1080); the captured-moon wiring calls `compose`, `draw_parent` and `fragment` |
| row | 6.119:CP1 | GEN.1083 | `captured_composition.FEEDING_ZONE_HILL_RADII`, read when the captured-moon wiring lands |
| row | 6.119:CP3 | GEN.1083 | `captured_composition.PARENT_LUMPS`, read when the captured-moon wiring lands |
| row | 6.119:CP2 | GEN.1083 | `belt_history.SURFACE_DENSITY_SLOPE`, read by the captured composition when the captured-moon wiring lands |
| module | planetgen.physics.belt_page | UX.241 | the data of a belt's information page (lock-free, UX.240); the page route and the API field call `build` |
| row | 6.120:BP1 | UX.241 | `belt_page.BIN_EDGES_KM, BIN_LOW_KM`, read when the page route lands |
| row | 6.120:BP2 | UX.241 | `belt_page.DENSITY_DIAMETERS_KM`, read when the page route lands |
| row | 6.120:BP3 | UX.241 | `belt_page.L2_FLOOR_KG`, read when the page route lands |
| row | 6.120:BP4 | UX.241 | `belt_page.SHARE_SAMPLES, MIN_SHARE`, read when the page route lands |
| module | planetgen.physics.trojans | GEN.1150 | what may sit at L4 and L5 of a pair: mass rule, libration, swarm split, map label (lock-free); the generator draw and the map icon call it |
| row | 6.121:TJ1 | GEN.1150 | `trojans.MU_CRIT`, read when the swarm draw lands |
| row | 6.121:TJ3 | GEN.1150 | `trojans.libration_period, TADPOLE_LIMIT_DEG`, read when the orbit overlay lands |
| row | 6.121:TJ4 | GEN.1150 | `trojans.LEADING_SHARE_RANGE`, read when the swarm draw lands |
| row | 6.121:TJ5 | GEN.1150 | `trojans.ECC_MAX, INCLINATION_MAX_DEG`, read when the swarm draw lands |
| module | planetgen.physics.trojan_bodies | GEN.1151 | which classes of world the mass rule admits at L4 and L5, and the formation gate (lock-free); the generator draw calls it |
| row | 6.125:TB1 | GEN.1151 | `trojan_bodies.CLASS_J_MIN_EARTH`, read when the whole-body draw lands |
| row | 6.125:TB2 | GEN.1151 | `trojan_bodies.CLASS_M_MASS_EARTH`, read when the whole-body draw lands |
| row | 6.125:TB3 | GEN.1151 | `trojan_bodies.FORMATION_EFFICIENCY`, read when the whole-body draw lands |
| row | 6.125:TB4 | GEN.1151 | `trojan_bodies.MIN_DISC_LIFETIME_MYR`, read when the whole-body draw lands |
| module | planetgen.physics.lagrange_geometry | MAP.354 | the five Lagrange points of a pair in the system frame and the label of each occupied one (lock-free); the System Map and the gravity layer (MAP.169) call it |
| module | planetgen.physics.climate_branches | GEN.972 | the ice-albedo law and every fixed point of the surface temperature with its stability, a pure function; GEN.972 puts the law into the fixed point and GEN.1051 hands both branches to pick_branch |
| row | 6.22:CC15 | GEN.972 | `climate_branches.ICE_FREE_ALBEDO, ICE_COVERED_ALBEDO`, read when GEN.972 wires the module |
| row | 6.22:CC16 | GEN.972 | `climate_branches.ICE_FREE_ABOVE_K, ICE_COVERED_BELOW_K, ICE_FRACTION_EXPONENT`, read when GEN.972 wires the module |
| module | planetgen.physics.roundness | GEN.1159 | the yield strength, round radius and roundness ratio of a body; Class D and C, the Asteroid type and the stored columns read it once GEN.1159 and GEN.1160 wire it (lock) |
| module | planetgen.physics.binding_energy | GEN.1159 | the gravitational binding energy of every kind of mass object, a black hole stored as -1; the stored columns read it once GEN.1159 and GEN.1162 wire it (lock) |
| data | roundness_strengths.json | GEN.1159 | read by `physics/roundness.py` only |
| row | 6.126:RN1 | GEN.1159 | `roundness.SPHERICITY_FACTOR`, not read until GEN.1159 wires the module |
| row | 6.126:RN2 | GEN.1159 | `roundness.ICE_STRENGTH_PA, ROCK_STRENGTH_PA`, not read until GEN.1159 wires the module |
| row | 6.126:RN3 | GEN.1159 | the volume-threshold mixing rule, not read until GEN.1159 wires the module |
| row | 6.126:RN4 | GEN.1159 | `roundness.METAL_STRENGTH_PA`, not read until GEN.1159 wires the module |
| row | 6.126:RN5 | GEN.1159 | the porosity rules, not read until GEN.1159 wires the module |
| row | 6.126:RN6 | GEN.1159 | `roundness.GRAIN_DENSITIES_KG_M3`, not read until GEN.1159 wires the module |
| row | 6.126:RN7 | GEN.1159 | `roundness.TEMPERATURE_RULES`, not read until GEN.1159 wires the module |
| row | 6.126:BE1 | GEN.1159 | `binding_energy.FULLY_CONVECTIVE_MASS_SOLAR`, not read until GEN.1162 wires the stars |
| row | 6.126:BE2 | GEN.1159 | `binding_energy.POLYTROPE_INDEX_GAS_GIANT`, not read until GEN.1159 wires the module |
| row | 6.126:BE3 | GEN.1159 | `binding_energy.POLYTROPE_INDEX_CONVECTIVE`, not read until GEN.1162 wires the stars |
| row | 6.126:BE4 | GEN.1159 | `binding_energy.POLYTROPE_INDEX_RADIATIVE`, not read until GEN.1162 wires the stars |
| row | 6.126:RN8 | GEN.1159 | `roundness.RUBBLE_RULES` (friction angle and cohesion), not read until GEN.1159 wires the module |
| row | 6.126:RN9 | GEN.1159 | `roundness.RUBBLE_RULES` (hot strength), not read until GEN.1159 wires the module |
| row | 6.126:RN10 | GEN.1178 | `roundness.RANGES_PA`, not read until the strength draw is wired into `draw.py` with GEN.1159 |
| module | planetgen.physics.ice_types | GEN.1159 | the strength of the other ices next to ice Ih and the ice type of a body from its components (lock-free); `roundness.py` takes the ratio as a factor on its ice strength when the types are wired |
| row | 6.126:RN11 | GEN.1201 | `roundness.ROCK_CREEP_COLUMNS`, not read until GEN.1201 gives each body a shape-set temperature to pass to the engine |
| module | planetgen.physics.ring_darkening | GEN.1055 | how bright a ring is for its age and the dust that falls on it (R6) and how much sunlight a ring shades from a latitude over the lit half of the day (R8); `rings.apply` and the climate read it when GEN.1055 wires it |
| row | 6.148:RGD1 | GEN.1055 | `ring_darkening.CLEAN_ICE_ALBEDO` and `HALF_DARKENING_SHARE`, not read until GEN.1055 wires the module |
| row | 6.148:RGD2 | GEN.1055 | `ring_darkening.DARK_ALBEDO`, not read until GEN.1055 wires the module |
| row | 6.148:RGD3 | GEN.1055 | `ring_darkening.SATURN_POLLUTION_FLUX_KG_M2_S`, not read until GEN.1055 wires the module |
| row | 6.148:RGD4 | GEN.1055 | `ring_darkening.FLUX_DISTANCE_EXPONENT` and `DUST_SPEED_AT_INFINITY_KMS`, not read until GEN.1055 wires the module |
| data | planet_class_fit.json | GEN.1207 | the percent-match table of every class; the shadow scorer reads it (GEN.1206) and the best-fit classifier makes it the way a stored letter is decided |
| module | planetgen.physics.class_fit | GEN.1207 | the percent-match scorer; the shadow scorer in src/tests reads it, the classifier replaces the first-match tree with it |
| module | planetgen.physics.planet_type | GEN.1219 | the terrestrial, transitional or gas type from the envelope share of the mass; the shadow report in src/tests reads it, GEN.1219 stores the type with GEN.1207's epoch |
| row | 6.156:PT1 | GEN.1219 | read by planet_type.type_of; the stored type waits for GEN.1207's epoch |
| module | planetgen.physics.life_readings | UX.251 | the abundance, trend and security words of each kind of life on a world, from its stored values; the Life summary and the Life details of a planet or moon page, and the API, read it when UX.251 wires it |
| row | 6.155:LR1 | UX.251 | `life_readings.ABUNDANCE_BANDS`, not read until UX.251 wires the module |
| row | 6.155:LR2 | UX.251 | `life_readings.TREND_SPAN_GYR`, `TREND_BAND` and `SPREAD_GYR`, not read until UX.251 wires the module |
| row | 6.155:LR3 | UX.251 | `life_readings.DYING_GYR`, `PRECARIOUS_GYR`, `MARGINAL_GYR` and `MARGINAL_SCORE`, not read until UX.251 wires the module |
| row | 6.155:LR4 | UX.251 | `life_readings.TEMPERATURE_EXPONENT`, not read until UX.251 wires the module |
| row | 6.155:LR5 | UX.251 | `life_readings.KIND_MILESTONES`, not read until UX.251 wires the module |
