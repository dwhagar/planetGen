# The wiring register: what is built or downloaded and not read yet (GEN.754)

Status: hand-kept, checked by `scripts/wiring_guard.py` and `src/tests/test_wiring_guard.py`. Every module, setup dataset, dataset reader, data file, stand-in registry row (GEN.725), stored column and public function (TEST.520) that nothing in the app reaches is named here with the open item that will wire it, or with the words `by design` and the reason. The test fails when something unread is missing from this page, when a line names something that is read now (delete the line in the same change that wires it), and when a named item is not open in `docs/TODO.md`. No shim: an unread piece is wired, filed or deleted.

What counts as read: the app is every module outside `physics`, `generation`, `galaxy`, `datasets` and `tuning` (the CLI, API, web, admin, queue, population, backup and describe packages). A library module is read when an import chain from the app reaches it, a lazy import inside a function counts, and so does a module named by a string (`importlib.import_module("planetgen.x")`, the registry's `"planetgen.x:function"`). A dataset reader is read when a reachable module outside `datasets` uses it. A data file is read when a reachable module names it. A registry row is read when one of its `module.CONSTANT` names is used by a reachable module. A `column` is a column of a table that holds objects (planets, moons, stars, systems, belts, comets, fields, rogue planets and the other tables of the `uid` group of `db/kind_tables.py`) and is read when a module outside the plumbing (`PLUMBING` in the guard: the schema, models, column units, store, storage checksum, checks, export and repair code and the generic Details sheet) loads its name as an attribute, a name, a dictionary key or a word of a SQL string; a module of physics or generation that loads the value is a reader, so a column read downstream is not listed. A `function` is a public module-level function of a reachable module of `physics`, `generation` or `galaxy` that no other module, no other function of its module and no string of the package names; only tests call it. The test is by name, so a name two modules share counts as called. A function with a decorator that could register it (anything but a cache or a property) is left out. Wire the piece, delete it, or add a line here with the open item that will; delete the line in the same change that wires it.

`python scripts/wiring_guard.py` prints what the guard finds; `--check` is the test.

| kind | name | item | reason |
|---|---|---|---|
| dataset | harris_clusters | GEN.297 | the Harris catalog of 157 globular clusters; the calibration of the cluster generator reads it |
| dataset | noaa_proton_events | GEN.805 | the 265 NOAA solar proton events; the solar particle event model fits its rate and slope to them |
| module | planetgen.physics.line_opacity | by design | the band opacities of the radiative-convective model, which the greenhouse table is fitted to and the checks hold it to; the generator reads the table, never the model |
| module | planetgen.physics.radiative_convective | by design | the model the greenhouse table (`greenhouse_table.py`) is fitted to by `scripts/fit_greenhouse_table.py`, too slow for the fixed point; run by the fit and the tests, not by generation |
| module | planetgen.physics.radiative_kernel | by design | the numba sweep of `radiative_convective`; the web app must not load numba, so only the fit and the tests reach it |
| module | planetgen.physics.clearance | by design | the independent overlap verifier (TEST.477); the tests and the redraw audit read it, the generator never does, so the proof cannot share a bug with the code it checks |
| module | planetgen.physics.photolysis | GEN.982 | the photolysis lifetime table, now with the mean UV cross sections (Hu 2012 table 6) and the photon-limited loss; part 3 of GEN.982 (the methane haze supply in `aerosols.haze_in`) reads it |
| function | planetgen.generation.remnant_orbits.pollution_rate_g_s | GEN.1231 | the white dwarf pollution flag reads the rate of its belt (the scattered mass over the delivery time) when belt birth sees the giants (GEN.1068); the census band of 25 to 50 percent reads it too |
| function | planetgen.physics.activity.far_uv_flux_ratio | GEN.982 | the star's 91.2 to 200 nm flux at the body relative to the quiet Sun at 1 AU, the `flux_ratio` of the photolysis rates; part 3 of GEN.982 hands it to the haze supply through `surface_air.Setting` |
| data | photolysis_lifetimes.json | GEN.982 | read by `physics/photolysis.py` only |
| module | planetgen.physics.mixture_trace | GEN.524 | the trace gases of an air drawn after the life stage; no caller yet |
| module | planetgen.physics.design_flags | GEN.524 | gates the trace rows `mixture_trace` draws |
| data | design_flags.json | GEN.524 | read by `physics/design_flags.py` only |
| module | planetgen.physics.equilibrium_edges | GEN.994 | the edges of the hot-Jupiter chemistry regimes, sourced and wired with the hot-giant thermochemistry |
| module | planetgen.physics.ground_albedo | GEN.684 | the colour of the ground by band, wired by the ground albedo item |
| module | planetgen.galaxy.influence | GEN.342 | the influence set of an object, read by the acceleration kernel of the orbital update |
| module | planetgen.galaxy.influence | GEN.105 | the influence set of an object, read by the orbital update run through the acceleration kernel |
| module | planetgen.galaxy.acceleration | GEN.105 | the acceleration kernel (influence set plus galactic gradient) and the velocity kick, read by the orbital update run |
| module | planetgen.galaxy.stepper | GEN.105 | the kick-drift-kick star stepper on arrays, read by the orbital update run in place of the turn about the pole |
| row | 6.159:ST1 | GEN.105 | `stepper.STEP_ETA` and `STEP_CAP_MYR`, read when the orbital update run steps stars |
| row | 6.159:ST2 | GEN.105 | `stepper.MAX_LEVEL`, read when the orbital update run steps stars |
| row | 6.158:AK1 | GEN.105 | `acceleration.MEMBER_SOFTENING_PC`, read when the orbital update run calls the kernel |
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
| module | planetgen.physics.planet_occurrence | by design | the reader of the observed occurrence bands; the architecture and the census read it (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
| module | planetgen.physics.kepler_multiplicity | by design | the Kepler transit multiplicity harness, read by the tuning tests until the wiring (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
| data | planet_occurrence.json | by design | read by `physics/planet_occurrence.py` only (read by the tuning tests and the census harness; the architecture was fitted to it, GEN.728 built) |
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
| module | planetgen.galaxy.rogue_density | GEN.1111 | the smoothed star density a rogue is drawn from (lock-free, GEN.1115); the virtual-donor sector draw takes its Poisson mean from it and the low-zoom map reads `old_density` |
| row | 6.102:RSM1 | GEN.1111 | `rogue_density.RADIAL_SMOOTH_PC`, read when the wiring step lands |
| row | 6.102:RSM2 | GEN.1111 | `rogue_density.VERTICAL_SMOOTH_PC`, read when the wiring step lands |
| row | 6.102:RSM3 | GEN.1111 | `rogue_density.YOUNG_AGE_MYR`, read when the wiring step lands |
| row | 6.102:RSM4 | GEN.1111 | `rogue_density.RING_STEPS`, read when the wiring step lands |
| module | planetgen.physics.dwarf_population | GEN.1092 | the count, round mass, orbit-clearing status and Class D of a belt's stored bodies (lock-free, GEN.1090); the wiring step lists the dwarf planets as members of their belt |
| row | 6.95:BR6 | GEN.1092 | `dwarf_population.SOTER_MU_PLANET`, read when the wiring step lands |
| row | 6.95:BR7 | GEN.1092 | `dwarf_population.MARGOT_COEFFICIENT_EARTH, MARGOT_C`, read when the wiring step lands |
| module | planetgen.physics.field_dispersal | GEN.1103 | the radius, age, retirement, freeze-out, shear and bodies of an interstellar asteroid field (lock-free, GEN.1100); the AsteroidField columns (GEN.1103) read the radius and the retirement, the bring-to-epoch step (GEN.1107) retires fields |
| module | planetgen.physics.field_mass | GEN.1103 | the mass of an interstellar field by its origin (a supernova-unbound belt by the belts' own collisional step, two colliding rogue planets by their unbound share) and its stored bodies by the belt ledger (lock-free pure part of GEN.1161); the field columns (GEN.1103) store the mass, the supernova host rule (GEN.1104) calls `unbound_belt_field` and the member rows and designations (GEN.1105) call `stored_bodies` |
| row | 6.161:FM1 | GEN.1103 | `field_mass.MIN_FIELD_MASS_EARTH`, read when the field columns land |
| row | 6.161:FM2 | GEN.1105 | `field_mass.COLLISION_LAW`, read when the field members land |
| module | planetgen.physics.shape_temperature | GEN.1201 | the temperature a body took its shape at (surface plus accretion plus the heat made in the conduction time, capped at the liquidus; lock-free pure part of GEN.1201); the wiring step passes it to the roundness engine in `body_strength` and stores it, once GEN.1159's columns are on main |
| row | 6.162:SHT1 | GEN.1201 | `shape_temperature.ICE_SPECIFIC_HEAT`, read when the wiring step lands |
| row | 6.162:SHT2 | GEN.1201 | `shape_temperature.WINDOW_SHAPE_FACTOR`, read when the wiring step lands |
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
| row | 6.116:RGP2 | GEN.1069 | `resonance_gaps.ORDER_WIDTH_EXPONENT`, read when the wiring step lands |
| row | 6.116:RGP3 | GEN.1069 | `resonance_gaps.RESONANCE_CLEARING_MYR`, read when the wiring step lands |
| row | 6.116:RGP4 | GEN.1069 | `resonance_gaps.MAX_ORDER, MAX_P`, read when the wiring step lands |
| row | 6.116:RGP5 | GEN.1069 | `resonance_gaps.FORCED_ECCENTRICITY_COEFFICIENT`, read when the wiring step lands |
| row | 6.116:RGP6 | GEN.1069 | `resonance_gaps.YARKOVSKY_COEFFICIENT_AU_KM_PER_MYR, YARKOVSKY_REFERENCE_AU`, read when the wiring step lands |
| row | 6.116:RGP7 | GEN.1069 | `resonance_gaps.ZONE_ECCENTRICITY_FLOOR`, read when the wiring step lands |
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
| row | 6.126:RN1 | GEN.1159 | `roundness.SPHERICITY_FACTOR`, not read until GEN.1159 wires the module |
| row | 6.126:RN2 | GEN.1159 | `roundness.ICE_STRENGTH_PA, ROCK_STRENGTH_PA`, not read until GEN.1159 wires the module |
| row | 6.126:RN4 | GEN.1159 | `roundness.METAL_STRENGTH_PA`, not read until GEN.1159 wires the module |
| module | planetgen.physics.ice_types | GEN.1159 | the strength of the other ices next to ice Ih and the ice type of a body from its components (lock-free); `roundness.py` takes the ratio as a factor on its ice strength when the types are wired |
| function | planetgen.physics.body_strength.ratio_for | GEN.1069 | the roundness ratio of a body from its mass and ice fraction; `body_status.roundness_ratio` and `dwarf_population` call it, and those modules are registered under GEN.1069 and GEN.1092 |
| function | planetgen.physics.body_strength.round_mass_for | GEN.1069 | the mass at which a body of this ice fraction turns round; `body_status.round_mass_kg` and `dwarf_population.round_mass_kg` call it, and those modules are registered under GEN.1069 and GEN.1092 |
| function | planetgen.physics.roundness.is_round | GEN.1159 | whether a stored roundness ratio is round (1 or more); the wiring step that stores the ratio and the pages that read it call it |
| row | 6.126:RN11 | GEN.1201 | `roundness.ROCK_CREEP_COLUMNS`, not read until GEN.1201 gives each body a shape-set temperature to pass to the engine |
| data | planet_class_fit.json | GEN.1207 | the percent-match table of every class; the shadow scorer reads it (GEN.1206) and the best-fit classifier makes it the way a stored letter is decided |
| module | planetgen.physics.class_fit | GEN.1207 | the percent-match scorer; the shadow scorer in src/tests reads it, the classifier replaces the first-match tree with it |
| module | planetgen.physics.planet_type | GEN.1219 | the terrestrial, transitional or gas type from the envelope share of the mass; the shadow report in src/tests reads it, GEN.1219 stores the type with GEN.1207's epoch |
| row | 6.156:PT1 | GEN.1219 | read by planet_type.type_of; the stored type waits for GEN.1207's epoch |
| column | binary_epoch_unix | by design | the daily orbit update (`store.advance_*`) reads and writes the pair's clock; nothing outside the store needs it |
| column | binary_next_update_due | by design | the daily orbit update (`store.advance_*`) schedules the pair's next step from it; nothing outside the store needs it |
| column | cluster_id | GEN.359 | `stars.cluster_id`, the cluster a star belongs to; nothing writes it until the sector fill links cluster stars to their cluster |
| column | dynamical_state | GEN.291 | `planets` and `moons`: every body is `bound`; the unbind and capture test of GEN.291 is what sets `unbound` |
| column | formed_at | GEN.346 | `asteroid_fields.formed_at`, the epoch an event-made field formed; written when a collision makes a field |
| column | host_node_id | GEN.332 | `planets.host_node_id`, the orbit node a planet orbits; NULL (its star) until circumbinary and S-type planets are placed at nodes |
| column | origin_event_id | GEN.346 | `asteroid_fields.origin_event_id`, the orbital event that made the field; written when a collision makes a field |
| function | planetgen.galaxy.cellkey.cell_keys | by design | the numpy twin of `cell_key`; the tests hold `sphere_runs` and the stored key to it on many points at once, the generator keys one point at a time |
| function | planetgen.galaxy.cellkey.sql_cell_key | by design | the text of the MySQL generated column's expression; `schema.sql` holds it written out and the schema test compares the two, so the column cannot drift from the Python key |
| function | planetgen.galaxy.potential.vertical_force | by design | the disk's vertical pull `K_z` at a radius and height, the check of the potential against the Oort limit and the scale height of a disk; read by tests |
| function | planetgen.galaxy.density.component_masses | by design | each density component's closed-form mass (GEN.379), the total `test_density_components` holds the numerical integral (`integrated_masses`) to; the potential (GEN.338) reads the component masses from `density_components` directly |
| function | planetgen.galaxy.density.predicted_star_count | by design | the sector gate of the density design (`predicted_star_count >= 1`) as one pure function of the position; `run_plan` and the sector profile apply the same product with the clamp at 0 inline, and the tests use this as the rule's definition |
| function | planetgen.galaxy.object_uid.is_on_demand | GEN.1161 | a belt or field member above the individually stored count is drawn on demand, never stored; the member ledger wires it |
| function | planetgen.galaxy.object_uid.captured_member | GEN.1083 | the designation and parent link of a belt member a planet has captured; the irregular-moon wiring prints it |
| function | planetgen.galaxy.object_uid.outside_reservoir_member | GEN.1083 | the designation of a captured body from a reservoir that is not a belt; the irregular-moon wiring prints it |
| function | planetgen.galaxy.object_uid.parse_captured | GEN.1083 | reads a captured body's designation and link back; the irregular-moon wiring and the designation cascade (GEN.1141) use it |
| function | planetgen.galaxy.region_run.next_units | GEN.1128 | the barrier of the region run (stages 0 to 2 for every unit, then the top-up); the command wires it when the top-up stage exists |
| function | planetgen.galaxy.spherical.cells_touching_sphere | GEN.105 | the sector search of the influence set (`influence.influence_set` calls it); it reaches the app with the acceleration kernel |
| function | planetgen.galaxy.spherical.shell_cells | by design | every cell of one shell in order, the independent reference the grid tests count, sum volumes and compare neighbours against; the generator reaches cells by band ranges and never walks a whole shell |
| function | planetgen.galaxy.star_formation.star_formation_rate_profile | by design | the present star-formation rate against radius that `test_star_formation_profile` compares with the observed profile; the generator reads `recent_factor`, never the profile |
| function | planetgen.galaxy.structure_draw.hot_component_bound | by design | the most the thick disc can weigh over the thin disc's (Comeron et al. 2021 Eq. 1), the bound `test_structure_draw` holds the drawn thick-to-thin ratio under; the draw does not read it |
| function | planetgen.galaxy.structure_draw.draw_neutron_star_unbound_share | by design | the kick share a seed draws, with the stellar mass pinned, so `test_natal_kick` and `test_halo_threshold` probe the halo threshold without building the whole structure; the structure carries the same value |
| function | planetgen.generation.bright_stars.cell_band_rate | by design | the stars a cell expects in a luminosity band, the independent form of the rate the block path of `backfill_cells` computes inline; `test_block_sampler` and `test_backfill_block_path` check the block path against it |
| function | planetgen.generation.possibility_table.concentration_from_spread | by design | the authoring helper that turns a published spread of a component into a Dirichlet row's concentration; a row written from a source that gives spreads calls it, the tests hold the formula |
| function | planetgen.generation.primary_mass_function.systems_per_star | by design | the systems born per star born (0.697 in the stellar mass function note) that `test_primary_mass_function` pins; generation reads `systems_per_solar_mass` |
| function | planetgen.generation.run_galaxy.backfill_bright_stars_around | by design | the bright-star backfill (GEN.23) as one call with its own progress bar and shell log (PERF.28); a run does the same draw in `scatter-backfill` units (`expand_backfill_bright_stars`), so this form is what the tests and `tests/sector_support.make_sector` call to backfill one sector |
| function | planetgen.generation.stats.job_shells_per_second | PERF.147 | the stored whole-job rate that starts the ETA of a task that takes over 15 seconds |
| function | planetgen.generation.tuning_harness.check_target | by design | the entry of the tuning harness: draws a shipped recipe at a target's subject and compares it; the tests hold every recipe to its target through it, the generator never runs it |
| function | planetgen.physics.belt_history.giant_event | GEN.1232 | the giant event and the companion limit of a belt, for the belt page |
| function | planetgen.physics.belt_history.history_text | GEN.1232 | the belt history sentences shown on the belt page |
| function | planetgen.physics.belt_history.is_polluted | GEN.1231 | the pollution flag of a white dwarf from the rate of its belt (GEN.1068 gives the belts that see the giants) |
| function | planetgen.physics.birth_radioisotopes.fraction_above | by design | the observed share of a fixed-stream sampler, read by test_birth_radioisotopes and the research note, never by the generator |
| function | planetgen.physics.brown_dwarf_types.type_probabilities | by design | the expected share of each class, read by test_brown_dwarf_types against the sampler it describes |
| function | planetgen.physics.clathrate_brine.hydrate_stable | GEN.471 | whether a clathrate hydrate is stable at a temperature and guest pressure; the surface, the moons and the methane supply read it when GEN.471 lands |
| function | planetgen.physics.clathrate_brine.dissociation_temperature_k | GEN.471 | the temperature above which a hydrate dissociates at a guest pressure; read by the surface, the moons and the methane supply when GEN.471 lands |
| function | planetgen.physics.clathrate_brine.eutectic_wt_percent | GEN.471 | the eutectic composition of a brine salt; read with the liquid range when GEN.471 lands |
| function | planetgen.physics.clathrate_brine.freezing_point_k | GEN.471 | the freezing point of a brine of a given strength; read by the surface and the moons when GEN.471 lands |
| function | planetgen.physics.clathrate_brine.liquid_range_k | GEN.471 | the temperatures at which some solution of a salt is liquid; read by the surface and the moons when GEN.471 lands |
| function | planetgen.physics.collision_outcome.conservation_error | by design | the mass and momentum check of a collision outcome; the tests and the census run it, the generator never needs it |
| function | planetgen.physics.comet_admission.observed_sampler | by design | the simulated observed sample of the comet admission selection model: the simulated side of the `comet_admission` selection model, run by test_comet_admission through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.comet_admission.expected_measured | by design | the expected count the comet admission model measures: the simulated side of the `comet_admission` selection model, run by test_comet_admission through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.compaction.taper_weight | GEN.1069 | the taper of the compaction law in log mass; `body_status.roundness` calls it, and that module is registered under GEN.1069 |
| function | planetgen.physics.elements.series | by design | a decay-chain lookup by parent, read by test_elements to check the table against the published chains |
| function | planetgen.physics.elements.isotopic_mass_fraction | by design | the isotope mass fraction the published specific activities are checked against, read by test_elements |
| function | planetgen.physics.field_kernel.accelerations | GEN.105 | the pull at an object's exact place; the update's motion vector reads the gravitational field there (the numba kernel `sector_path` loads) |
| function | planetgen.physics.giant_tilt.large_tilt_chance | by design | the expected share of ice giants tilted over 60 degrees that `MERGER_SHARE` gives, the check of the tilt draw; read by test_giant_tilt, never by the generator |
| function | planetgen.physics.giant_tilt.final_period_h | by design | the spin after the impact that tilts an ice giant, the check that the impacts that turn a Uranus also spin it as fast as Uranus; read by test_giant_tilt, the draw keeps the giant's own spin |
| function | planetgen.physics.habitability_world.harshness | POP.38 | one minus the worst core domain score, the 0 to 1 harshness the tech draw and the facility ratings read (population-and-politics.md); POP.38 makes the facility ratings read the host, which is where it is called |
| function | planetgen.physics.hosting_gate.moon_survives_tides | GEN.550 | the tidal-survival gate of a moon that would host a submoon; a moon hosts only through this gate and a sourced satellite-share row, which the recursive generation of GEN.550 adds |
| function | planetgen.physics.impact_spin.mean_rms_breakup_fraction | by design | the rms birth spin of a mass averaged over the distribution of the mass ratio, the expected value the draws are checked against; read by test_impact_spin and test_tidal_spin, never by the generator |
| function | planetgen.physics.jwst_selection.observed_sampler | by design | the simulated observed sample of the JWST target selection model: the simulated side of the `jwst_selection` selection model, run by test_jwst_selection through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.jwst_selection.expected_observed | by design | the expected count the JWST target selection observes: the simulated side of the `jwst_selection` selection model, run by test_jwst_selection through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.kepler_recovery.observed_sampler | by design | the simulated Kepler survey sample of the recovery model: the simulated side of the `kepler_recovery` selection model, run by test_kepler_recovery and test_kepler_multiplicity through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.kepler_recovery.expected_found | by design | the expected count the Kepler pipeline finds: the simulated side of the `kepler_recovery` selection model, run by test_kepler_recovery and test_kepler_multiplicity through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.king.mass_within_pc | GEN.359 | the share of a cluster's mass inside a radius; the sector fill that draws a cluster's stars reads it |
| function | planetgen.physics.king.density_per_pc3 | GEN.359 | the density at a radius in a cluster's own units; the sector fill that draws stars from a cluster's profile reads it |
| function | planetgen.physics.mass_radius_sample.pass_probability | by design | the chance a planet passes the mass and radius measurement cut: the simulated side of the `mass_radius_sample` selection model, run by test_mass_radius_sample and test_iron_branch through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.mass_radius_sample.observed_sampler | by design | the simulated mass and radius sample of the iron branch tuning: the simulated side of the `mass_radius_sample` selection model, run by test_mass_radius_sample and test_iron_branch through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.mass_radius_sample.population_pass_fraction | by design | the share of a population that passes the measurement cut: the simulated side of the `mass_radius_sample` selection model, run by test_mass_radius_sample and test_iron_branch through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.materials.class_composition | by design | the rows a planet class text names; read by test_materials_table, which fails when a class text of `tuning.PLANET_CLASSES` resolves to no row |
| function | planetgen.physics.moon_composition.origin_probabilities | GEN.1083 | the probabilities of a moon's origin with the captured row weighted by its mass (`captured_weight`, GEN.1081); the captured-moon wiring weights the origin draw with them |
| function | planetgen.physics.moon_orbits.drift_rate_m_s | by design | the instantaneous `da/dt` of one moon at a given spin; `evolve` integrates all moons together and does not call it, the tests use it to hold the Earth-Moon recession (3.8 cm a year), Phobos and Triton to the torque |
| function | planetgen.physics.neutron_star.radius_limits_km | by design | the radius bounds of the equation-of-state spread, read by the plausibility gate (phenomena_plausibility, test_phenomena_plausibility, its batch CLI) and test_neutron_star_radius |
| function | planetgen.physics.neutron_star_structure.integrated_maximum_mass_solar | by design | the integrated maximum mass of an equation of state, read by scripts/build_neutron_star_table.py (run by hand when the table is rebuilt) and test_neutron_star_structure_real |
| function | planetgen.physics.neutron_star_structure.integrated_curve | by design | the integrated stable branch, read by test_neutron_star_structure_real to check the shipped table against a fresh integration |
| function | planetgen.physics.neutron_star_structure.curve_with_crust | by design | the arrays scripts/build_neutron_star_table.py writes (run by hand when the table is rebuilt), never read at generation |
| function | planetgen.physics.orbits.classify_encounter | GEN.105 | the guard for very close passes (GEN.108); the orbital update calls it before it moves a body |
| function | planetgen.physics.pair_history.widest_forced_axis_au | by design | the widest present axis the history rule forces into a post-common-envelope pair, found by bisection; the tree draw asks `must_be_post_common_envelope` at the drawn axis, and `test_pair_history` pins the rule's limit with this |
| function | planetgen.physics.photometry.blackbody_colour | VIEW.11 | the VIEW.11 photometry library is the reader: a star page or map colour reads it when VIEW.11 is built |
| function | planetgen.physics.planet_architecture.survival_probability | by design | Pu and Wu's survival law, whose inverse (`survival_decades`) is what the build draws; `test_planet_architecture` pins the inversion against it |
| function | planetgen.physics.planet_architecture.ejection_yield | by design | the expected number and mass of ejected bodies per star, the yield set against microlensing; `test_architecture_losses` runs it, the generator never needs it |
| function | planetgen.physics.planet_survival.expected_passes | by design | the ionising passes N(a) whose root `planet_survival_limit_au` solves in closed form; `test_planet_survival` pins the root against it |
| function | planetgen.generation.belt_birth.is_inner | by design | the inner-belt test (a belt ending within 1.2 snow lines) the census of GEN.1156 counts with; the prevalence suites and the tests call it, generation does not need the word |
| function | planetgen.physics.belt_origin.conservation_error | by design | the check that a belt's formed mass equals what is left plus what the giants cleared and what was scattered; the tests hold every generated belt to it |
| function | planetgen.physics.belt_origin.destroyed_planet | GEN.1072 | the `Debris` of a planet the star's history destroyed; the engulfment and white dwarf belt events of GEN.1072 pass it to `belt_origin.build` |
| function | planetgen.physics.planet_zones.flux_at | GEN.129 | the light on a planet of a window's host from every star of the tree; the multi-star planet climate reads it (design 5.6, 6.5) |
| function | planetgen.physics.proton_range.csda_range_g_cm2 | GEN.444 | the proton range table of the shielding of GEN.444 (the data file is read through this function) |
| function | planetgen.physics.rprocess_enrichment.effective_duration_gy | by design | the inverse of constant_input_factor by bisection, read by test_rprocess_enrichment to state a history as a duration |
| function | planetgen.physics.sector_path.box_cell | by design | a box domain for the path integrator; the tests integrate through it because its edges are exact, the generator's only domain is `sector_inside` |
| function | planetgen.physics.sector_path.ball_cell | by design | a ball domain for the path integrator; the tests integrate through it because its edge is exact, the generator's only domain is `sector_inside` |
| function | planetgen.physics.sector_path.relevant_masses | by design | the selection step of `integrate_path` run alone, so the tests and the timing in math-and-library-efficiency.md can check and time it; the integrator calls the kernel's selection itself |
| function | planetgen.physics.small_body_composition.makeup_from_components | GEN.473 | the makeup of a body from rows of the materials table; `belt_member`, `field_composition` and `captured_composition` call it, and their wiring lands under GEN.473, GEN.1103 and GEN.1083 |
| function | planetgen.physics.small_body_composition.draw_ring | GEN.1055 | the composition of ring material from where it formed in its parent's disk (ice fraction by distance); the stored ring ice share and its darkening by pollution read it when GEN.1055 lands |
| function | planetgen.physics.solar_particles.events_per_year_above | by design | the expected event rate of the NOAA list, read by test_solar_particles against the sampled record |
| function | planetgen.physics.state_vectors.equinoctial_from_state | GEN.105 | the modified equinoctial elements (Walker, Ireland and Owens 1985) of the orbital update for circular and equatorial orbits |
| function | planetgen.physics.state_vectors.state_from_equinoctial | GEN.105 | the inverse of `equinoctial_from_state`, the update's propagated state |
| function | planetgen.physics.state_vectors.mean_anomaly_from_true | GEN.105 | the mean anomaly of an orbit the update reads from a state vector |
| function | planetgen.physics.stellar_abundances.branch_probabilities | by design | the branch table probabilities for a test or a report, read by test_stellar_abundances |
| function | planetgen.physics.stellar_spectrum.band_flux_w_m2 | by design | the band sum of a spectrum, read by test_stellar_spectrum and test_uv_shares_wiring to check the UV shares |
| function | planetgen.physics.surface_dose.share_past_cutoffs | GEN.444 | the shielding of the surface dose from the column and its composition (the spectrum past the cutoffs); radiation.py still uses the single attenuation length |
| function | planetgen.physics.surface_dose.effective_air_column_g_cm2 | GEN.444 | the effective stopping column of a mixed air, for the shielding of GEN.444 |
| function | planetgen.physics.surface_dose.neutron_weight | GEN.444 | the ICRP 103 neutron weight, for the secondary dose of GEN.444 |
| function | planetgen.physics.surface_liquids.stable_liquids | GEN.1218 | the second liquid and the mixed seas of the composition model (Raoult sea of an air, its density, the liquids stable at a state); `surface_liquid_world` picks one liquid with `dominant` today |
| function | planetgen.physics.surface_liquids.ideal_sea | GEN.1218 | the second liquid and the mixed seas of the composition model (Raoult sea of an air, its density, the liquids stable at a state); `surface_liquid_world` picks one liquid with `dominant` today |
| function | planetgen.physics.surface_liquids.sea_density_kg_m3 | GEN.1218 | the second liquid and the mixed seas of the composition model (Raoult sea of an air, its density, the liquids stable at a state); `surface_liquid_world` picks one liquid with `dominant` today |
| function | planetgen.physics.surface_liquids.stable_liquids | TEST.520 | only tests call it; wire it into the code that needs it or delete it (TEST.520 first sweep) |
| function | planetgen.physics.surface_liquids.ideal_sea | TEST.520 | only tests call it; wire it into the code that needs it or delete it (TEST.520 first sweep) |
| function | planetgen.physics.surface_liquids.sea_density_kg_m3 | TEST.520 | only tests call it; wire it into the code that needs it or delete it (TEST.520 first sweep) |
| function | planetgen.physics.unit_registry.constant_in | by design | the number of a physical constant in another unit system (G in AU, year and solar mass); `relations.py` says a law in those units takes its G from it, and the Gaussian-constant test holds it |
| function | planetgen.physics.unit_registry.units_of | by design | the registry's listing of a kind's units in declaration order; the registry audits and the page-side export tests walk it, the app looks units up by key |
| function | planetgen.physics.unit_registry.parts_dimension | by design | the dimension a compound unit's parts give it; the registry's self-check holds every compound's stated dimension to it |
| function | planetgen.physics.unit_registry.parts_factor | by design | the factor a compound unit's parts give it; the registry's self-check holds every compound's stated factor to it |
| function | planetgen.physics.unit_registry.from_base | by design | the inverse of `to_base`; the round-trip property test holds the two together, the app converts unit to unit with `convert` and `converter` |
| function | planetgen.physics.wd_pollution.observed_sampler | by design | the simulated polluted white dwarf sample: the simulated side of the `wd_pollution` selection model, run by test_wd_pollution and test_tuning_harness through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.wd_pollution.class_probabilities | by design | the detection chance by class of a polluted white dwarf: the simulated side of the `wd_pollution` selection model, run by test_wd_pollution and test_tuning_harness through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.wd_pollution.expected_sample | by design | the expected polluted sample of the detection model: the simulated side of the `wd_pollution` selection model, run by test_wd_pollution and test_tuning_harness through `tuning_harness.SELECTION_MODEL_MODULES`; the generator never needs an observed sample (TEST.527) |
| function | planetgen.physics.magnetism.induced_moment | TEST.520 | only tests call it; wire it into the code that needs it or delete it (TEST.520 first sweep) |
| function | planetgen.physics.radiation_belts.dose_at_orbit_gy_yr | TEST.520 | only tests call it; wire it into the code that needs it or delete it (TEST.520 first sweep) |
| function | planetgen.physics.radiation_belts.equatorial_field_t | TEST.520 | only tests call it; wire it into the code that needs it or delete it (TEST.520 first sweep) |
| function | planetgen.physics.photometry.bolometric_correction_v | VIEW.11 | the VIEW.11 photometry library is the reader: apparent magnitude on a star page reads it when VIEW.11 is built |
| function | planetgen.physics.species_thermo.log10_k | by design | the reference the equilibrium tests check the gas-equilibrium code against; only tests run it (TEST.524) |
| function | planetgen.physics.clearance_ledger.log_uniform_quantile | GEN.1068 | the belt placer draws a belt's distance log-uniformly over the spans the ledger leaves open; read when belt birth (GEN.1068) samples from the ledger |
