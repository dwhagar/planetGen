# The wiring register: what is built or downloaded and not read yet (GEN.754)

Status: hand-kept, checked by `scripts/wiring_guard.py` and `src/tests/test_wiring_guard.py`. Every module, setup dataset, dataset reader, data file and stand-in registry row (GEN.725) that nothing in the app reaches is named here with the open item that will wire it, or with the words `by design` and the reason. The test fails when something unread is missing from this page, when a line names something that is read now (delete the line in the same change that wires it), and when a named item is not open in `docs/TODO.md`. No shim: an unread piece is wired, filed or deleted.

What counts as read: the app is every module outside `physics`, `generation`, `galaxy`, `datasets` and `tuning` (the CLI, API, web, admin, queue, population, backup and describe packages). A library module is read when an import chain from the app reaches it, a lazy import inside a function counts, and so does a module named by a string (`importlib.import_module("planetgen.x")`, the registry's `"planetgen.x:function"`). A dataset reader is read when a reachable module outside `datasets` uses it. A data file is read when a reachable module names it. A registry row is read when one of its `module.CONSTANT` names is used by a reachable module.

`python scripts/wiring_guard.py` prints what the guard finds; `--check` is the test.

| kind | name | item | reason |
|---|---|---|---|
| dataset | harris_clusters | GEN.297 | the Harris catalog of 157 globular clusters; the calibration of the cluster generator reads it |
| dataset | noaa_proton_events | GEN.747 | the 265 NOAA solar proton events; the solar particle event model fits its rate and slope to them |
| dataset | baraffe_isochrones | GEN.936 | the Baraffe 2015 isochrones; the 0.08 to 0.1 solar mass main sequence (GEN.936) and the brown dwarf radius by age (GEN.948) read them |
| dataset | hitran_cia | GEN.850 | the HITRAN collision-induced absorption is reduced at setup and read by `line_opacity`, which the greenhouse step does not call yet |
| dataset | ir_ktables | GEN.850 | the ExoMol correlated-k tables are reduced at setup and read by `line_opacity`, which the greenhouse step does not call yet |
| reader | neutron_star_matter.equations_of_state | GEN.968 | the 21 equations of state are read by a test only; the TOV integration reads them |
| module | planetgen.physics.line_opacity | GEN.850 | the band opacities of the greenhouse model, called only by `radiative_convective` |
| module | planetgen.physics.radiative_convective | GEN.850 | the greenhouse column, replaces the `infrared_opacity.py` plateau code once wired |
| module | planetgen.physics.photolysis | GEN.455 | the photolysis lifetime table for the exosphere and the coma losses; no caller yet |
| data | photolysis_lifetimes.json | GEN.455 | read by `physics/photolysis.py` only |
| module | planetgen.physics.mixture_trace | GEN.524 | the trace gases of an air drawn after the life stage; no caller yet |
| module | planetgen.physics.design_flags | GEN.524 | gates the trace rows `mixture_trace` draws |
| data | design_flags.json | GEN.524 | read by `physics/design_flags.py` only |
| module | planetgen.physics.equilibrium_edges | GEN.485 | the edges of the hot-Jupiter chemistry regimes, sourced and wired with the hot-giant thermochemistry |
| module | planetgen.physics.ground_albedo | GEN.462 | the colour of the ground by band, wired by the ground albedo item |
| module | planetgen.physics.outside_hours | GEN.768 | the hours a person can be outside in a day, stored per planet by the item |
| module | planetgen.physics.surface_dose | GEN.444 | the ground-level dose kernels, wired by the surface radiation item |
| module | planetgen.physics.radiation_belts | GEN.800 | trapped radiation belts of a magnetised body, read when the stored magnetosphere reaches dose |
| module | planetgen.galaxy.influence | GEN.342 | the influence set of an object, read by the acceleration kernel of the orbital update |
| module | planetgen.physics.closest_approach | GEN.197 | the pair closest-approach library, read by the orbital event writer |
| module | planetgen.physics.collision_outcome | GEN.197 | what a collision does, read by the orbital event writer |
| module | planetgen.physics.body_positions | MAP.253 | where every body of a system is at a time, read by the scene and the maps once they read the parent chain |
| module | planetgen.physics.light_travel | VIEW.1 | light-travel positions for the view from a planet |
| module | planetgen.physics.radius_valley | GEN.722 | the simulated Kepler survey of the radius valley; the generator reads the result once the follow-up settles the mechanism |
| module | planetgen.physics.iron_branch_tuning | by design | the simulated sample the iron branch table is tuned to; run by the tuning tests, not by generation |
| module | planetgen.generation.class_sample | by design | a sample that holds every class, a helper of the sheet audit tests |
| module | planetgen.generation.plausibility | by design | the anomaly finder engine; run by the plausibility tests and their command line wrapper in `src/tests` |
| module | planetgen.generation.phenomena_plausibility | by design | the exotic-phenomena counterpart of the anomaly finder, run by its tests |
| data | ngpps_formation_water.json | by design | the binned Bern synthesis shares that `test_water_delivery` and `test_atmosphere_recipes` hold the delivered water to |
| row | 6.7:AM18 | GEN.524 | `mixture_trace.ABUNDANCES`, not drawn until the trace gases are wired |
| row | 6.8:OH8 | GEN.768 | `protection_limits.GARMENTS` insulation, read by `outside_hours` |
| row | 6.8:OH9 | GEN.768 | `protection_limits.GARMENTS` cooling and heating power, read by `outside_hours` |
| row | 6.10:RV2 | GEN.722 | `radius_valley.CORE_MASS_SIGMA`, a width of the simulated survey, read by the follow-up |
| row | 6.10:RV3 | GEN.722 | `radius_valley.ENVELOPE_SIGMA`, a width of the simulated survey, read by the follow-up |
| row | 6.15:RB1 | GEN.800 | `radiation_belts.WIND_EDGE_FRACTION`, read when the belts are |
| row | 6.15:RB2 | GEN.800 | `radiation_belts.ROTATION_EDGE_FACTOR`, read when the belts are |
| row | 6.15:RB3 | GEN.800 | `radiation_belts.WIND_PEAK_FRACTION`, read when the belts are |
| row | 6.15:RB4 | GEN.800 | `radiation_belts.INNER_PEAK_FRACTION`, read when the belts are |
| row | 6.15:RB5 | GEN.800 | `radiation_belts.ATMOSPHERE_EDGE_SCALE_HEIGHTS`, read when the belts are |
| row | 6.15:RB6 | GEN.800 | `radiation_belts.EARTH_OUTER_PEAK_GY_YR`, read when the belts are |
| row | 6.15:RB7 | GEN.800 | `radiation_belts.OUTER_RISE_EXPONENT`, read when the belts are |
