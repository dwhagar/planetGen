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
| module | planetgen.physics.neutron_star_structure | GEN.969 | the TOV integration and the reduced curves; the neutron star radius and crust will read them in place of the Miller et al. radius band |
| data | neutron_star_curves.json | GEN.969 | the mass-radius curves and crusts of the 21 equations of state, read by `physics/neutron_star_structure.py` only |
| dataset | stellar_spectra | GEN.971 | the ATLAS9 and PHOENIX-ACES spectra reduced to band shares; the star's UV shares and map colour will read them in place of the blackbody |
| module | planetgen.physics.stellar_library | GEN.971 | the band lookup and the corona's far ultraviolet; replaces the seeded UV-C share of `activity.generate_activity` once wired |
| row | 6.32:NS1 | GEN.969 | `neutron_star_structure.CRUST_CORE_DENSITY_FM3`, not read until the radius and crust are wired |
| row | 6.32:NS2 | GEN.969 | `neutron_star_structure.NEUTRON_DRIP_DENSITY_FM3`, not read until the radius and crust are wired |
| row | 6.33:SL2 | GEN.971 | `stellar_library.FUV_OVER_XUV`, not read until the UV shares are wired |
| row | 6.33:SL3 | GEN.971 | `stellar_library.FUV_IN_UV_C`, not read until the UV shares are wired |
| module | planetgen.physics.alpha_tracks | GEN.744 | the equivalent [Fe/H] of an alpha-enhanced star; the track lookup of a high-alpha star will read it |
| module | planetgen.physics.age_metallicity | GEN.744 | the disc's age-metallicity relation by branch; the star's [Fe/H] draw will follow its age |
| row | 6.34:AL1 | GEN.744 | `alpha_tracks.ALPHA_SHARE`, not read until the track lookup is wired |
| row | 6.34:AR1 | GEN.744 | `age_metallicity.THIN_BIRTH_FE_H`, not read until the [Fe/H] draw follows the age |
| row | 6.34:AR2 | GEN.744 | `age_metallicity.THICK_ANCHOR_FE_H`, not read until the [Fe/H] draw follows the age |
| row | 6.34:AR3 | GEN.744 | `age_metallicity.SIGMA`, not read until the [Fe/H] draw follows the age |
| module | planetgen.physics.stellar_state | GEN.786 | the state word of a star and the age since it stopped burning (GEN.785); `Star.__init__` fills `stars.state` and `stars.decay_age_yr` once Star 1's GEN.935 phases are in the star epoch |
| module | planetgen.physics.line_opacity | GEN.850 | the band opacities of the greenhouse model, called only by `radiative_convective` |
| module | planetgen.physics.radiative_convective | GEN.850 | the greenhouse column, replaces the `infrared_opacity.py` plateau code once wired |
| module | planetgen.physics.photolysis | GEN.982 | the photolysis lifetime table for the exosphere and the coma losses; no caller yet |
| data | photolysis_lifetimes.json | GEN.982 | read by `physics/photolysis.py` only |
| module | planetgen.physics.mixture_trace | GEN.524 | the trace gases of an air drawn after the life stage; no caller yet |
| module | planetgen.physics.design_flags | GEN.524 | gates the trace rows `mixture_trace` draws |
| data | design_flags.json | GEN.524 | read by `physics/design_flags.py` only |
| module | planetgen.physics.equilibrium_edges | GEN.485 | the edges of the hot-Jupiter chemistry regimes, sourced and wired with the hot-giant thermochemistry |
| module | planetgen.physics.ground_albedo | GEN.462 | the colour of the ground by band, wired by the ground albedo item |
| module | planetgen.physics.outside_hours | GEN.768 | the hours a person can be outside in a day, stored per planet by the item |
| module | planetgen.physics.surface_dose | GEN.444 | the ground-level dose kernels, wired by the surface radiation item |
| module | planetgen.physics.radiation_belts | GEN.983 | trapped radiation belts of a magnetised body, read when the stored magnetosphere reaches dose |
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
| module | planetgen.physics.condensate_shares | GEN.950 | Lodders 2003 condensates by temperature and the water share of the solids (lock-free); the wiring step reads the snow line temperature, the ice boost of the surface density and the water range of an icy body from it |
| row | 6.35:CS1 | GEN.950 | `condensate_shares.ICE_CONDENSATION_K`, read when the wiring step lands |
| row | 6.35:CS2 | GEN.950 | `condensate_shares.TABLE_11_MASS_PERCENT`, read when the wiring step lands |
| row | 6.35:CS3 | GEN.950 | `condensate_shares.ROCK_OXYGEN_SHARE`, read when the wiring step lands |
| module | planetgen.physics.planet_architecture | GEN.728 | the planet architecture (inner system, thinning, outer chain, gas, migration); `system.py` calls it at the wiring on merge train 4 |
| module | planetgen.physics.planet_occurrence | GEN.728 | the reader of the observed occurrence bands; the architecture and the census read it |
| module | planetgen.physics.kepler_multiplicity | GEN.728 | the Kepler transit multiplicity harness, read by the tuning tests until the wiring |
| data | planet_occurrence.json | GEN.728 | read by `physics/planet_occurrence.py` only |
| row | 6.36:AR1 | GEN.728 | `ARCH_DISC_SOLIDS_ENHANCEMENT, ARCH_OUTER_GROWTH_*`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR2 | GEN.728 | `ARCH_GAS_*`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR3 | GEN.728 | `INNER_INITIAL_COUNT_EXTRA, INNER_COUNT_SIGMA_LN`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR4 | GEN.728 | `INNER_PLANET_MASS_*`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR5 | GEN.728 | `INNER_NEIGHBOUR_SIGMA_DEX`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR6 | GEN.728 | `*_SEPARATION_*`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR7 | GEN.728 | `OLIGARCH_INNER_EDGE_*`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR8 | GEN.728 | `SCULPT_AGE_YR and the survival law`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR9 | GEN.728 | `GIANT_MIGRATION_*`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR10 | GEN.728 | `GIANT_INWARD_MIGRATION_SIGMA_LN, GIANT_CLEARING_SPACING`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR11 | GEN.728 | `PLANET_MIN_MASS_*`, read by `planet_architecture`, which `system.py` calls at the wiring |
| row | 6.36:AR12 | GEN.728 | `DISC_SOLIDS_SCATTER_DEX, INNER_SYSTEM_FRACTION_ANCHORS`, read by `planet_architecture`, which `system.py` calls at the wiring |
