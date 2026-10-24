# Surface and crust composition (GEN.441)

Code: `physics/surface_composition.py`. Tests: `tests/test_surface_composition.py`. Research base: part E of the body composition research (`body-composition-parts/E-ground-radiogenic-shielding.md`, sections 2.2, 2.3, 2.6, 3.2 and 3.5) and the gated possibility-table rule (`composition-model.md` section 4b).

## What it reads and gives

Reads (a `SurfaceInputs`): class, surface temperature and pressure, the air's oxygen and other partial pressures, the bulk fractions of GEN.428 (water, other ices, rock, metal, rock FeO, carbon rock), the hydrosphere state, ocean fraction and class of GEN.461, the internal and tidal heat flow, the age and the plate-tectonics chance.

Gives (a `SurfaceComposition`): the areal `mix` (ocean, ice and the ground types, summing to 1), the mass `components` (rows of the shared materials table, summing to 1) and the grain state. A giant has no ground (`None`).

## The order of the draw

1. Covered share from the hydrosphere: liquid ocean (with `RANGES["ocean_salt"]` of salts) or ice. Only ices that stand at the surface temperature and air are kept (`ice_stable`: water by the repo's 273.15 K under air and 220 K in vacuum rule, other ices by their saturation curve against their partial pressure or the vacuum floor, which is water's vapour pressure at 220 K). A frozen body whose water and ices reach half the rock's mass is ice everywhere.
2. Dominant ground of the land from `surface-dominant-ground`: metal (metal share 0.7 or more), carbon (carbon rock), basaltic (heat flow over the drawn floor), else ultramafic or anorthositic by a draw whose chance rises with the rock's FeO.
3. Modifiers, in the fixed order felsic, clay, oxidised, evaporite, sulfur, each a uniform areal share of the land when its inputs allow it, scaled together to at most 0.9 of the land.
4. Grain: slab if molten (1400 K) or slab or crushed if under 0.1 Gyr, else powder.

The recipes (`surface-ground-types`) are Hu et al. 2012 Table 2 on materials rows. Stand-ins are named on the rows: kaolinite for Mg-smectite, troilite and kamacite for pyrite, `basalt_powder` for glass, `orthoclase` for K-feldspar.

## Parametrisations (all in `RANGES`, none a measurement)

The heat-flow floor 1e-4 to 1e-3 W/m2, the overturn scale 0.2 of FeO, the modifier share ranges, the 0.1 Gyr young-surface limit, the ice-rich ratio 0.5, the stripped-mantle share 0.7, the free-oxygen floor 20 Pa and the ocean salt share 0.035 (Earth's). Each stands for a fit against data and is on the part E gap list.

## Reads back

`band_albedo` gives the lab reflectance by band area-weighted over the ground types that have a measured spectrum (not a Bond albedo; the final albedo is GEN.462). `crust_reservoir` gives K, Th and U for ultramafic, basaltic and felsic ground from the reservoir models of part E table 2.8b, for the radiogenic heat of GEN.459 to replace its basalt and felsic stand-ins.

## Still to do (waits on GEN.428/429, GEN.459, GEN.461)

Call it from `generate_planet_properties` after the hydrosphere; store `ground_type`, `ground_grain` and the areal mix (columns and migration, riding the chain's epoch); the child table of component shares; the API and describer; the invariants of E 3.6 against stored values.
