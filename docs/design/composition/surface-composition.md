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

## Stored (second slice, GEN.441)

`generate_surface(body)` runs right after the hydrosphere and sets `ground_type`, `ground_grain` and `ground_<type>_fraction` (twelve areal shares that add up to 1; schema 159). A giant's fields are `None`. The non-water ices' proportions are a stand-in (the median comet volatiles, `comet_population.median_volatiles`) until the species are stored with the bulk composition; the tidal heat flow is not yet an input, so the sulfur modifier waits on it.

The ground dose reads the stored ground: `land_rock_ppm(body)` is the area-weighted K, Th and U of the land (`RESERVOIR_OF`: clay and iron-oxidised ground follow basalt, being altered basalt; land with no reservoir takes basalt's), and `radiation._ground_dose` raises the bulk silicate's content to it. `radiogenic.SURFACE_RESERVOIRS` and `tuning.RADIOGENIC_FELSIC_FRACTION_MAX` are gone.

## Rock types from the oxides (GEN.474, lock-free step)

`GROUND_TYPES` is a library: each ground type is a fixed mineral list taken from the Earth's and the Moon's crusts (Hu, Ehlmann and Seager 2012 Table 2). `physics/rock_types.py` (registry RT1 to RT3, section 6.45) derives the ground type from the oxides of the body's own rock instead:

- **The norm.** The CIPW norm allocates the oxides to ideal minerals by mass balance in the method's fixed order (apatite, ilmenite, orthoclase, albite, anorthite, magnetite, diopside, hypersthene; quartz from the silica left, olivine, nepheline, leucite and periclase when it falls short). The norm adds up to the oxides given (tested); a peralkaline rock or one without the silica for any norm is refused.
- **What melts.** A mantle becomes a crust by partial melting. The minerals of the norm melt in three stages: the low-melting feldspars, quartz and diopside together in the proportions the rock holds them (an invariant eutectic, so the liquid keeps one composition while the stage lasts: a basalt for the bulk silicate Earth up to a fifth melted), then hypersthene, then olivine (a magnesian liquid, komatiite-like, over a third melted). K, Th and U are incompatible, so a batch melt holds them `1 / (D + F (1 - D))` times the source (D taken as zero); N-MORB over depleted mantle (22.8, 31.9 and 25.9 times for K, Th, U) is one melt of 3 to 4.5 percent for all three, which the test pins as a consistency check and never uses as an input.
- **The ground type** is anorthositic at 90 percent plagioclase in the norm (IUGS), ultramafic under 45 percent silica or over 18 percent magnesia (komatiite), then basaltic (45 to 52), intermediate (52 to 63) or felsic (63 and over) by the silica bands of the total alkali-silica diagram. The upper continental crust (Rudnick and Gao 2014 Table 11) comes out felsic, the lower crust intermediate, a lunar highland anorthosite anorthositic and the bulk silicate Earth ultramafic, all as test points.

Limits: the staged melt is not a phase-equilibrium liquid (MELTS or pMELTS would give one: the stage mix is too aluminous and calcic for a MORB, 23 percent Al2O3 against about 15, and its high-degree liquids rise in silica as the orthopyroxene joins, so 25 to 30 percent melts are called intermediate); the norm's iron split is proportional; there are no hydrous, carbonate or sulfur rocks. The wiring step is the second step of GEN.474: the surface draw takes the primary and secondary crust from the stored mantle through the melt (the melt fraction from the body's heat flow, the volcanism floor of `RANGES`), a new `intermediate` ground type is added to the recipes, and the library of Hu's assemblages stays as the calibration. It changes the stored ground type of rocky bodies (generation fill stage) and the crust K, Th and U, so it rides a merge train.

## Still to do

Ground albedo with weathering (GEN.462's follow-up reads `band_albedo`), the child table of component shares for the describer (GEN.430), the API and the describer lines, tidal heating as an input, and the invariants of E 3.6 against the stored mass balance.
