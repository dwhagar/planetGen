**Headline:** Every rocky planet and moon now has a stored ground: the minerals, ices and ocean that cover its surface.

### Added
- **The surface composition is drawn and stored (GEN.441, second slice).** After the hydrosphere, `surface_composition.generate_surface` gives each rocky planet and moon `ground_type`, `ground_grain` and the areal share of each ground type (ocean, ice, ultramafic, anorthositic, basaltic, felsic, clay, oxidised, metal, carbon, evaporite, sulfur), which add up to 1 (schema 159, migration 0159; NULL for a giant and for rows saved before). It reads the bulk fractions of GEN.428, the hydrosphere's cover, the surface temperature and the air: ices stand only where the temperature and the air let them, clay only where liquid water stood, oxidised ground only under free oxygen, felsic ground only with plate tectonics.

### Changed
- **The ground dose reads the stored ground.** The gamma and radon dose used a basaltic crust with a drawn felsic share (`RADIOGENIC_FELSIC_FRACTION_MAX`, `radiogenic.SURFACE_RESERVOIRS`, both gone); it now reads the area-weighted K, Th and U of the body's own ground (`surface_composition.land_rock_ppm`), so an ultramafic world has a lower ground dose than a basaltic one and a granitoid one a higher.
- generation-output: changed (generator epoch 108: the ground dose and the stream after the hydrosphere move for every rocky body).
