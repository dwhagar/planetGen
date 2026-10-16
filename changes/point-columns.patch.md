**Headline:** The columns that hold an object's position, mass and movement are defined in one place.

### Changed
- **One definition of the stored point (GEN.335).** `planetgen/db/point_columns.py` names, per table, the position, velocity, mass, epoch and due columns with their frame and units, and converts a row to a `SpatialPosition3D` and back. The planets, moons, bright stars, scatter points and star systems read and write their point through it (the nightly orbit update, the sector load, the path and stone lookups), with the same behaviour as before.
- **A test lists every table that holds an object.** A table must be a full layout (a point, a mass, a movement vector and an epoch), be on the named allow-list of tables still to convert (stars, black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae, remnants, belts, the quasar, comets, facilities and scatter members, each with the item that fixes it), or be a non-object. A new kind of object that lacks a part now fails the suite.
- **A galaxy fuzz test no longer fails on rounding (TEST.192).** The skeleton-bound test subtracts the halo from the model density, and where the halo is nearly all of it the remainder is smaller than the rounding of the subtraction; the check now allows four ulps of the model.

TODO items resolved: GEN.335, TEST.192. Added: none.
