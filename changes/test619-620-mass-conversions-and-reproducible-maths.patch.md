**Headline:** The mass-conversion and reproducible-maths guards are green again.

### Fixed
- **Hand-written mass conversions (TEST.619).** `galaxy/keepout.py`, `galaxy/resident_mass.py`, `generation/rogue_topup_build.py` and `physics/remnant_composition.py` convert a mass through `unit_registry.converter`. The last three give the same bits as before (checked on 200,000 random values); the sphere of influence of a supermassive hole in `keepout.py` may differ in the last bit.
- **Maths that depends on the Python or BLAS library (TEST.620).** `db/field_clouds.py`, `db/obstacles.py` and `galaxy/unresolved_crossings.py` use `detmath.dist`, and `galaxy/avoid.py` and `galaxy/sky_band.py` use a new `detmath.dot` (a vector's dot products summed left to right, the bits of `detmath.matmul`) instead of `@` and `numpy.dot`.
