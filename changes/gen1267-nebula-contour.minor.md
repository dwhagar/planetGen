**Headline:** A nebula is the outline the map draws: a system is inside a nebula only where the galaxy's density holds it, and a hot-star nebula in the outer disk, which the map drew as nothing, now has a body.

### Changed

- GEN.1267: the generator and the stored-row refresh decide which nebula holds a system with the same density the map draws (`galaxy/nebula_density.py` `inside`: inside the nebula's sphere, and the density read for it at least the form level, 3 octaves, rounded to 1e-4). Before, a system was "inside" a dark cloud by the GEN.75 metaball shape, which overlapped the drawn outline by 14 percent (measured on 77 clouds). `inside_nebula_id` of generated systems changes, so a galaxy is made again at the epoch of this change. Field-cloud object IDs change too: a cloud is born in the sector holding the center of its sphere (the metaball centroid it used is gone).
- MAP.423: the map and the test read a hot-star nebula (emission, reflection, diffuse, planetary, remnant) with its density anchored on its own sphere (`YOUNG_ANCHOR_FLOOR`, GEN.725 row NC1). The model's young-star density alone is under the form level outside about 6 kpc (0.6 percent of points at the solar circle, none at 10 kpc), so a 40 pc emission nebula near the Sun drew 15 points of strength 6 of 255. It now draws a body with a full-strength core at every galactic radius.
- `nebula_volume.sample_group`'s density callable takes the owner's center and radius as a fourth argument (`anchors`); `density_at` takes `anchors=None`.
- The generator battery case `nebula-shape` is now `nebula-contour` (which of 50 points a toy-galaxy nebula holds, by each density).

### Removed

- The GEN.75 metaball shape: `galaxy/nebula_shape.py` (the gradient noise moved to `galaxy/gradient_noise.py`), `Nebula.get_shape`, `query.nebula_shape`, the route `GET /api/nebulae/<id>/shape` (the maps draw `.../volume`), `store._store_nebula_shape`, the twelve `nebulae.shape_*` columns and the table `nebula_shape_balls` (migration 0241, schema v241). No shim.

### Tests

- `tests/test_nebula_contour.py` (the contour at every galactic radius, the anchor, the map and the test agree, a dark cloud, rounding, a database never planned). `test_nebula_field.py`, `test_gen_cloud_at_generation.py`, `test_db_persistence.py`, `test_generator_battery.py` follow the new rule; the shape endpoint tests are gone with the route.

generation-output: changed
generation-affects: fill
