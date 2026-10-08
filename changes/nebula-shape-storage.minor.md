### Added
- Each nebula now keeps a seeded irregular shape (schema v56: `nebulae.shape_*` columns and a `nebula_shape_balls` table), and `GET /api/nebulae/<id>/shape` serves it as a triangle mesh at a low-poly or full level of detail. Nebulae already saved get the same shape drawn from their own properties (GEN.75, part 2 of 3). **Run `sudo ./update.sh` to migrate.**
