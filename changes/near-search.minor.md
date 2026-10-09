### Added
- `GET /api/near`, `planetgen.db.near.objects_within` and `python -m planetgen.cli.query near` (NAV.43): everything generated within up to 50 parsecs of a place (an object reference such as `system:12`, or a galaxy-frame point), of every kind: systems with their stars, planets, moons, belts and comets, facilities and the standalone phenomena. Nearest first, paged, with a kind filter and a count of the sectors in range that are not generated yet. The search never generates anything.

### Removed
- `GET /api/systems/<id>/near` and `queryDb.systems_within_radius` (same-sector systems only, in light-years): `GET /api/near` replaces them.
