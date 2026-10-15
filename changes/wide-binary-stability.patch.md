### Changed
- Wide binaries keep both stars' orbit regions apart (GEN.228). Each star's stable limit is now the tighter of the Holman & Wiegert fit and its own Hill sphere in the pair at periapsis (`orbits.hill_radius_in_pair_au`), the clearance between the two stars' outermost planets is measured at the pair's closest approach, not its average separation, and every comet a wide binary's star gets is an elliptical one whose whole orbit stays inside that limit (re-rolled up to `COMET_CEILING_ATTEMPTS` times, dropped when none fits). Before, comets ignored the limit: with comets forced on, 1,266 of 1,500 generated wide binaries had a comet reaching past it, toward the companion's planets.
- `validation.check_star_system` reports a comet that leaves its star's stable limit, and `stabilize_star_system` with removal drops it.
- `planetgen check-db --deep --wide-binaries` validates only the wide binaries, to count the stored pairs that break the rule.

TODO items resolved: GEN.228 (the count of existing broken pairs is `check-db --deep --wide-binaries`; repair waits for Boss) / Added: none
