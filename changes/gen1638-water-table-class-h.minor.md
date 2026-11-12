**Headline:** a world with a little liquid water at its surface now counts as having a water table under its land when the classes are judged, so class H (the desert world with some water) can be met; it could not be before.

### Changed

- `class_domains.derived`: `groundwater_present` is 1.0 for a world whose hydrosphere holds liquid (`habitability_world.LIQUID_WATER_STATES`: surface ocean, brine seep, hycean, ice-covered ocean) and covers part of the surface, as well as for a dry world with an aquifer. Before, only a dry, icy or vapor world could be given an aquifer, and at 275 to 315 K and 50 kPa or more any world with water has a liquid surface, so H was unreachable (0 in 40-trial probes at bulk water 1e-10 to 5e-5, reports/gen1585-class-h-reachability-2026-11-12.md).
- The stored `groundwater` column is not changed: it means the water a world holds only underground, which the habitability score and the water source read (changing it would score a sea world as an aquifer world). Only the class reading moves. `solvent` reads the same state tuple from one place.
- The class page's wording for the H condition now says a water table under the land of a world with a small sea.
- Changes the class of oxygen-air worlds with under a fifth of the surface covered by liquid water (they were F or G), so it is a lock item: the golden and epoch are the train's.

generation-output: changed
generation-affects: fill, population
