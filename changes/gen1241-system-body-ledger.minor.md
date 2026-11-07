**Headline:** A system page now shows the mass of the whole system, where its center of mass lies against the primary star, and how far its asteroid belts pull that star (GEN.1241).

### Added

- `GET /api/systems/<id>/bodies`: the body ledger of a system built from the typed tables with no new storage (stars, planets, moons, rings, comets, belts, the bodies each belt keeps and the belt's remainder), its total mass, barycenter, the central star's offset, the belt wobble and each belt's rows against its stored mass.
- The "Mass and barycenter" panel on the system page, with the mass and distance of every body from the center of mass and a table of the belts (`web/system_ledger.py`).
