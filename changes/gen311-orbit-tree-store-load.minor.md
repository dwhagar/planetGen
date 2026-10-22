**Headline:** A built star system now keeps its tree of orbits in the database, and the system scene places its stars from that tree.

### Added
- **The orbit tree is stored and loaded with the system (GEN.311).** `insert_star_system` writes the system's tree into `orbit_nodes`: a generated binary is one node made from the pair it was generated with (separation, eccentricity, plane, phase held at the database's orbit epoch, period, mass), and a system that carries its own tree has it stored as it is. `load_star_system` and the JSON form (`StarSystem.to_dict`, serialization version 4) bring it back as `StarSystem.orbit_nodes`. Each star row says where it hangs: `parent_node_id` and its offset from the system's centre of mass. The old binary columns stay filled and equal the node.
- **The Roche classification is stored on the pair (GEN.269).** The node of a built pair carries its `interaction` (`detached`, `semi-detached`, `contact`, `accreting` or `post-common-envelope`), worked out from the stars' radii at the pair's closest approach by the same classifier the member-tree generator uses (`multistar.classify_pair`).

### Changed
- `GET /api/systems/<id>/scene` carries the system's `tree`, and `positions_at` places every star from it with the same function that places the members of a scatter row. A circular pair lands where the stored phase put it; a wide pair now follows its own eccentric Kepler orbit (the browser's copy of the star orbits still draws it as a circle, until the scene's JavaScript moves to the tree).

generation-output: changed
