**Headline:** A star system keeps its stars and what orbits them in one list of hosts, so the code that visits every orbiting body no longer rebuilds the same four lists.

### Changed

- `StarSystem.hosts` is the one list of `Host(star, role, orbiters)` (primary, wide secondary, members, nodes); `system.primary`, `system.secondary` and `system.hosts_of(role)` read it. The saved keys (`planets`, `secondary_planets`, `member_planets`, `node_planets`) are unchanged, and so is every generated value (GEN.1664).

generation-output: unchanged
