**Headline:** the GEN.1554 zone-edge sweep script is back in scripts/ and works on the current generator.

### Added
- **`scripts/gen1554_zone_edge_sweep.py`.** It draws seeded M5V, M2V, K5V, K0V and G2V systems and tables their terrestrial worlds by incident flux and by position in the habitable zone. The earlier copy was not in the repository and read `system.planets`, which a system no longer has, so a run drew every system and kept nothing. It now lists the primary star's planets through the host's orbiters, prints progress every 10 systems and keeps a partial file. The rerun on main is in `reports/gen1554-rerun-2026-11-27/`.

generation-affects: none
generation-output: unchanged
