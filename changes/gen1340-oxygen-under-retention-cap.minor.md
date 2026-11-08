**Headline:** The oxygen life puts in the air no longer lifts the air above the pressure its body can keep, so thin airs are not 75 to 95 percent oxygen at 23 to 45 kPa (GEN.1340).

### Changed

- The air solve holds the air to the pressure its retention regime keeps (row `d6`, 1 kPa for a tenuous regime, 10 kPa for a moderate one), but the biotic oxygen was added afterwards, on top of that cap. In the census after train 14, four of the ten oxygen worlds with a living sea had 23 to 45 kPa of air that was 75 to 95 percent oxygen over 1 to 10 kPa of other gas.
- `surface_air.Air.retention_cap_pa` carries the cap out of the solve (kept on the planet like the surface water, not stored); `biogenic_oxygen.apply` scales the whole air to it when the oxygen would take the total over it, keeping the shares of the gases, as the pressure stage does. A body without a cap is as before.
- Measured: see the report `gen1340-oxygen-under-retention-cap-2026-11-07.md`. Every thin or moderate-regime world that makes oxygen changes (epoch).
