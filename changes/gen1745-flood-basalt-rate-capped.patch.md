**Headline:** A body that tidal heat has melted no longer asks for billions of flood basalts a gigayear, which ran a system's build out of memory.

### Fixed
- GEN.1745 slice 4b: `air_events.flood_basalt_rate_per_gyr` stops at `FLOOD_BASALT_MAX_PER_GYR` (1,000 a gigayear: one flood basalt takes about a million years to erupt, so a body whose heat asks for more is erupting all the time and not more often). A tidally heated rocky body (25 MW/m2) gave 6.5e10 flood basalts over its history, the draw never ended and `tests/test_aquifers.py` was killed at 12 GB on the train; the same system now builds in the memory it had before the air record.

generation-output: unchanged
generation-affects: none
