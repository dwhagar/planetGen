**Headline:** The multi-star light of a planet has one implementation: the unused `planet_zones.flux_at`, which repeated `radiation_field`, is gone.

### Changed

- GEN.1559 (wiring burn-down, found by `scripts/wiring_guard.py`): `planet_zones.flux_at` was built for the multi-star planet climate and read by nothing; `radiation_field` (GEN.1441) already sums the light of the host and of every star outside it from `planet_zones.companion_orbits`. `host_members` stays, now marked by design in the wiring register (the radiation field's host choice is tested against it). The function, its tests, its wiring-register line and the mentions in the design notes and registry row RF2 are removed or point at `radiation_field`.

generation-output: unchanged
generation-affects: none
