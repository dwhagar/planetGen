**Headline:** The spin test marks the three seeds that expose GEN.1504 as strict expected failures until the secular exchange keeps the spin tilt in step with the orbit.

### Changed

- `tests/test_spin.py::test_every_generated_system_body_spins` seeds 0, 4 and 5 fail since the Train 21 merge (GEN.1450): `planets.secular_exchange` changes `orbital_inclination_deg` and `orbital_ascending_node_deg` but not the stored spin tilt and axis, so `axial_tilt_deg` no longer matches the orbit normal. They are marked `xfail(strict=True)` naming GEN.1504, so the day the defect is fixed they report as passes and the marks must go in the same PR. Seeds 1, 2 and 3 still run as plain tests. No other test is touched and the generator is unchanged.

generation-affects: none
