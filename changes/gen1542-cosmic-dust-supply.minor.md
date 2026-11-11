**Headline:** Whether a planet has the cosmic dust that the Venus-cloud hypothesis V8 needs is now worked out from its size and place, not left to the caller.

### Changed

- New `physics/cosmic_dust.py` (GEN.1542): the cosmic dust that reaches a body, kg a year, from the Earth's 28 tonnes a day scaled by cross-section and focusing (`ring_origin.focusing`, the same encounter speed as the comets), the star's distance to the power -1.8 and the star's mass to the power 0.5. It gives Venus 31 tonnes a day and Mars 1.7.
- `sulfur_cloud_set.fill` passes `carbon_dust` to the chain from `dust_supply`: the body's influx per area at least `dust_min_share_of_reference` (0.1) of Venus's. Before it was always none, so V8 could never appear.
- V8 stays a single-source hypothesis and is off by default (`INCLUDE_HYPOTHESES`), so no stored value changes with the defaults. On 86 systems (663 rocky bodies) with the hypotheses switched on, one body gets V8 with the gate and one with the gate forced open: the acid cloud is the rare gate, not the dust. The gate closes a rocky body beyond about 2.5 AU of a Sun-like star.
- `ring_origin.focusing` is the focusing that `meteoroid_flux_kg_m2_s` computed inside; same numbers. Finding: the ring law's Saturn-anchored flux (exponent 1.5) gives the Earth 0.7 tonnes a day, 40 times under the measured 28, so the inner system is anchored at the Earth here.
- Stand-ins CV11 to CV14 in section 6.229; CV8 updated. New reference [R5228].

generation-output: unchanged
generation-affects: none
