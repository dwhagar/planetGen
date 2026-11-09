**Headline:** The atmosphere sweep's cloud rule counts the sulfuric acid cloud that the stored cloud depth and cover already include.

### Changed

- `tests/audit/sweep_atmosphere.py` rule `aerosol-cloud` compared the stored `tau_cloud` and `cloud_cover` with the condensate cloud alone. Since `draw_aerosols` adds the acid cloud (GEN.1019, `sulfuric_cloud.of_draw`: its depth adds to the condensate cloud's and its cover is the larger of the two), a hot wet-carbon-dioxide world with no condensate (sweep body `F5V/Nh9`, 1950 K) was reported as "no condensate but a cloud of depth 22.8". The rule now takes the acid cloud from the stored values (`sulfuric_cloud.of_body`), checks the condensate part of the depth against its limits and accepts the acid cloud's cover. Its breaker adds the acid depth so it still breaks the rule. The generator is unchanged.
- `test_every_body_of_the_sweep_obeys_every_rule` stays red on the `temperature` rule for eleven bodies until the GEN.1450 secular exchange refreshes the air of a planet whose eccentricity moves by any amount (a separate defect, handed to its owner).

generation-affects: none
