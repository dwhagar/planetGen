**Headline:** The equipment tier a world stores now matches what its description asks for: an environment suit where it is too hot, too cold or too rich in oxygen, not a mask or nothing.

### Changed
- `habitability.equipment` names an environment suit when the temperature is past the Green heat or cold edge or the inspired oxygen is above the chronic limit, as `protection` and the world's page already say. Such worlds' stored `equipment_tier` rises (and so do the search page's equipment tags).
- `tests/test_wording_agrees_with_the_index.py`: the gear class and the equipment label now agree for every world, so the expected failure is gone.

generation-output: changed

TODO items resolved: GEN.561 / Added: none
