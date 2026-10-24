**Headline:** The population model, the tuning table and the body-look pressure convert milliparsecs, megayears and pascals through the one unit registry, with the same numbers.

### Changed
- `population/model.py` (system and capital positions from milliparsecs to parsecs), `tuning.py` (the open-cluster age windows from megayears to gigayears) and `web/maps/bodylook.py` (pressure from pascals to kilopascals) call `unit_registry.converter` instead of a typed `/ 1000.0` (GEN.571, fifth file group). Each replaced expression gives the same float bit for bit, so no generated value moves. The `web/maps/starmap.py` sites use `physics.units` helpers and wait for the epoch-changing pass.
