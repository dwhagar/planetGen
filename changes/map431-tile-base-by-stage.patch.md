**Headline:** an epoch that did not touch the galaxy's plan or scatter no longer makes the map draw every tile again.

### Changed
- MAP.431: the map's tile cache key (`tile_base`) counts the epoch a stored shape or scatter was made under as the last epoch that changed that stage (`generator_epoch.last_change`, from each epoch's `affects`), not the raw epoch number. A stamp that moves the galaxy up an epoch that changed neither the plan nor the scatter (a refill of the derived values, a new fill) keeps every stored tile; an epoch that changed the scatter still draws the scatter's tiles again. The first load after this release draws the tiles once, as the key's make-up changed.

generation-output: unchanged (no golden, no epoch, no lock item)
