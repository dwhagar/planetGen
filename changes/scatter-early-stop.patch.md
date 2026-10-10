### Changed
- The layer-walking scatter passes (stars and phenomena) stop and go on to the next stage once 100 layers in a row, in walk order, produced nothing (`tuning.SCATTER_DRY_LAYERS`, 0 walks every layer). The stop is logged and stored as a stage result (`layers_walked`, `stopped_early`); the nucleus and hypervelocity phenomena are not layers and are never skipped.
