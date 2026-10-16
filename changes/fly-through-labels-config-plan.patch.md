**Headline:** The fly-through camera pieces that MAP.146 still lacked, map labels and the single configuration file are now planned.

### Added

- MAP.225 to MAP.228: the wheel zooming toward the cursor, the double-click flight, the container named from the camera position and the camera in the address, each a subitem of MAP.146 in 8.1 phase 2.
- MAP.229 to MAP.233: a map label engine (priority, no overlaps, show and hide by zoom), names for stars, systems and sectors, the Sector and System Maps on the same engine, names for the big features (8.8) and for nebulae (8.8).
- ADM.96: one configuration file, the user's own `config.json`, with no `settings.json` overlay; the example is never written and updates never overwrite the user's copy.

### Changed

- ADM.43, ADM.44, ADM.83 to ADM.85, ADM.87, USR.3, GEN.303 and the settings design note now say one file: secrets can be replaced but never read back, archive downloads are Owner-only, a web restore puts back only admin-changeable options, and the writer is `planetgen.util.settings_writer`.

TODO items resolved: none / Added: MAP.225, MAP.226, MAP.227, MAP.228, MAP.229, MAP.230, MAP.231, MAP.232, MAP.233, ADM.96
