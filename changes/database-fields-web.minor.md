### Added

- **New database fields on the web pages.** A black hole's page shows its Class (stellar, intermediate or supermassive) and a rogue planet's its Mass Class. A runaway or hypervelocity star shows its speed as a badge on its system page and in its sector's Contents, and `GET /api/systems/<id>` returns `runaway_class` and `runaway_speed_kms`. The sector page adds its star count and the estimated number of interstellar comets and planetesimals drifting through it, and folds two or more rogue planets into one expandable Contents row.
