**Headline:** The admin pages lose the Compact checkbox and the Show on Galaxy Map button on generate jobs.

### Changed
- **Compact checkbox removed from the admin shell.** The compact layout is now always on; the switch, its script and its remembered browser setting are gone.
- **"Show on Galaxy Map" removed from generate job pages.** The button, the `made_url` field, `GET /api/galaxy/made`, the Galaxy page's `?made=` view and the map's ring layer for a run's sectors are removed, since the map cannot be loaded while a run is in progress.
