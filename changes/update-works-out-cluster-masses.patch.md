**Headline:** The update works out every cluster's mean system mass once, so no web request ever has to, and the Galaxy Map's first load after a restart is not slow.

### Changed
- `update.sh` runs the new `planetgen.cli.warm_shape` after the migration (step 5). It reads the distinct ages of the stored clusters and satellites of every galaxy database, works out `system_mass.mean_system_mass_at_age_solar` for each (a minute or so for tens of thousands, nothing to stop it) and writes them to `system_mass_means.json` in the data cache folder, world-readable, for the release that wrote it. A web process looks the ages up there. A missing, unreadable or other-release file is ignored and the means are worked out as asked, as before; a failure of the step only warns.
