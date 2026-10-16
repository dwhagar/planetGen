**Headline:** The galaxy settings files are gone; the galaxy keeps one internal record of how it was made.

### Changed
- `planetgen plan` no longer writes a JSON settings file, and the Stats page and `GET /api/admin/galaxy-settings[/<name>]` no longer offer one. The mismatch warning compares the running release, Python, platform and version key with the record that `galaxy_shape` already holds, and still names the sectors made by another version.
- `update.sh` deletes the old `galaxy-settings` folder from the jobs directory once.
