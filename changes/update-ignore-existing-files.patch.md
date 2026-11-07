**Headline:** The update script no longer warns about files that already exist.

### Changed
- **`update.sh` stopped warning when `config.json` exists in both the checkout and the state folder.** Existing files are normal on every update after the first, so `set-permissions.sh` now prints a plain note and leaves both files alone.
