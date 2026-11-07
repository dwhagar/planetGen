**Headline:** The update script no longer warns about files that already exist.

### Changed
- **`update.sh`'s pull step stopped listing untracked files as "uncommitted local changes that will be overwritten".** Untracked files (`config.json`, logs, caches) always exist and the reset never touches them, so only edits to tracked files are warned about now.
- **`set-permissions.sh` prints a plain note, not a warning, when `config.json` exists in both the checkout and the state folder.** Both files are left alone, as before.
