**Headline:** An update.sh that was already running when it pulled the web-user ownership change says to run it again instead of only refusing.

### Fixed
- `set-permissions.sh` now says, when it is handed `src/html` (what the update.sh that pulled the OPS.52 change still passes), that it takes the checkout and that `./update.sh` should be run again. The second run is complete.
- `update.sh`'s warning about uncommitted local changes ignores file-mode changes, which `set-permissions.sh` makes on every script.

TODO items resolved: none / Added: none
