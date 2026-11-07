**Headline:** The server side of the admin Update button: update.sh runs as root from a request file, with no root shell.

### Changed
- **SEC.43: root-owned systemd units for the admin pages' Update button.** The web user only creates an empty `update.request` file in the state folder. New `planetgen-update.path` sees it and starts the existing `planetgen-update.service`, which removes the request and then runs `update.sh` as root. The file is never read and the command takes no argument, so nothing the web user can write is run as root.
- **A press while an update runs waits and runs after it.** The update guard's own rules still decide whether running work is cancelled or the update refuses.
- **`install.sh` and `update.sh` keep the units right** (`ensure_update_request_units`: the checkout path and state folder put in, the path unit enabled; Linux with systemd only, a failure only warns). `install-maintenance-timer.sh` installs the path unit too, and `--skip-update-button` leaves it out.
- The admin page, route and status view are Core 3's part of SEC.43. Tests: `test_update_button_units.py`.
