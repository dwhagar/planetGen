**Headline:** The admin page and the map rebuild can switch maintenance mode off again on a real server.

### Changed
- The maintenance switch (`maintenance.json`) is kept in `/var/lib/planetGen/state`, a folder `install.sh` and `update.sh` give to Apache's user (`examples/apache/create-cache-dir.sh`), instead of beside `settings.json` in the install folder, which is root's: turning maintenance off from the admin page, and the map rebuild's own switch-off after `update.sh`, failed there with "Permission denied", so the site stayed closed. `PLANETGEN_STATE_DIR` names another folder. The old `/var/lib/planetGen/maintenance.json` is no longer read, so a server that was closed by it opens on the next update.
- A switch that can't be written is no longer silent: the admin page says so with what to do (and shows a warning while the folder isn't the web user's), the command prints it and exits 1, and `planetgen.cli.warm_map` ends with an error saying the Galaxy Map is built but the site is still closed.

TODO items resolved: none (a bug found on Boss's server; filed with TODO) / Added: none
