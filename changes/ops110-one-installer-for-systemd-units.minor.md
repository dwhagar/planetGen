**Headline:** `install.sh` and `update.sh` now install every maintenance timer and the Update button's units with the real checkout and virtual environment paths, and `install-maintenance-timer.sh` is a thin caller of the same code.

### Changed

- `scripts/deploy-common.sh` (OPS.110): one function, `install_maintenance_units`, renders the monthly orbit update, the unattended update timer, the Update button's path unit, the nightly tile rebuild and (with `--with-smoke-timer`) the smoke run for the checkout the script sits in and the virtual environment (`PLANETGEN_VENV`, default `/usr/local/planetgen/venv`), writes them to `/etc/systemd/system`, reloads systemd and enables the timers. It replaces `ensure_update_request_units`. A path with a space or a special character is refused and nothing is written. A unit already installed that differs from the rendering is replaced and the log says so; a deploy enables only a unit it installs for the first time, so a timer an operator disabled stays disabled.
- `install.sh`, `update.sh`: call it on every run on a systemd host. A failure only warns.
- `examples/maintenance/install-maintenance-timer.sh`: calls the same function (no copy of the rendering), keeps its options, and still writes the credentials file; it enables every timer it was asked for.
- The gunicorn unit stays a manual copy. `examples/systemd/planetgen-gunicorn.service`, `docs/deployment/nginx.md` and `docs/deployment/caddy.md` now say which two paths to edit: the checkout and the virtual environment.
- `INSTALL.md`, `docs/reference/database-schema.md` and the unit file headers describe the automatic install.
- `src/tests/test_systemd_units_installer.py` (new) renders every unit under a temporary checkout and finds no `/var/lib/planetGen` or `/usr/local/planetgen` in the installed copies.

generation-output: unchanged
