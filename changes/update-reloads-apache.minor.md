### Changed
- `update.sh` now reloads Apache itself when it is running (and restarts it when the update just enabled a module), instead of printing the command; if Apache isn't running or the command fails it prints the command as before (OPS.8). macOS and Windows are unchanged.
