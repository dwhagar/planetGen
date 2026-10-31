**Headline:** The installer and updater now tell you how big MySQL's redo log should be, from the writes the server has actually seen.

### Changed

- `scripts/redo_log_advice.py` reads how much redo MySQL has written and its redo log capacity, and prints the capacity that holds 30 minutes of those writes (1 to 16 GB) with the config line and the `SET GLOBAL innodb_redo_log_capacity` that changes it without a restart on MySQL 8.0.30 and later. `install.sh` and `update.sh` run it after the Redis check and only print; with `--seconds N` it measures a live window, and it counts the checkpointer warnings (MY-014084) in the error log. INSTALL.md says how to read it. Nothing on the server is changed. generation-output: unchanged
