**Headline:** A backed-up sector repair file can be put back and is checked against the database straight away.

### Changed

- `python3 -m planetgen.cli.backup restore --part parity ARCHIVE` (and any restore that names the `parity` part) puts the backed-up repair file back (DB.41). The copy is checked first (whole records only), the galaxy database is reached before anything is changed, and the restore is refused at once, before the safety archive is taken, while a save is writing parity. The file is written beside the live one, synced, and renamed over it with the live file locked.
- The restored file is then checked against the database (`db/parity.reconcile`, the check of DB.8 with the stale-record rule of DB.9): a record that no longer matches the rows because the database changed after the backup is made again from the rows and is not read as damage; a record that fails its CRC is made again when the sector verifies; a member whose sector is gone is let go. The result says how many sectors the copy holds no record of (unprotected until their next save or `planetgen refresh-parity`) and lists the sectors that are damaged or cannot be trusted, for `planetgen repair-db`.
- A backup part may now refuse a restore before anything is changed (`Part.check_restorable`), and `RestoreContext.database` names the galaxy database a part checks against; `ParityFile.locked(wait=False)` raises `ParityFileBusy` instead of waiting.

generation-output: unchanged
