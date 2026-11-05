**Headline:** A full backup now holds a copy of the sector repair file, so a galaxy whose live repair file is lost can still be repaired.

### Changed

- New backup part `parity` (DB.40, `backup/parts.py`): a copy of the parity file of the database repair (`db/parity_file.py`, DB.9) in `parity/sector-parity.pgpr`. The copy is taken with the file locked against every save (the lock the saves take) and cut at the end of its last whole record, so it is a state the file really had, never half a record. A file damaged behind its last record is refused with the reason (run `planetgen repair-db`), because that copy could repair nothing. With no parity file yet the part is empty and says so.
- The manifest keeps what a part says about its dump as `info` (`Part.dump` may return a dict); for the parity part it is `bytes`, `sectors` covered and `groups`, or `present: false` with the reason. `planetgen.backup.PartInfo.info` holds it.
- New option `settings.backup_parity` (default on, Admin settings, "Back up the repair file"): off leaves the copy out of full archives, for a very large galaxy. The settings backup made before each save still holds the settings part alone.
- Restoring the part is DB.41; until then a restore of a whole archive leaves the parity copy in the archive and says so.

generation-output: unchanged
