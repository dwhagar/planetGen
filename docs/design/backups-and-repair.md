# Backups and repair: one archive, and the ways to use it

**Status (2026-10-16):** plan; not built. The archive and the settings part are 8.1 (ADM.87, ADM.83 to ADM.85); the control database and parity parts, the restores and the repair chain are targeted for 9.0 (umbrella ADM.91).

Boss (2026-10-16, 1:59 PM to 2:00 PM EDT): "let's add backups of the control DB and backups of the parity data for DB repair", "All those should be bz2 tarballed together", "we also have to build a way to use that data so build that out and let's make this targetted for v9.0 ... Settings page is v8.1".

## 1. The archive

Every backup is one `backups/<UTC timestamp>.tar.bz2` taken at one moment (stdlib `tarfile`, mode 0600, web user). It holds `manifest.json` (format 1, time, who, version key, and per part its name, directory, sizes and SHA-256) and one directory per part: `settings/` (the user's `config.json`, secrets included), `control/` (a `mysqldump` of the control schema without sessions) and `parity/` (the parity file cut at its last complete record). The layout is fixed in ADM.87 for all three, so the later parts add providers and not a format. Retention is the option `settings.backups_kept` (default 20, 1 to 200) and counts the manual and settings-save archives, while nightly archives (OPS.55) are counted by `settings.nightly_backups_kept` (default 7, 1 to 60), so one never pushes out the other; the newest archive that holds every registered part is never pruned. A save on the settings page writes a small archive with the settings part alone; "Back up now" and the command write every registered part.

What it does not hold: the galaxy database. It is too large, can be regenerated from its seed (DB.17) and is repaired from parity (DB.9).

## 2. Using it

| Need | Item | What happens |
|---|---|---|
| Settings back | ADM.85 | restore the whole archive or the settings part, with a diff |
| Control database back | ADM.88, ADM.89 | dump into the archive; restore through a scratch schema, then replace, behind the update guard |
| Parity file back | DB.40, DB.41 | copy under the writer's lock; restore, then the DB.8 check recomputes stale records |
| Galaxy repaired | DB.42 | `repair-db --from-archive` uses the archive's parity part when the live file is missing or damaged; what parity cannot repair goes to the seed regeneration of DB.17 |
| One place for it | ADM.90 | the check, the damaged list, repair and part restores as queue jobs |
| Nightly | OPS.55 | a `backup` step in the daily maintenance run takes every part; schedule `settings.nightly_backup_hour`, own retention `settings.nightly_backups_kept` (default 7) |
| Proof | TEST.180, DOC.56 | end-to-end scenario; the deployment guide and help page |

## 3. Decisions and their reasons

- **One mechanism.** The settings page lists the archives and restores them whatever parts they hold, so no second backup screen exists.
- **Layout first.** ADM.87 reserves the directories and manifest entries now, and an unknown part is skipped, reported and kept when an archive is copied, so 9.0 adds providers without rebuilding the archive.
- **Restores are guarded.** Every restore is refused while maintenance, an update or a run is active (OPS.41, OPS.53, OPS.54) and goes through a confirmation that names what will be lost.
- **Secrets.** The control part holds password and key hashes: admins may download an archive now, only the Owner from 8.2 (ADM.86).
- **Parity can be switched off** for a very large galaxy (`settings.backup_parity`, default on); the manifest records its size and coverage.

## 4. Before 9.0

DB.9 (the parity file and `repair-db`), DB.8 (the check, done), DB.17 (seed regeneration), the update guard (OPS.41, OPS.53, OPS.54, done), and the settings page with ADM.87 (8.1). No migration is needed; the control schema is dumped and restored as it is.
