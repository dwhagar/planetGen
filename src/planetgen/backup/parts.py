# planetgen/backup/parts.py

"""
The backup parts this build provides (ADM.87): the configuration, `config.json`, and the DB repair parity data (DB.40). The
control database (ADM.88) registers here when it is built.
"""

import os
import shutil
from dataclasses import dataclass, field

from planetgen.backup.core import BackupError, Part, fsync_folder, match_owner, register

CONFIG_FILE = "config.json"
PARITY_COPY = "sector-parity.pgpr"
COPY_CHUNK = 8 * 1024 * 1024


class SettingsPart(Part):
    """`settings/config.json`: the user's own configuration file, secrets included (so the archive is 0600 and its
    download is for the Owner). With no config.json the part is empty."""

    name = "settings"
    label = "Settings"

    def dump(self, directory):
        from planetgen.util import settings_writer
        data = settings_writer.read_bytes(settings_writer.config_path())
        if data is not None:
            with open(os.path.join(directory, CONFIG_FILE), "wb") as out:
                out.write(data)

    def restore(self, directory, context):
        from planetgen.util import settings_writer
        path = os.path.join(directory, CONFIG_FILE)
        if not os.path.isfile(path):
            raise BackupError("This archive's settings part is empty: there was no config.json when it was taken.")
        with open(path, "rb") as handle:
            data = handle.read()
        return settings_writer.restore_bytes(data, context)


register(SettingsPart())


@dataclass
class ParityRestore:
    """What restoring the parity part did: the file put back (`bytes`, the `sectors` and `groups` it holds) and `found`, the `db.parity.Reconciled`
    of its check against the database."""
    bytes: int
    sectors: int
    groups: int
    found: object = field(default=None)

    def lines(self):
        found = self.found
        out = [f"Put back the repair file: {self.bytes:,} bytes, {self.sectors:,} sector(s) in {self.groups:,} group(s).",
               f"Checked against the database: {found.checked:,} sector(s), {found.healthy:,} matching, {found.recomputed:,} records made again "
               f"from the rows (the database had changed since the backup; not damage), {found.released:,} let go (sector gone)."]
        if found.unprotected:
            out.append(f"{found.unprotected:,} sector(s) are not protected until their next save (or `planetgen refresh-parity`): the backup holds no record of them.")
        if found.damaged:
            shown = ", ".join(f"{sector_id} ({state})" for sector_id, state in found.damaged[:10])
            out.append(f"{len(found.damaged):,} sector(s) are damaged or cannot be trusted: {shown}{' ...' if len(found.damaged) > 10 else ''}. "
                       "Run `planetgen repair-db`.")
        return out

    def __str__(self):
        return "\n".join(self.lines())


class ParityPart(Part):
    """
    `parity/`: a copy of the sector parity file (`db/parity_file.py`, DB.9), the one thing a damaged galaxy database can be repaired with. A repair
    works only with a parity file that matches the database, so the file is a backup part (DB.40).

    The copy is taken with the file locked against every writer (`ParityFile.locked`, the lock the saves take), and cut at the end of the last
    whole record, so it is a state the file really had and never half a record, whatever is being saved. Nothing is indexed or stored beside the
    records: the file is an append-only log whose records describe themselves, and a header scan finds them again. The dump returns what the manifest
    keeps as the part's `info`: `bytes` copied, `sectors` the copy covers (the sectors the latest record of each group lists), `groups`, and
    `present` (false, with `reason`, when there was nothing to copy: no file yet, or the option `settings.backup_parity` is off).
    """

    name = "parity"
    label = "Sector repair file"

    def dump(self, directory):
        from planetgen.db import parity, parity_file
        from planetgen.util.settings import get_settings
        if not get_settings().settings.backup_parity:
            return {"present": False, "reason": "settings.backup_parity is off"}
        path = parity.parity_path()
        if not os.path.isfile(path):
            return {"present": False, "reason": "there is no parity file yet"}
        source = parity_file.open_file(path)
        target = os.path.join(directory, PARITY_COPY)
        with source.locked():
            if source.damaged_at is not None:
                raise BackupError(f"The parity file {path} is damaged at byte {source.damaged_at}: it would make a backup that cannot repair "
                                  f"anything. Run `planetgen repair-db` (it makes the file again from the database) and back up again.")
            copied, wanted = 0, source.valid_end
            with open(path, "rb") as reading, open(target, "wb") as out:
                while copied < wanted:
                    block = reading.read(min(COPY_CHUNK, wanted - copied))
                    if not block:
                        break
                    out.write(block)
                    copied += len(block)
                out.flush()
                os.fsync(out.fileno())
            sectors = sum(1 for entry in source.entries.values() for member in entry.members if member.sector_id != parity_file.EMPTY)
            groups = len(source.entries)
        if copied != wanted:
            raise BackupError(f"The parity file {path} is shorter than its index says ({copied} of {wanted} bytes): it changed under the copy.")
        return {"present": True, "bytes": copied, "sectors": sectors, "groups": groups}

    def check_restorable(self, context):
        """Refuses at once, before the safety archive is taken, while a save is writing parity (the safety copy would otherwise wait for the save and
        the restore go on behind it)."""
        from planetgen.db import parity, parity_file
        path = parity.parity_path()
        if not os.path.isfile(path):
            return
        try:
            with parity_file.open_file(path).locked(wait=False):
                pass
        except parity_file.ParityFileBusy as exc:
            raise BackupError("A save is writing the repair file right now, so it was not restored: try again when the save is done.") from exc

    def restore(self, directory, context):
        """
        Puts the backed-up repair file back (DB.41). The copy is checked first (a whole file of whole records, nothing damaged), the database
        is reached before anything is changed, and a save that is writing parity right now refuses the restore. The file is written beside the
        live one, synced and renamed over it with the live file locked, so a reader sees the old file or the new one. Then `db.parity.reconcile`
        checks it against the database (DB.8): a record that no longer matches the rows is made again (not read as damage), and the result says
        how many sectors are unprotected until their next save and which are damaged.

        Returns:
            ParityRestore: What was put back and what the check found.

        Raises:
            BackupError: The part is empty, its file is damaged, the database cannot be reached, or a save is writing parity; nothing was changed.
        """
        from planetgen.db import parity, parity_file, store
        copy = os.path.join(directory, PARITY_COPY)
        if not os.path.isfile(copy):
            raise BackupError("This archive's parity part is empty: there was no repair file to copy when it was taken (the manifest's info says why).")
        looked = parity_file.ParityFile(copy)
        with looked.locked():
            if looked.damaged_at is not None or looked.valid_end != os.path.getsize(copy):
                raise BackupError("The repair file in this archive is damaged (it is not a whole number of whole records); it was not restored.")
            groups = len(looked.entries)
            sectors = sum(1 for entry in looked.entries.values() for member in entry.members if member.sector_id != parity_file.EMPTY)
        try:
            conn = store.get_connection(context.database or store.MySQLConfig(), ensure_schema=False)
        except Exception as exc:  # noqa: BLE001 -- say what it was, change nothing
            raise BackupError(f"The galaxy database could not be reached ({exc}), so the restored file could not be checked against it: "
                              "nothing was changed.") from exc
        try:
            path = parity.parity_path()
            live = parity_file.open_file(path)
            temp = path + ".restore"
            try:
                with live.locked(wait=False):
                    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
                    shutil.copyfile(copy, temp)
                    with open(temp, "rb") as written:
                        os.fsync(written.fileno())
                    os.chmod(temp, os.stat(path).st_mode & 0o777 if os.path.exists(path) else 0o640)
                    match_owner(temp, os.path.dirname(os.path.abspath(path)))
                    os.replace(temp, path)
                    live._unlock_after_replace()
                fsync_folder(os.path.dirname(os.path.abspath(path)))
            except parity_file.ParityFileBusy as exc:
                raise BackupError("A save is writing the repair file right now, so it was not restored: try again when the save is done.") from exc
            finally:
                if os.path.exists(temp):
                    os.unlink(temp)
            found = parity.reconcile(conn, parity.store())
        finally:
            conn.close()
        return ParityRestore(os.path.getsize(path), sectors, groups, found)


register(ParityPart())
