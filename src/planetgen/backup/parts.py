# planetgen/backup/parts.py

"""
The backup parts this build provides (ADM.87): the configuration, `config.json`, and the DB repair parity data (DB.40). The
control database (ADM.88) registers here when it is built.
"""

import os

from planetgen.backup.core import BackupError, Part, register

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

    def restore(self, directory, context):
        return "Not restored: putting the repair file back is DB.41; the archive keeps it."


register(ParityPart())
