# planetgen/backup/parts.py

"""
The backup parts this build provides (ADM.87): the configuration, `config.json`. The control database (ADM.88)
and the DB repair parity data (DB.40) register here in 9.0.
"""

import os

from planetgen.backup.core import BackupError, Part, register

CONFIG_FILE = "config.json"


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
