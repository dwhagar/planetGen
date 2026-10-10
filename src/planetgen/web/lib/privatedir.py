# planetgen/web/lib/privatedir.py

"""
A private directory in a shared place (the system temp directory), for the
tile cache's and the Generate jobs' fallback when their configured
directory can't be created.

`/tmp` is shared with every local user, so a directory found there under
a well-known name may have been planted by someone else -- as a symlink
to somewhere else, or as a directory they can read or write. Such a
directory is refused, never reused.
"""

import os
import stat


def ensure_private_dir(path):
    """
    Creates `path` with mode 0o700, or checks that an existing one is
    safe to use: a real directory (not a symlink), owned by this
    process's effective user, and writable by nobody else.

    Returns:
        str: `path`, when it is safe to use.

    Raises:
        OSError: When it can't be created or isn't safe (refused, never
            fixed up: whoever planted it may still hold it open).
    """
    try:
        os.mkdir(path, 0o700)
    except FileExistsError:
        pass
    info = os.lstat(path)
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise OSError(f"{path} is not a directory (a symlink or a file); refusing it")
    if info.st_uid != os.geteuid():
        raise OSError(f"{path} is owned by uid {info.st_uid}, not this process (uid {os.geteuid()}); refusing it")
    if info.st_mode & (stat.S_IWGRP | stat.S_IWOTH):
        raise OSError(f"{path} is writable by other users (mode {stat.S_IMODE(info.st_mode):o}); refusing it")
    return path
