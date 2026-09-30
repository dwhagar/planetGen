### Changed
- **`update.sh` no longer reinstalls anything.** It used to re-run
  `install.sh` after every pull that brought new commits, which
  force-reinstalled the Python package with pip, rebuilt the fallback
  venv and re-fetched the NLTK corpus. Now it checks instead: each Python
  library is imported with the system Python and compared with
  `setup.py`'s floor (`scripts/install-python-deps.sh --check`), and only
  one that is missing, too old or broken is installed, the same way
  `install.sh` would on that host (apt then the venv on an externally
  managed Python, pip on an ordinary one). The NLTK corpus and mod_wsgi
  are installed and Apache's modules enabled only when missing. It
  prints one line per library (present, installed, upgraded, repaired or
  failed), and finishes by importing the web app as Apache's user, so a
  library www-data can't use fails the update instead of the site.
- The NLTK and Apache module steps now live in `scripts/deploy-common.sh`,
  shared by `install.sh` and `update.sh`, and `install.sh` installs
  mod_wsgi when apt can instead of only warning about it.
- `/usr/local/bin/planetgen` is now the checkout wrapper on every host
  (pip's console script ran the pip-installed copy, which would go stale
  once updates stopped reinstalling it).
