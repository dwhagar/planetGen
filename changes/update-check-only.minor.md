### Changed
- **planetGen's Python libraries live in their own venv,
  `/opt/planetgen/venv`, and the system Python is left to apt.** Nothing
  is pip-installed into `/usr/lib` or `/usr/local` any more, and
  `--break-system-packages` is never used. The venv is isolated and built
  from the system `python3` (the Python mod_wsgi embeds). Apache runs it
  through `python-home=/opt/planetgen/venv` on `WSGIDaemonProcess`
  (added to `examples/apache/planetgen.conf.example`).
  `/usr/local/bin/planetgen` and `planetgen-orbits@.service` run the
  venv's Python too. On an existing server, the next run removes the old
  `planetgen-venv.pth` from the system Python and rebuilds the old
  shared-packages venv as an isolated one. Until the vhost has the
  `python-home` line (install.sh and update.sh print it), `wsgi.py` adds
  the venv to its own path, so the site keeps working.
- **`update.sh` no longer reinstalls anything.** It used to re-run
  `install.sh` after every pull, which force-reinstalled the package and
  rebuilt the venv. Now it imports each library with the venv's Python and
  installs only one that is missing, below `setup.py`'s floor or broken
  (`scripts/install-python-deps.sh --check`). It rebuilds the venv only
  when a distribution upgrade changed the Python under it. It prints one
  line per library (present, installed, upgraded, repaired or failed), and
  finishes by importing the web app as Apache's user.
- The NLTK corpus and mod_wsgi are installed, and Apache's modules
  enabled, only when missing. These steps live in
  `scripts/deploy-common.sh`, shared by `install.sh` and `update.sh`. Both
  scripts warn when mod_wsgi is built for a different Python version than
  the venv, and both accept `PYTHON=` to build the venv from another
  interpreter.
- The admin Stats page shows the web app's Python prefix and the
  directory it imports its libraries from.
