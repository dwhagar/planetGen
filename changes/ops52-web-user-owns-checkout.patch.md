**Headline:** The whole planetGen install belongs to the web user (www-data:www-data), set the same way by install.sh and update.sh, and the example Apache2 site refuses to serve anything but the app and static files.

### Changed
- `examples/apache/set-permissions.sh` now takes the checkout (and optionally a state folder) instead of `src/html` and `db`: it gives everything in it to Apache's user and group, directories 750, files 640, scripts and `*.py` 750 (files that already have the owner execute bit keep it), makes `/var/lib/planetGen/state` if it is missing, and never follows a symbolic link. `install.sh` and `update.sh` both call it on the whole checkout, so a fresh install and an update end with the same ownership. The virtual environment and the NLTK corpus stay root's, because the root-run scripts execute that Python.
- `update.sh` trusts its own checkout for git (it is no longer root's), and says when an update changed `examples/apache/planetgen.conf.example` or `examples/redis/planetgen.conf`, with the `diff` command to compare. It still never edits your Apache site or your Redis config.
- `examples/apache/planetgen.conf.example` denies the whole install tree (no listings, no handlers, `AllowOverride None`), lets back in only `wsgi.py` and `static/`, and refuses source, configuration, state, logs, `.git`, dotfiles and the virtual environment by name even if a later `Alias` points at them. Update your own site file by comparing it with the example.

### Fixed
- A web process can write its own state folder on a fresh install as well as after an update (the maintenance flag could not be written when the install tree belonged to root).

TODO items resolved: OPS.52 / Added: none
