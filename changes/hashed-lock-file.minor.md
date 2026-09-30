### Security
- **pip installs only locked, hash-checked files.** `requirements.lock`
  pins every runtime library and its dependencies to one version (per
  Python range) with the sha256 hashes of its files, and
  `scripts/install-python-deps.sh` installs from it with
  `--require-hashes` on both the ordinary and the externally managed
  path. apt-provided libraries are still used as they are.
  `scripts/lock-requirements.sh` (uv) regenerates it; a test fails when
  it no longer meets `setup.py`'s floors, and CI audits it with
  pip-audit.
