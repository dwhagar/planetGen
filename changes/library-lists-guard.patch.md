**Headline:** A test now fails when install.sh or update.sh could miss a library that setup.py requires.

### Changed
- **`test_install_python_deps.py` checks the one list both scripts use.** install.sh and `update.sh --check` both run `scripts/install-python-deps.sh`; the test fails if either stops doing so, if the list differs from setup.py's requirements, if `psutil` is missing, or if a listed library is not pinned in `requirements.lock` and `requirements-server.lock`.
