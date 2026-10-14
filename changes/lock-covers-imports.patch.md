### Changed
- A test now fails when the code imports a library that `requirements-server.lock` lacks, so a venv built from the lock alone cannot be missing one. A fresh venv from the lock was checked: every planetgen module imports, `migrate` and `check-math` run, and `/api/health` answers.
