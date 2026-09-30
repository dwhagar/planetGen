### Changed
- **`install.sh` now works on an externally managed Python (PEP 668),
  such as Ubuntu 26.04 LTS's.** Where pip used to fail with
  "externally-managed-environment", the installer now detects the
  `EXTERNALLY-MANAGED` marker and installs the libraries as apt packages
  instead, pip-installing only what the distribution lacks (or packages too
  old) into a venv at `/opt/planetgen/venv`. On an unmanaged Python it still
  uses pip exactly as before. It prints which path it took; see
  `docs/apache-deployment.md`'s "Managed Python".
