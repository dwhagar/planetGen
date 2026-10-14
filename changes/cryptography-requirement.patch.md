### Fixed
- `update.sh` and `install.sh` on a venv no longer stop at the database step with "'cryptography' package is required for sha256_password or caching_sha2_password auth methods". `cryptography` (what pymysql needs to log in to MySQL 8's default accounts) is now a requirement and in the locks, and the library check fails a venv where pymysql cannot use it, so the library step installs it.
