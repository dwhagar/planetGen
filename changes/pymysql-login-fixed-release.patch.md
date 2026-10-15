### Changed

The MySQL 8 login hook (`planetgen/db/login.py`) stays in place when PyMySQL releases its fix and is tested against both the current pymysql step and a fixed one on a fake MySQL 8 server without TLS, so the login works on either release (OPS.49 follow-up).
