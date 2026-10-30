**Headline:** INSTALL.md now says to raise MySQL's redo log capacity on a machine that runs the full test gate.

### Changed

- DOC.73: INSTALL.md names `innodb_redo_log_capacity`, the value (2 GB instead of the 100 MB default, in `mysqld.cnf` or with `SET GLOBAL` on MySQL 8.0.30 and later) and why (under the parallel gate the log checkpointer lags and tests can lose their connection), in the requirements and in the development setup. Documentation only.
