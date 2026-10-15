### Fixed
- **The fail2ban jail and its guide named the wrong log file** (DOC.41). Both said `/var/log/planetgen/planetgen-log`, which is the logrotate file's name; the activity log is `planetgen.log`, so the shipped jail watched a file that does not exist and the guide's `tail` and `fail2ban-regex` commands failed. A test now checks the jail's `logpath` against the name the program writes.
