### Added
- **A debug log for the whole program.** `"debug": true` in `config.json`
  (off when missing) makes the generator CLI, the maintenance scripts, the
  web pages and the API all write to `/var/log/planetgen.log` (`log_file`
  to move it): every decision the generator makes and why, every random
  roll with the source line that asked for it and the probabilities it was
  compared against, every SQL statement, web request, API call and admin
  access check, and every error with its traceback, all timestamped to the
  millisecond. A seeded run generates the same result with it on or off.
- **Log rotation for it.** `install.sh`/`update.sh` create the log file
  when debug is on (writable by Apache and shell users alike) and install
  `/etc/logrotate.d/planetgen`: daily, or as soon as it passes 100 MB
  (checked hourly), keeping 7 compressed copies.

### Changed
- **Web pages no longer show tracebacks when debug is on.** They went to
  the page itself before; now they go to the debug log, and the 500 page
  says so.
