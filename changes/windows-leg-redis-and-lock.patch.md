### Changed
- TEST.115: the Windows CI job keeps Redis in WSL alive (it ran in the foreground of a wsl.exe it holds open), waits until Windows can connect and passes the working URL on to the tests. A job's lock is removed with a few retries, and three Windows-only test failures now say what state they saw.
