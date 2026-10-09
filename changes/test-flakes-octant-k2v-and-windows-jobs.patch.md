### Fixed
- A Windows job whose runner was slow to start (more than 15 seconds) was called interrupted, because the process-creation check used the short grace period; it now uses the start limit.
- Releasing the job lock on Windows failed silently (the lock file was removed while still open), so a job that failed to start kept the lock.
- Tests: the octant test finds a point inside a warped nebula even when no metaball centre is inside it; the end-to-end K2V test reads the primary star, not whichever star the database returns first; the Windows no-Redis job test no longer swaps in POSIX process options.
