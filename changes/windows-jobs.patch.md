### Fixed
- **The admin Generate page's jobs work on native Windows.** Cancel now
  writes a `cancel` file that the job runner checks while a step runs,
  and stops the step's whole process tree (`os.killpg` on POSIX,
  `taskkill /T /F` on Windows), so a cancelled generation step no longer
  keeps writing to the database. Liveness uses `OpenProcess` and
  `GetExitCodeProcess` on Windows (with a creation-time check against
  reused pids), the runner starts detached in its own process group
  (breaking away from IIS's job object where allowed), `state.json`
  writes retry while the page has the file open, and the private
  fallback directory check no longer needs `os.geteuid`. CI runs the job
  tests on a Windows runner.
