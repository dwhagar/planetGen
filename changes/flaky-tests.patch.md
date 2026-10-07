### Fixed

- A Generate page job whose runner is slow to start (a busy server) is no longer shown as interrupted while the runner is still alive, which could let a new job start and then find the old one running again (TEST.90).
- Several intermittently failing tests no longer depend on random draws or timing (TEST.82, TEST.84, TEST.86, TEST.88, TEST.89, TEST.71).
