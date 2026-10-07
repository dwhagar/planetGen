### Fixed

- A Generate page job whose runner is slow to start (a busy server) is no longer shown as interrupted while the runner is still alive, which could let a new job start and then find the old one running again (TEST.90).
- A test now checks that bright stars spread far above the galactic plane, by population (GEN.117: the thin band seen on a server was the code from before 7.194.608).
- Several intermittently failing tests no longer depend on random draws or timing (TEST.82, TEST.84, TEST.86, TEST.88, TEST.89, TEST.71).
