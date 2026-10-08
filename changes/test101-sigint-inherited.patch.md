### Fixed
- The parallel-galaxy interrupt test starts its run with Ctrl+C at its default, since a pytest worker can leave it ignored and an ignored Ctrl+C is inherited, so the run finished with status 0 in the full suite (TEST.101).
