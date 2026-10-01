### Added
- **Test suite markers and more database engines in CI.** Every test is
  marked `db`, `slow` or `browser` as it applies, so `pytest -n auto -m
  "not db and not slow"` is a one-minute loop. CI's test job now runs on
  MySQL 8.0, MySQL 8.4 and MariaDB 11.4 (local runs cover MariaDB 10.11).
  `test_todo_tags.py` takes its category list from `bump_version.py`, so
  `TODO(TEST.N)` tags are accepted.
- **A guide to running CI on your own computers** (`docs/ci-runners.md`):
  what each runner needs, how to register it, security settings,
  troubleshooting. Pull requests from forks now always run on
  GitHub-hosted runners, never on self-hosted ones, and the browser job
  no longer needs passwordless sudo.
