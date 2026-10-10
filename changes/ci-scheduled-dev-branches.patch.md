### Changed
- GitHub CI (`ci.yml`, every test leg) no longer runs on a push to `main` or on any push or pull request. It runs on a daily schedule (07:00 UTC) on `main`, where the same run also starts CI on every `dev-*` branch, and by hand (Actions > CI > Run workflow). A branch named `claude/*` or `claude-*` never runs the test jobs.
