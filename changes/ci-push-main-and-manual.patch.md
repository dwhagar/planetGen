### Changed
- GitHub CI (`ci.yml`, every test leg) now runs on its own only for a push to `main`, and by hand on any branch (Actions > CI > Run workflow). It still does not run on a pull request or on a push to another branch. A newer run on the same branch cancels the older one that is still going, so a merge's run can be cancelled by the "Release x.y.z" commit pushed right after it.
