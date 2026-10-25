**Headline:** The run log download answers a plain 403 to someone who is not signed in as an admin, like the other file downloads, instead of a login redirect.

### Fixed
- **`/admin/queue/<run_id>/log` redirected to the login page for an admin still on the seeded credentials.** The view is marked `json_only` (a download, not a page), and the access sweep in `test_api_auth_sweep.py` expects those to answer 403; it was red on main. The view now aborts with 403 where `/admin/settings/export` already does.
TODO items resolved: none (the test was red on main with no item of its own). Added: none.
