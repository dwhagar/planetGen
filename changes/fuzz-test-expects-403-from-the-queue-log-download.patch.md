**Headline:** The route fuzz test no longer expects a redirect from the queue run's log download, which has answered 403 to a signed-out visitor since SEC.36.

### Fixed
- **`test_fuzz_web_routes.py::test_admin_pages_bounce_anonymous_visitors` failed on `/admin/queue/<run>/log` (TEST.396).** SEC.36 made that route answer 403 to a signed-out visitor, as the other admin JSON routes do, because it is a download and not a page (a redirect to the login page would be saved as the log file); `test_api_auth_sweep` holds that, and the fuzz test's list of such routes (`ADMIN_JSON_ENDPOINTS`) was not updated with it. I ran the sweep and not the fuzz test with SEC.36. The route is added to the list.
TODO items resolved: TEST.396 (TODO closes it after reading this). Added: none.
