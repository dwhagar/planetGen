**Headline:** One page test no longer fails on a machine whose MySQL has no `planetgen` login.

### Fixed
- **`test_web_pages.py::test_in_process_transport_skips_http_inside_a_request` (TEST.370)** asked the app for `/databases`, which lists the schemas on the configured server, with the default configuration's `planetgen` login and no password. On a test machine that has only the test login it failed with `DatabaseLoginError` (an `internal server error` from the API). The test now builds its app with the test server's login.
- Not changed: the runner's check of the test environment (TEST.376) already stops a lane that has no `PLANETGEN_TEST_MYSQL_USER` set, which is how the same gap shows up in other tests.
TODO items resolved: TEST.370, and TEST.390 (TODO closes them after reading this: each of Aux 10's tests passes on main; the runner self-tests, test_stages, test_web_fresh_after_cli and the browser map and table tests ran 116 passed, and the other four were fixed by #3189, #3280 and #3291). Added: none.
