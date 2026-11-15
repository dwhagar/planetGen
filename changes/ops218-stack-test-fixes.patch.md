**Headline:** The job engine stack's tests start no supervisor of their own and the route snapshot knows the supervisor route.

### Fixed
- The fixtures a module makes (a sample galaxy for the page tests) start their run through one turn of the supervisor's duty, as a test's own fixtures do, so they no longer wait for a supervisor that is not there.
- The API route snapshot lists `/api/admin/work/supervisor`.
- The debug-run test puts debug logging back after itself, so the tests that follow it in a process do not run slowly.
