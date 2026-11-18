**Headline:** The resume-across-builds runner tests no longer depend on the machine's MySQL and Redis, so one of them stops failing or hanging when those are missing or down.

### Fixed
- **`test_main_says_how_many_tests_were_reused_and_from_which_run`** (`src/tests/test_resume_across_builds.py`): it ran `scripts/test.py quick --resume` with the real service check, which exits 2 when the test environment is not set and waits on a MySQL or Redis that is down when it is (a hang on main, seen by Lane 16 at 1:23 AM EST on 2026-11-18). The runner these tests build now has a quiet services check, as `test_test_runner.py` already does; the tests are about reusing passes, not about services.
