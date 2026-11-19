**Headline:** A test that runs a stage no longer waits for ever when the test box sets only the `PLANETGEN_TEST_MYSQL_*` account (TEST.801).

### Fixed

- The test queue's pump (`tests/inline_queue.run_ready_rows`) opened the control database with the process default account, built from `PLANETGEN_MYSQL_*`. On a box that sets only `PLANETGEN_TEST_MYSQL_*` the login was refused, the pump swallowed the error and ran no unit, and a stage's follower (`units.follow_stage`) waited for a stage end that never came (test_habitability_fill hung). The pump now logs in with the test server's account, and only an unknown control database ends it quietly.
