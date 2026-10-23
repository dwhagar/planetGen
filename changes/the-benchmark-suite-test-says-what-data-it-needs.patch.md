**Headline:** A benchmark suite test that needs the downloaded data now skips, saying so, on a machine without it.

### Fixed
- **`test_a_real_run_of_the_plan_writes_the_json_and_drops_its_databases` failed on a machine without the static downloads (TEST.285).** The suite starts the real plan step as a process, which refuses to run while a downloaded data file is missing (GEN.590, and Boss's missing-files rule), and the test then found exit status 1 where it expected 0. The suite times the real generator, so it takes no stand-in data; the test now looks for the missing files first and skips with their names, URLs and the setup command (`python3 -m planetgen.cli.datasets fetch`), and runs as before where the data is there.

TODO items resolved: TEST.285. Added: none.
