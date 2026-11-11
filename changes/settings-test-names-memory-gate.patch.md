**Headline:** The settings test names `PLANETGEN_MEMORY_GATE` as a one-process switch, as it does the other environment variables that are not site options.

### Changed

- `tests/test_settings.py::test_the_code_reads_no_variable_the_model_does_not_name` was red on main: `queue/memorygate.py` reads `PLANETGEN_MEMORY_GATE` (PERF.254; `0` turns the gate off for the process so every claim is admitted at once), a variable that is neither a `config.json` option nor in `NON_OPTION_VARIABLES`. It is a switch for one process (a benchmark or debugging run), not a site setting, since the gate is the rule on every install, so it joins the list with that reason. The code is unchanged.

generation-affects: none
