**Headline:** tests only: a fixture that hung once (the browser tests' server, say) is no longer waited for again by every later module.

### Changed

- TEST.190: a session fixture that fails is raised again for its dependents at once, but a module- or class-scoped one is set up again by each module, so a server that never came up cost every module its whole 900 s timeout. The new `tests/dead_fixtures.py` plugin (loaded in `pytest.ini`) remembers a fixture of a scope wider than one test whose setup `pytest-timeout` stopped, and fails the later tests that need it at once, naming the test that waited first. A function-scoped fixture is left alone. `PLANETGEN_TEST_RETRY_DEAD_FIXTURES=1` turns it off. `tests/test_dead_fixtures.py` shows three modules waiting once, not three times. No program code changes.
