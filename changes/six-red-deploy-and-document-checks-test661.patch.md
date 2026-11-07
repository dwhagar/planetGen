**Headline:** Six checks that were red on main pass again: a bash 4 expansion in the deploy script that macOS's bash 3.2 cannot run, three test harnesses that missed the log library, seven item codes in text users see, and a variable the settings check did not know.

### Changed

- `scripts/deploy-common.sh` (TEST.661): the cache warm-up's skip key was compared with `${WARM_SKIP_KEY^^}`, a bash 4 expansion; the script says it runs under macOS's bash 3.2 (no `mapfile`, no associative arrays), so a Mac would have stopped at the line. The upper case spelling is a second constant, `WARM_SKIP_KEY_UPPER`. This was a real portability bug; the `^^}` check of `test_install_python_deps.py` was right.
- `src/tests/test_install_python_deps.py`: the harnesses that run one function of `install-python-deps.sh` (`install_checkout`, `install_venv`) source `scripts/logrecord.sh` first. OPS.98 made the script print through `log_info`, which the extracted function calls and the harness did not define, so three tests ended with `log_info: command not found` (127). The code was right; the tests were behind it.
- `src/tests/test_no_todo_codes_shown.py` found item codes in text a user sees: the HYG dataset's `replaces` and `fill` (shown as `hyg_stars (VIEW.30) is missing`), the sheet audit's excuse for the cloud composition columns, the charting log line, a debug line of the life bar, the circular speed check's note and the `sky_limit` setting (and so `docs/reference/config.md` and `config.schema.json`, regenerated). Each is now plain words; `test_datasets_readers.py` names the new `fill`.
- `src/tests/test_settings.py`: `PLANETGEN_SYSTEMD_UNIT_DIR` (SEC.43's Update button, read in `admin/webupdate.py`, set by its tests) joins the variables that are not options, with its reason. It was red on main with the other six.
- The plain-words document test passed on main when run alone (69 passed with `test_no_todo_codes_shown.py`), so nothing is changed for it; the failure the Train Manager saw is not reproduced.

generation-output: unchanged
