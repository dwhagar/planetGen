**Headline:** The command line drops a workaround for a Python 3.9 argument-parsing bug that no supported Python has.

### Changed
- **`reject_lost_option_values` is gone from `planetgen/cli/generate.py` (DOC.24).** Python 3.9 turned `--workers=--` into an empty list; Python 3.12 and later keep the text and the option's own type check rejects it ("invalid int value"), which `tests/test_fuzz_cli.py` still covers for every such case.
