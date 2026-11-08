**Headline:** A test now stops bash 4 syntax from reaching the install scripts macOS runs.

### Fixed
- **No bash 4 in the scripts macOS runs (OPS.120).** `deploy-common.sh` used `${var^^}`, which macOS's bash 3.2 cannot parse, and three install tests died on `log_info: command not found`. The line is already gone from main; a new test scans `deploy-common.sh`, `install-python-deps.sh`, `update.sh` and `install.sh` for case-change expansions, associative arrays, `mapfile`, `readarray`, `|&` and `coproc`, so it cannot return. `enable-swap.sh` keeps `${SIZE^^}`: it runs only on Linux, where it needs `numfmt` and `df --output`.

generation-output: unchanged

TODO items resolved: OPS.120.
