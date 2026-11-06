**Headline:** Every line the update and install scripts print now carries a time, and both end with a summary that repeats their warnings.

### Changed

- `update.sh`, `install.sh` and the scripts they run (`scripts/deploy-common.sh`, `install-python-deps.sh`, `enable-swap.sh`, `examples/apache/setup-debug-log.sh`, `set-permissions.sh`, `create-cache-dir.sh`) print through `scripts/logrecord.sh`: about 200 plain `echo` lines became `log_info`, `log_notice`, `log_warning` and `log_error` records with a local time, its zone, a severity and a component (OPS.98). A heredoc of instructions or a tool's output goes through the new `log_lines`; `git fetch` and `git reset` run quiet. `install.sh`'s `== n/8 ==` headers are `Stage n of 8` dividers.
- A command line's output that is not a terminal (a pipe, a file, `tee`, an update run) is now display lines with the time, not bare text: `log.apply_debug_arguments` asks for it and `log.times_on_console` decides; a terminal keeps the plain text. `PLANETGEN_LOG_TIMES=always` or `never` overrides.
- `update.sh` and `install.sh` end with `log_summary`: how long it took, the code version, the schema, when the site was reloaded and answered, the warm-up result, then every warning raised, repeated with its headline (OPS.100). The reload's "Done." lines are "Site reloaded: ..." so "done" appears once, last. Long warnings (the health check, the warm-up) are a headline and an indented fix line.
- `docs/design/operations/logging-standard.md` has a section for the shell scripts and the summary; `test_log_shell.py` fails a bare `echo` in these scripts and tests `log_lines`, `log_duration` and `log_summary`.

generation-output: unchanged
