### Added
- One settings model describes every `config.json` option (ADM.42): type, default, help text, unit, secret, restart and editable flags, and the environment variable that overrides it. `python -m planetgen.cli.config check` validates the files (a misspelt option name is an error there) and `config docs` writes the option table in `docs/config.md`, `config.json.example` and `config.schema.json`. A web-owned `settings.json` overlay is read when present, holding only options marked editable from the web. Tests fail when the generated files drift or when code reads a `PLANETGEN_*` variable the model doesn't name.

### Changed
- A wrong value in `config.json` (a port out of range, a negative proxy count, an unknown `log_rotation`) now stops the program with a message naming the field, instead of being read as given. `PLANETGEN_ADMIN_COOKIE_INSECURE` reads text like `PLANETGEN_DEBUG` does (`false`, `0`, `no`, `off` and empty mean off) instead of only `1`.

### Removed
- `planetgen.util.appconfig`; the log locations moved to `planetgen.util.logpaths` and every other option to `planetgen.util.settings`.
