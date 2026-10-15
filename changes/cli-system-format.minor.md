### Changed

- `planetgen system --output` writes plain text by default and takes `--format text|markdown|wikitext` (it used to write wikitext, or Markdown with `--markdown`; that flag is gone from `system`). The admin one-off system page passes the format it was asked for. `sector`, `galaxy` and the phenomenon commands keep `--markdown` until they are described the same way.
