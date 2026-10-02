### Changed

- The `+name`/`-name` forcing options (`+habitable_world`, `-planets`,
  `+comets` and the rest) now work only for a single system: `generate.py
  system` and the one-off system page (GEN.51). `generate.py sector` and
  `galaxy` refuse them with an error that names the option and points to
  `system`, so a saved or queued command line that still has one gets a
  clear message. `--min-habitable` still works for sectors. Prevalence
  controls for sector and galaxy runs come later (GEN.52).
