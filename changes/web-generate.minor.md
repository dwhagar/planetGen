### Added
- **Generate, plan and reset the galaxy from the web interface.** A new
  admin-only page, `/admin/generate` (the Generate link in the header),
  runs `generate.py plan`, `generate.py galaxy` (every mode: random start,
  whole ring, around a sector, one address) and `resetDb.py` as
  background jobs, plus a one-click "New galaxy" that resets, plans and
  generates a first neighborhood. Reset and New galaxy ask for the
  database name to be typed back. The running job shows its step, a
  progress bar, elapsed time and live output, and can be cancelled; the
  last 20 jobs keep their full output. New `jobs` section in
  `config.json` (`docs/config.md`).
- `generate.py` writes its progress to `$PLANETGEN_PROGRESS_FILE` when
  that is set (`stellarObjects/progressFile.py`).
