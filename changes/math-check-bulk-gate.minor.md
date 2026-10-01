### Added
- **Bulk generation checks the math first (TEST.68).** `generate.py
  check-math` runs the math check by hand (`-v` lists every check). Every
  bulk run (`galaxy`, `plan`, `population`, and `sector` with more than one
  sector), every Generate page job (its new first step, "Check the math")
  and the Sector page's neighbourhood generation run it first and refuse
  to start if a check fails, naming the failed checks and writing nothing.
  `update.sh` and `update.ps1` run it after updating (a new step 4) and
  warn if it fails, skipping the population pass.
