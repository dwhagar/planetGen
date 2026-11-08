**Headline:** `planetgen system --output` with an unwritable path now prints a one-line error instead of a traceback.

### Fixed
- **`planetgen system --output <bad path>` printed a Python traceback (UX.288).** The page write in `run_system` now catches the operating system's refusal, prints `error: the system page could not be written to <path>: <reason>` and exits 1.
- **The system-bodies API used `math.dist` and `math.hypot` (GEN.1389).** `api/system_bodies.py` now uses `detmath.dist` and `detmath.hypot` like the rest of the package, so the barycentre distance and the central offset in km do not depend on the Python release. The sites the bug hunt named in `db/sky_tier_ledger.py`, `generation/sky_charting.py` and `galaxy/stellar_table.py` were already on `detmath` (TEST.663); the one place `test_reproducible_draws` still refused was this API. The displayed values can differ from the old ones in the last binary digit (about 1 call in 5); nothing stored reads them.

generation-output: unchanged
TODO items resolved: GEN.1389, UX.288
