**Headline:** The update's reload check no longer fails just because the site is closed for the update.

### Fixed
- **`update.sh` reported "does not answer" for a healthy site (OPS.75).** The reload check asked `/api/health` while the update's own maintenance mode was still on, and got a 503. `/api/health` now stays open while the site is closed and reports `maintenance` (cause and reason). The warning names the failing check, the address and the answer, and says whether the site is in maintenance mode and why. The failing step is named for the reload, not "step 9/9".
- **`update.sh` listed the `state/` folder as an uncommitted local change that would be overwritten (OPS.76).** That folder is the update's own (the maintenance switch and lock), and is no longer listed.
