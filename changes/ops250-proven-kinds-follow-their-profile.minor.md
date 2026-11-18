**Headline:** A job kind whose profile holds a full window of finished units from this machine is now sized by what it measured, and the typed figure is only its default and its clamp.

### Changed
- **OPS.250 (absorbs OPS.160).** `profileclaim` marks a kind *proven* when its own history (this machine shape, this generator epoch) holds `PROVEN_SAMPLES` (30) finished units. A proven kind's memory claim and expected unit seconds (`profiletime`) are no longer held to 25 percent a run from the last figure or the typed one: the learned figure is the figure, held between the typed figure over `CLAMP` (4) and times it where there is one. A memory kill still wins above the clamp. A kind with fewer samples, from another machine shape or from another epoch behaves as before.
- The stage cut's worker count (`splitter.resources_now`) and the "largest worker held" line of a stage (`wrapper._log_peak_memory`) read the learned figure, not `workers.STAGE_WORKER_MEMORY_MB` alone; the line no longer tells an admin to raise the typed constant when the profile leads.
- The Queue page's learned-per-kind table has a Leads column (Fixed, Profile (proven), Profile, Typed).
- Reuse: `profileclaim.settle`, `claim` and `target`, `profiletime.target` and `unit_seconds`, `workers.stage_worker_mb` and `stage_worker_cap`, and the existing overview and change log; the new part is the proven flag and the clamp.
- Output unchanged, no epoch.
