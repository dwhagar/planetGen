**Headline:** The Queue page lists the site's waiting light jobs, by lane, with why they wait, how long and when they are dropped.

### Added
- **The "Why jobs wait" panel lists the waiting light jobs (UX.320).** An API job, a slow page's answer and a visitor's prefetch run as light jobs since OPS.230: a record in Redis, no row in the ledger, so the panel of OPS.233 did not see them. `queue/waiting.py` gains `light_overview` (the pending list and each record, one pass, nothing written, the 500 oldest) and `GET /api/admin/work/waiting` a `light` key: for each lane (`interactive`, `long`) how many wait, how long the oldest has, the least time any has left before it is dropped for want of anyone asking (the deadline), the reasons the dispatch rule gave and what the jobs are. A job the dispatch rule has not judged yet says so. Times show as hours, minutes and seconds. Redis not answering leaves the ledger's half of the panel as it was.

generation-output: unchanged
generation-affects: none
