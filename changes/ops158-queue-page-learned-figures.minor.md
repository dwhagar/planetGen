**Headline:** The Queue page now shows what the queue has learned for each kind of job, next to the typed figure, with the latest changes it made, and lets an admin fix a claim, release it or make the queue forget a kind.

### Changed
- **OPS.158.** A panel "What the queue has learned" on `/admin/queue` lists, per job kind, the typed claim, the claim worked out from the job profile, a figure an admin fixed, the samples behind it and the newest change, and below it the latest 20 changes the queue and the admins made (each with its reason). New API: `GET /api/admin/work/tuning` and `POST /api/admin/work/tuning` (`pin` with `mb`, `unpin`, `clear`), audited as `work.tuning.*`. A fixed claim wins over the profile until it is released (`PLANETGEN_SELF_TUNING=0` does not lift it); forgetting a kind deletes its settled claim, its largest-peak memory, its starting limit and its profile rows, so its next run claims from the typed table and measures again.
No drawn value changes, so there is no epoch or golden. The API gains routes only (no version bump).
generation-output: unchanged
