**Headline:** The Queue page's "What the queue has learned" table lists only the job kinds the code still emits, and folds the retired ones away where they can still be forgotten.

### Fixed

- The table listed every job kind the job profile had ever seen, including `cut:` and `stage:` kinds from a run weeks ago whose command or stage no longer exists. A kind is now current while its command is still registered (a `unit:` or `cut:` kind), while it has a typed claim (an `api:` or `fn:` kind), while an admin has fixed its figure, or, for a `stage:` kind or any other name, while it had a sample in the last fortnight. Every other kind is retired.
- The retired kinds sit under a folded "Retired kinds" line below the table, with the same Forget control, so none is lost from view and the admin can still delete its samples.
- `GET /api/admin/work/tuning` gives each kind a `current` flag.
