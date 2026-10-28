**Headline:** The Generate page has a task that brings a stored galaxy to the running generator epoch, showing what it would redo before anything runs.

### Changed
- **"Bring the galaxy to the running generator epoch" on the Generate page (GEN.764, step 4).** The task page reads the galaxy (changing nothing) and says what a run would do: the stages that are behind, the first sectors it would generate again and how many more, the sectors left alone because an admin edited them, and the species to make again. Review shows the one command (`planetgen bring-to-epoch`) and Run starts it as a job like any other, so it appears in the job list with its log and a stopped run is resumed by running the task again. The task refuses to start, and says why, when the plan was made under an epoch that changed the galaxy's shape (a full regeneration is needed), when a scatter is behind (draw it again first), when the admin audit log cannot be read, or when the galaxy is already at the running epoch.
- `galaxy.bring_to_epoch.status` is the read-only look the page uses; `report_lines` takes a `limit` so the page names a few sectors of each kind and counts the rest.

TODO items resolved: GEN.764 (step 4). Added: none.
