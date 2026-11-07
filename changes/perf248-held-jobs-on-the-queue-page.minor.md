**Headline:** The Queue page lists the jobs that are waiting for a job they require.

### Changed
- A new Held jobs panel on the admin Queue page (shown only while a job is held) lists each held job with its queue, its run and stage when it is a unit of a run, and the jobs it still waits for. It reads a new `GET /api/admin/work/held`.

generation-output: unchanged
