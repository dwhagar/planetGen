### Fixed

- When the database disk has no room for a run, the Generate page and a sector's "Generate more sectors around this one" now offer "Generate anyway" instead of only refusing (ADM.33). Choosing it starts the job and records the override in the activity log (`job.generate_anyway`).
