### Fixed

- A parallel run no longer stops when the control database drops mid-run: the task rows are only a record, so the run finishes its sectors and its lease goes stale on its own (TEST.20).
- When a worker process dies, the tasks it never got are recorded as cancelled and the ones in flight as failed, instead of being left as running (TEST.20).
- Cancelling a parallel run from the Generate page (or any SIGTERM to its process group) now ends it as cancelled: each worker rolls back its unfinished sector and the pool stays whole, instead of the workers being killed and the run recorded as failed (TEST.21).
