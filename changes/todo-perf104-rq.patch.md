**Headline:** The cache warm-up item now makes every kind of warm-up a visible task in the admin queue.

### Changed

- Widened PERF.104 on Boss's follow-up (2026-10-15, 9:57 PM EDT): each kind of warm-up is its own RQ job shown on the admin runs page with its status and log, the command-line script enqueues those jobs by default, and an option such as `--no-rq` runs them inline instead.

TODO items resolved: none / Added: none
