**Headline:** the rule that decides which waiting jobs start now is written for the whole backlog at once, so a visitor's page job waits under a second however long a stage's backlog is.

### Added

- OPS.227: `queue/dispatch.py`, the dispatch rule the job engine's supervisor will ask each turn: the site's lanes first, oldest first, against the places kept for them; then every other lane ranked and judged by `scheduler.judge`; only the top 64 of each group judged; each waiting job's last reason returned and written onto its ledger row. A reservation made for a job that has waited 120 s lasts until that job runs (`scheduler.judge(reserve_until=)`). In the stage 0 simulator (300 units) a visitor's job waited over half an hour at the 99th percentile when the whole backlog is ranked as one, and under 1.5 s with this rule. Nothing calls it from a worker yet: OPS.228 and the other stage 4 items move the lanes onto it, and the old per-worker tickets go with them.
