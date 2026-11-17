**Headline:** The Galaxy Map's background node build (the nodes around a visit, and the nodes a change touched) is now handed out by the supervisor as light jobs, and a chunk nobody started within 15 minutes is dropped with its claims.

### Changed
- **The node prefetch and rebuild are light jobs on the maintenance lane (OPS.230 slice B).** `tileprefetch.schedule` and `rebuild` record each chunk of 16 nodes as a light job (`lightjobs.add_job`, the API function `tileprefetch:build_nodes`) instead of queuing it on the `planetgen-maintenance-nodes` RQ queue, and ask for the supervisor. The scheduler still weighs them as the site's (`scheduler.SITE_KINDS`, now keyed `api:planetgen.web.lib.tileprefetch:build_nodes`), so a visitor's nodes go before a generation run's units and the upkeep. The limit a chunk may run is what it was (ten times its estimate, at least ten minutes), now ended by the supervisor.
- **A chunk nobody started is dropped, not run late.** A chunk stays wanted for as long as its nodes stay claimed (`tileprefetch.CLAIM_SECONDS`, 15 minutes); one that has not started by then is dropped and its node claims are deleted, so the next visit asks again. `prefetch_max_queued` now counts the pending chunks of the prefetch itself (`lightjobs.pending_count`).
- Deleted: `redisqueue.NODES_QUEUE` (the maintenance lane has the two queues `planetgen-maintenance` and `planetgen-maintenance-check`), the prefetch's `workers.ensure_workers` call and its entries in `allowlist.FUNCTIONS` and `claims.FUNCTION_CLAIMS` (the function is on `allowlist.API_FUNCTIONS` and `claims.API_CLAIMS`).

generation-affects: none

generation-output: unchanged
