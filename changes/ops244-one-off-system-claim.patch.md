**Headline:** The one-off system page's light job now claims the memory of one drawn system, so the memory gate no longer runs it unaccounted.

### Fixed
- **A claim for `system_page.generate_in_process` (OPS.244, slice A).** The function the one-off system page runs in the warm light-job process had no entry in `claims.API_CLAIMS`, so `test_every_function_the_api_dispatcher_may_call_has_a_claim` failed on main. It claims what a `generate-system` unit does (`workers.BULK_WORKER_MEMORY_MB`, one core), the same work done in a light job.

generation-output: unchanged
generation-affects: none
