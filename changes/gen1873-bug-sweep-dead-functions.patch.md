**Headline:** A bare local name's search job now has a memory claim, and six functions nothing called are gone.

### Fixed

- `planetgen.api.lookup:name_search_job` (API.38) had no line in `claims.API_CLAIMS`, so `test_every_function_the_api_dispatcher_may_call_has_a_claim` failed and the job ran on the default 1000 MB figure. It now claims 700 MB and one core, the figure of the other read-only light job.

### Removed

- GEN.1873: `memorygate.kill_note`, `units.open_job_id`, `wrapper.run_kind`, `region_stats.generated_mass_kg`, `corridor.unknown_space_flags` and `run_galaxy._neighborhood_batch`, which had no caller, test or wiring-register line.

generation-affects: none
generation-output: unchanged
