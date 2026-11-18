**Headline:** The Queue page has no worker count or Ludicrous Speed any more: nothing reads it, because the job engine sizes every stage by memory.

### Removed
- **OPS.235 slice 1c.** The Queue page's Workers panel (the count form, Ludicrous Speed, the host worker cap line and its refusals), its `workers-set` and `workers-ludicrous` actions, `apiclient.admin_work_workers`, `GET` and `POST /api/admin/work/workers` (with the `WorkerCount` schema) and `workers.wanted`, `set_wanted`, `ludicrous_count` and `WANTED_KEY` with their tests. The launcher sizes every stage by memory (`workers.stage_worker_cap`) and the gate, as it did since OPS.231; the admin's count had not started a worker since slice 1b.
- Reuse: nothing new is built; the profile-led stage cap (`stage_worker_cap`) is the one count.
- API version 14: the route snapshot, `docs/reference/api-versions.json`, `api.md` and the compatibility page carry the removed routes.
- Still to go in slice 2 with `cli/worker.py`: the "workers active" tile, `host_worker_cap`, `host_worker_processes` and the `queue.max_host_workers` setting, which only RQ workers used.

api-version: changed
