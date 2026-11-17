**Headline:** The Queue page says why each job has not started, by class, with how long the oldest has waited.

### Added
- **A "Why jobs wait" panel on the admin Queue page, and `GET /api/admin/work/waiting` (OPS.233).** For each job class it shows the jobs waiting on another job and the jobs ready that the gate has not let start, the age of the oldest of each, the next retry, the reasons the gate wrote for the ready ones (`engine_jobs.waited`) with how many jobs gave each, how long the jobs that started in the last hour waited from their making to their start, and the memory gate's figures in words (`memorygate.gate_words`). One read of the control database in the new `queue/waiting.py`; it writes nothing. What the controller last did stays in the supervisor panel's list of acts (OPS.214), which the memory controller of OPS.231 will fill. Times show as hours, minutes and seconds. The route is new and nothing is removed, so `API_VERSION` stays; the route snapshot gained the line.

generation-output: unchanged
generation-affects: none
