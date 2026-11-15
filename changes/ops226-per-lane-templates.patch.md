**Headline:** One warm template per lane, a path for hosts with no kernel handle on a job, and a memory proof on the real star tables; built and tested, not yet used.

### Added
- `jobhost.TemplateSet` starts one template per warm set when the first job that wants it is asked for: the lane whose units draw stars (`units.TABLE_KINDS`) gets the sampling tables built once, every other lane gets none (`WARM_SETS`, `template_key`). A host with no star-drawing work never builds the tables.
- A slow test (`tests/test_template_memory.py`) forks three star-drawing jobs from one template that built the real tables and checks that the template and its jobs hold well under what as many workers would each hold.

### Changed
- Where the kernel gives no handle on a job (macOS; `PLANETGEN_TEMPLATE_PIDFD=0` tries it on Linux), the template's report of each exit is the way, and a template that has died leaves the supervisor to look for each job's pid (OPS.226, the third slice).
