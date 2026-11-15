**Headline:** A warm template process and the jobs forked from it exist, built and tested, not yet used.

### Added
- `queue/template.py` is the template process: it builds a warm set, holds no connection and one thread, forks one job per request, and reports each job's exit code or killing signal at once. A job has a pipe to the supervisor and finishes its piece and exits (`EXIT_SUPERVISOR_GONE`) when the pipe closes.
- `queue/jobhost.py` is the supervisor's side: it starts the template, forks jobs, holds each job's pipe and a kernel handle (`pidfd`, opened by the template while the job is still its child) and reports every job's end once, even when the template died first.
- Nothing runs on it yet; the bulk lane moves to it in OPS.228 (OPS.226, the first slice).
