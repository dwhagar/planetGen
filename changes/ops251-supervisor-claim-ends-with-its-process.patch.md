**Headline:** The claim to be starting the supervisor ends when that process ends, not after 20 seconds.

### Changed
- `supervisor.ensure_started` writes its claim with no time limit. `redisqueue.start_detached` gives the exit notice of the process it started (`on_exit`), which ends the claim; a claim whose caller has gone is taken at once; a spawn that fails leaves none (OPS.251). The worker counter's 90 seconds goes with RQ's workers (OPS.271) and is left.

generation-output: unchanged
