**Headline:** The self-tuning log prints worker counts as whole numbers again.

### Fixed
- **Self-tuning log line for workers and start limits.** `profileclaim.log_change` wrote "start limit of unit:link-neighbors was 8.0, now 3.0"; worker counts and start limits now print with no decimal ("was 8, now 3"), which `test_stage_worker_profile` checks. Memory claims keep one decimal.
generation-output: unchanged
generation-affects: none
