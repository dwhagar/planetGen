**Headline:** Six quiet error handlers in the job profile code now leave a debug line.

### Fixed
- **OPS.168 on the self-tuning jobs (OPS.155 to 158).** `jobprofiles` (buffer pool read, connection close), `memorygate` (usage note at the heartbeat), `profileclaim` (pin read) and `usage` (processor time, resident memory) went on past a failure with no trace, which `test_every_broad_except_reraises_or_logs` flags. Each now calls `log.swallowed`, one debug line and nothing at normal level; what each handler returns is unchanged.
generation-output: unchanged
generation-affects: none
