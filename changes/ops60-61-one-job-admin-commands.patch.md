**Headline:** The admin commands each do one job: listing lockouts, lifting one, forgetting devices and resetting two-factor are four commands, and two obsolete one-off tools are gone.

### Changed
- **`cli.lockouts` only lists** (OPS.60). Lifting a lockout is `python3 -m planetgen.cli.lift_lockout` with `--ip`, `--user` or `--all`; revoking an admin's trusted-device cookies is `python3 -m planetgen.cli.forget_devices USER`; turning off two-factor is `python3 -m planetgen.cli.reset_two_factor USER`. The old `--ip`, `--user`, `--all`, `--forget-devices` and `--reset-two-factor` options of `cli.lockouts` are gone, with no shim.
- **`cli.version_history` only records** the version key (OPS.60). Its `--list` is removed; `planetgen versions` lists the history.

### Removed
- **`cli.render_parity` and `cli.dedupe`** (OPS.61), the one-off check for schema v29 and the name backfill for v24, which migrated long ago, with their tests.
