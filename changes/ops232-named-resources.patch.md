**Headline:** a job can now hold a named resource besides memory, so two heavy writers or two builds of one warm set never run at once, and the unit that runs alone is a claim like any other.

### Added

- OPS.232: `queue/tokens.py` and `Claim.tokens`: the gate lets a claim go only when no claim held holds one of its names (`alone`, `mysql_write`, `warm:<set>`), and a job that waits is told which job holds the token. A command names its tokens in `CommandSpec.tokens` and the claim of its unit holds them. The stage 0 simulator checks that no token is held twice and that the job that runs alone has the bulk and maintenance lanes to itself.

### Changed

- OPS.232: the claim of the restarted unit that runs alone (OPS.145) holds the `alone` token instead of a `solo` flag; the gate's rows carry `tokens` in place of `solo`. `pause.hold_solo` stays as the way workers on main hold the queue until their lanes move onto the supervisor (OPS.228 and OPS.229).
