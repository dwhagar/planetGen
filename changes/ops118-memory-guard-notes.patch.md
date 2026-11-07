**Headline:** `install.sh` and `update.sh` now say when a memory guard (earlyoom, systemd-oomd, nohang) is running, and a warm-up that fails prints what the kernel and the guards logged meanwhile; INSTALL.md says how to keep a guard from ending a work unit.

### Changed
- **A work unit of Boss's warm-up (unit 11) ended with status 15 (SIGTERM), and nothing in the program sends that to a work horse (OPS.118).** The usual sender is a memory guard on the host: earlyoom and systemd-oomd end the biggest process with SIGTERM when memory runs low, before the kernel's own killer would end it with SIGKILL, and a work unit of the warm-up is often the biggest. The cause on the test server cannot be named from the program; the next occurrence is in `error.log` (OPS.116) with its signal, unit and run.
- `install.sh` and `update.sh` warn when `earlyoom`, `systemd-oomd` or `nohang` is active and say what a unit ended that way looks like.
- A warm-up stage that fails now prints (and writes to `update.log` and `error.log`) the kernel and memory-guard lines about ending a process from the minutes it ran, when `journalctl` is there, so a status 15 comes with its cause.
- INSTALL.md has a "Memory guards" section next to "Swap": keep swap on, give the guard something to protect (`earlyoom --avoid` for mysqld, redis-server, apache2, httpd), or lower `queue.max_host_workers`.

generation-output: unchanged
