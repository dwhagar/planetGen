**Headline:** A burst of map requests can no longer start over a hundred queue workers and run the server out of memory.

### Fixed
- **PERF.181: the server ran out of memory and the kernel killed MySQL during a New galaxy run.** Every Galaxy Map request that queues a tile asked the worker starter for upkeep workers, and requests arriving together each counted "none starting" before any of them recorded its own start, so each launched the whole set. On 2026-10-24 about 130 workers of 125 MB took all 16 GB of a host with no swap. Only one caller at a time now counts and starts the workers of a queue, and a started worker counts as starting for 90 seconds (it was 15) so a host that is slow to load a worker is not given a second one for the same need.

### Added
- **`scripts/enable-swap.sh`** makes a swap file (8 GB by default), switches it on and adds it to `/etc/fstab` so it stays on after a reboot; it does nothing when swap is already on. `install.sh` and `update.sh` warn when the host has no swap, and INSTALL.md has a Swap section and a requirements row.
