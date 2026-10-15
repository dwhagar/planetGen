### Added

Added `examples/redis/planetgen.conf` (append-only file with `appendfsync everysec`, `noeviction`, a `maxmemory` guard and a loopback bind), and `install.sh` and `update.sh` now warn, with the exact line to add, when the running Redis server's append-only file is off, its memory policy evicts, it has no memory limit or it listens beyond loopback. They never edit your Redis config. The install guides describe it (OPS.48).
