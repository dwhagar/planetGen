### Fixed
- The disk-space check and the admin Stats tile measure the drive that actually holds the database's data directory (asked of the server with `SELECT @@datadir`, symlinks and mounts resolved), not the boot drive. They show which path and drive were measured, and say "unknown" when a remote server's disk can't be reached.
