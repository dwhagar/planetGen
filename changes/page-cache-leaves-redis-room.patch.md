**Headline:** The page cache can no longer fill Redis and stop a galaxy run; the Redis memory limit is now 4 GB.

### Changed

- The shared page cache stops adding answers to Redis once Redis is 70 percent full, so the job queue and run records always have room. (A 512 MB Redis was full of cached pages, about 365 MB of them, when the phenomena scatter died.)
- `./update.sh` now warns when the Redis `maxmemory` is under 2 GB and `examples/redis/planetgen.conf` ships `maxmemory 4gb`.
- `planetgen reset` drops that database's cached pages from Redis.
