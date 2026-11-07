**Headline:** A new galaxy no longer fails at stage 3 because the compiled kernels could not be cached.

### Fixed
- **OPS.112: the galaxy job failed at stages 3 to 14 with "cannot cache function: no locator available".** Since SEC.37 the code tree is read-only to the service user, so numba could not write its cache beside the source. `planetgen` now points `NUMBA_CACHE_DIR` at `/var/cache/planetgen/numba` (made by `create-cache-dir.sh`, so by `./update.sh`), or at a private per-user temporary folder when that is missing, before any kernel loads.
