### Changed

- The Generate page's size and time estimate and the one-off system page's generator now run on the Redis queue and the page waits for them (PERF.24, step 4d). Without a Redis server (Windows without WSL's Redis) they run in the web process as before.
