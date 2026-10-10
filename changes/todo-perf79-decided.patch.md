### Changed
- The TODO list records Boss's decision on PERF.79: the scatter stays one job on the Redis queue, never bypassing it, with parallelism only where it measurably helps.
