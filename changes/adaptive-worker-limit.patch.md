**Headline:** A generation run lowers its bulk worker count when jobs wait for memory and raises it again when memory is free, instead of leaving jobs to wait 300 s at the memory gate.

### Changed

- Adaptive worker control (`queue/adaptive.py`): a bulk job that has waited 30 s at the memory gate for memory, or a host whose available memory is under three quarters of the gate's reserve, lowers the bulk worker limit by one (never below 1); the workers that wait at the gate leave first and put their jobs back at the front of the queue. After 120 s with no memory wait and room for one more worker, the limit goes up by one. Each change is logged and kept in Redis. See `redis-and-the-work-queue.md` section 5.17.

generation-output: unchanged
generation-affects: none
