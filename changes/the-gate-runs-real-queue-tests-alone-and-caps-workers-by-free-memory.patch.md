**Headline:** The test gate runs the tests that start real queue workers one at a time and cuts its worker count to the memory and idle cores the machine has.

### Changed
- **`scripts/test.py` runs three phases instead of two (TEST.344).** The tests that need no browser and no real queue worker run on every core; the `real_queue` tests, which start real worker processes and read the shared queues by name (so two at once kill each other's workers, TEST.302), run on one worker; the browser tests run last on two. The phase lines show the new phase.
- **The parallel phases are capped by free memory and load.** At the start of a phase the runner counts one worker per 1000 MB of available memory (browser workers 2500 MB) and, in the first phase, no more than the cores the load leaves idle, never fewer than one, and says so on the phase line when it cuts. A lane machine that runs other lanes' gates now gets fewer workers instead of timeouts and OOM kills. `--worker-mb`, `--browser-worker-mb` and `--no-worker-cap` tune or turn it off.
TODO items resolved: TEST.344 (its measuring half, which tests fail only in a parallel phase, stays with the load-flake items it names). Added: none.
