**Headline:** The train gate ends a file's pytest that has gone idle, with a stack dump of every process in the log, instead of waiting out the whole limit.

### Fixed

- A file the gate ran under `pytest -n 2` sometimes left its controller idle for 50 minutes or more after a test passed its 1500 s limit, with no worker running, and the lane lost the time. The gate now runs each file's pytest in a session of its own and looks at it every 10 seconds: when no process of it has used processor time and it has written nothing for 600 seconds, the gate writes the stacks of every process (py-spy) and the wait channel of each into the file's log and ends the run (code 125), so the next stall names what it waits on.
- A file that passes the gate's `--limit` is now ended together with every worker and child process it started (code 124). Before, only the controller was ended and the xdist workers kept running through the next file.
