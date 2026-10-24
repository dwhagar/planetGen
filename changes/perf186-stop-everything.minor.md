**Headline:** A real Stop everything: one button, one command and one update option end every run, queue and worker process on the server and keep it quiet until you resume.

### Added
- **PERF.186: Stop everything.** The admin Queue page has a Stop everything button (and `planetgen queue stop-all` the same on the command line). It cancels the Generate page's job and every live run, empties every queue, and ends every planetgen worker, job runner and command process on the server together with the processes they started. It also writes a hold that survives a restart of Redis, the web server and the host: while it is on, no worker, run, Generate job, orbital update or cache warm-up starts, the map queues no background work, and a banner on every admin page says so. Resume (the button, or `planetgen queue resume`) lifts it. It works when Redis is not answering.
- **PERF.186: `./update.sh --kill-all`.** Runs the same stop and clear before the update, only when you give it. It leaves no hold behind.

### Fixed
- **PERF.186: Cancel left work running.** A job whose runner had died (killed by the memory killer, a restart) showed as interrupted and Cancel did nothing, while its step and the processes of that step went on using the host. Cancel now ends them, cancels the job's runs and frees the job lock. A running step that ignores the request to stop (it caught the signal or hung in a database call) is now killed after 20 seconds.
- **PERF.185: the process list was read cut off to the terminal's width** where the terminal sets one, so a worker whose command line was long could be missed by the host cap. It is now read in full, here and in the memory watch.
