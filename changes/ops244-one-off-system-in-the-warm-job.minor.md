**Headline:** The one-off system page makes its system inside the supervisor's warm light-job process instead of starting a command-line process, so it needs no process, imports or run of its own.

### Changed
- **The one-off system page (OPS.244, slice A).** `system_page.run_generator` hands the page's options to a light job (`api_jobs.call_and_wait`, new) that calls `system_page.generate_in_process` in the process the supervisor forked from its warm template: the generator's own parser reads and checks the options (`generate.process_args` takes the words as an argument now), the system is drawn by the piece a `system` unit draws it with (`run_system.one_off`, new; the page text is `_output_text`, shared with the unit's finish piece), and a rogue planet or comet by `run_phenomenon`, which makes no run. The page's output is the generator's warning, what came of a class request and, with Debug on, the lines the process logged (the areas asked for kept). Before, each request started `planetgen system --output` as a process of its own (28.7 s of start-up and a run for one system), and on the way to a ledger run for that command (OPS.271) it would have started a supervisor and a template per request. A page that finds no queue (no Redis, or an update running) says so instead of running the generator in the web process.

generation-output: unchanged
generation-affects: none
