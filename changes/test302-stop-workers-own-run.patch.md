**Headline:** The benchmark suite stops only its own queue workers, so it no longer kills the workers other tests started on the same machine.

### Changed
- **TEST.302: `benchmark_suite.stop_workers` is scoped to the run's own control database.** It ended every process whose command line held `planetgen.cli.worker`, so a benchmark run beside the queue tests killed their workers and those tests failed with exit -9. It now ends only the workers started with the run's `PLANETGEN_CONTROL_DATABASE` (the scratch one its steps hand on), and a test starts two stand-in workers to show a foreign one survives. Nothing the generator makes changes.

generation-output: unchanged

TODO items resolved: TEST.302.
