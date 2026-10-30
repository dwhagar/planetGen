**Headline:** The sector and system commands now run as stages of the shared run, so no generating command starts a run of its own.

### Changed
- **`planetgen sector` and `planetgen system` follow a stage of the shared run (PERF.202).** The sectors of a `sector` run and the one system of a `system` run were each started by a one-step run of their own; they are now a stage of the shared run like every other generating command, and a test keeps it so (only the one-unit maintenance commands start their own). The output is unchanged. A `sector` or `system` run belongs to no galaxy, so it records no galaxy seed or build.
