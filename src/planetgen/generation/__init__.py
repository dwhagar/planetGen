"""Generation: building the objects a run creates, and the commands that run it.

Objects: `star`, `binary`, `wide_binary`, `planet`, `belt`, `comet`, `life`,
`evolution` and `system` (a whole star system). The galaxy: `bright_stars`,
`star_population`, `star_labels`, `luminosity_floor`, `phenomenon_scatter`
and `object_first` (the scatters), `directives`, `prevalence` and `config`
(what a run asks for). The commands: `run_system`, `run_sector`, `run_galaxy`,
`run_plan`, `run_phenomenon` and `run_population`, with `run_common`, `stages`
and `steps` (progress and shared rules). Checks and estimates: `limits`,
`validation`, `plausibility`, `phenomena_plausibility`, `stats` and
`object_timing`."""
