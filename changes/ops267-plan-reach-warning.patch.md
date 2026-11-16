**Headline:** `planetgen plan --reach-shells 3` no longer ends with an error when the galaxy reaches past the shells it was allowed to walk.

### Fixed

- OPS.267: when the outline's walk was cut short and no `--edge` was given (`--reach-shells 3` alone), the plan's closing warning tried to write the Edge dial's value and failed with a TypeError. It now says the plane still qualified at the reach (the shell and `--reach-shells`) and to re-run with a larger one; with `--edge` given it says what it said before.
