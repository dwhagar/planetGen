# planetgen/datasets/__init__.py

"""
Static data the generator reads from files fetched at setup (GEN.548).

A dataset is a small set of files from a published source (isochrone grids, white-dwarf layers, an equation of
state, line lists). They are not kept in the repository. `update.sh` and `install.sh` run `python -m
planetgen.cli.datasets fetch` once, which downloads each file with the standard User-Agent (API.31), checks its
sha256 and keeps it in a cache folder outside the checkout (`fetch.cache_dir`). Nothing is fetched while the site or
the generator runs: a reader calls `fetch.locate`, which only looks at the folder, and gets either the files or a
labelled stand-in (docs/design/gap-filling-plan.md 6.5 and the registry of section 6) when they are missing.
"""
