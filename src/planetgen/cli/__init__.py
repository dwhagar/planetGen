# planetgen/cli/__init__.py

"""
The command-line tools. `generate.py` is the `planetgen` command (its
subcommands `system`, `sector`, `galaxy`, `plan`, `phenomenon`, `population`,
`check-math`, `fingerprint`, `check-db`, `versions`, `warm` and `benchmark`; `benchmark.py` and
`version_history.py` back two of them, `warm.py` the `warm` one). The others run on their own with
`python3 -m planetgen.cli.<name>`: `orbits` (the positional update),
`worker`, `query`, `migrate`, `reset`, `lockouts`, `lift_lockout`,
`forget_devices`, `reset_two_factor` and `config`, described in `docs/cli.md`.
`stage_progress.py` is the progress bar of the maintenance commands.
"""
