# The dynamic queue chain: jobs declare, the engine splits (OPS.161)

Status: skeleton and profile reader built (declaration, engine, reader of the OPS.155 rows, logging); no stage uses it yet. Asked by Boss, 2026-11-12, 3:00 PM EST ("a full dynamic queue chain where the queue engine will decide the split based on resources profiles and suggested break points by the job entering the queue"), with his rulings of 3:01 to 3:04 PM EST. It is the second half of [self-tuning-jobs.md](self-tuning-jobs.md) (DOC.120: what is recorded and learned) and uses the plan-then-execute job model of [job-graph-and-memory-broker.md](job-graph-and-memory-broker.md) (DOC.113).

## What a job declares

Every job entering the queue declares three things, as plain data (`queue/breakpoints.py`, `Declaration`, named by `CommandSpec.breaks`):

1. **Suggested break points.** The cut the job would make today (`suggested`), the item positions a cut may fall at (`boundaries`, any position when unset) and spans of items that must stay in one unit (`whole`). What is left is `legal_cuts`. The typed split of a stage stays: it is the default, the clamp (a unit holds a quarter to four times the suggested unit unless the job sets `floor_items` and `ceiling_items`), and the cut whenever the engine has nothing to go on.
2. **A split flag** (Boss, 3:02 PM EST: "dynamic split is byte identical, or use the suggested split in the job information"). `byte-identical` says the job's output is the same at any legal cut, so the engine may choose the cut; `suggested` says the output depends on the cut, so the engine uses the suggested cut exactly. A job that declares nothing is `suggested`. PERF.260's test holds a job flagged `byte-identical` to identical output at two cut sizes; until a job has passed it, it stays `suggested`.
3. **Its dependencies** (Boss, 3:01 PM EST: "what depends on what"). `after` are the kinds of command that must be done first (a command that names none takes `commands.REQUIRES`, the edges of PERF.247, so nothing is declared twice); `reads` and `writes` are the tables, files and named facts the job uses and changes. `breakpoints.may_overlap` says two jobs may run together when neither writes what the other reads or writes. A job that declares neither is unknown, not free: it runs in plan order, as today. Whether a pair really is independent is still proven byte for byte by OPS.131/OPS.132 before it overlaps; the declaration is what that proof starts from. A dependency inside a job (one item needs the one before) is a span that stays whole, so the engine never cuts across it.

A command that is `whole` (one unit by nature) declares itself: one item, no cut. A command that declares nothing is `None` in `commands.declaration` and keeps its typed split. The migration of the registry to declarations goes with PERF.261, one stage at a time.

## What the engine does

`queue/splitter.py`, `choose(kind, declaration, reading, resources)`, is pure:

| Case | Cut | Reason logged |
|---|---|---|
| flag `suggested` (or nothing declared) | the suggested cut | its output depends on the cut |
| flag `byte-identical`, no usable reading (fewer than 3 samples, no time per item, no trust) | the suggested cut | the profile has nothing usable for the kind yet |
| flag `byte-identical`, usable reading | units of about `TARGET_SHARE` (a fifth) of `queue.unit_ceiling_seconds` at the profile's seconds per item; fewer items if `mb_base + mb_per_item x items` would not fit `MEMORY_SHARE` (a half) of the host's free memory over the workers the stage may have; kept between the job's floor and ceiling; at least one unit per worker where the floor allows; at most `MAX_UNITS` units; every cut snapped to a legal position (the largest at or before the aim, else the next after it) | which of time, memory, floor, ceiling or workers limited the unit |

The decision (`Decision`) carries the flag, the old and new unit counts and sizes, the samples behind it, the figures it was made from and the limit that bound. The same inputs give the same cut.

**The cut is fixed at plan time.** `commands.Expansion.split` is the decision an `expand` function made with `splitter.decide(kind, declaration)`; `wrapper._open` writes it into the stage's plan (`plan["split"]`, with the chosen cuts) and into the run's log (`kind=split.decided`: the job, the flag, `by`, the old and new units and items, the samples and the limit as fields) before the units are queued. A restarted unit reads the stored plan, so it has the cut it started with; the planner never cuts a stage twice.

**Everything is logged** (Boss, 3:03 PM EST: "ALL of this is logged including decisions when the queue changes values"): a stage whose engine kept the suggested split says so with the reason, and a stage whose split changed says what it changed from and to and why. The Queue page shows the run's log, so the recent decisions are on it; OPS.158 adds the per-kind table.

## The profile contract (OPS.155)

`queue/splitprofile.py` reads the `job_profiles` rows of `jobprofiles` (`unit_rows`, kind `unit:<command>`, finished units with seconds per item). `lookup(kind)` gives a `Reading(seconds_per_item, mb_base, mb_per_item, samples, trust)`: the median seconds per item, and a memory line through the lowest and highest size bucket's median peaks (one bucket gives the raw peak over the items and no base, which overstates the slope and so never cuts a unit too large). Rows from this machine's shape and epoch (the rule of `profileclaim._shape_of`, reused) count in full at 3 samples or more; otherwise rows from elsewhere count at half trust (0.5); fewer than 3 of either gives `None`, so the job keeps its suggested split. `PLANETGEN_SELF_TUNING=0` or `PLANETGEN_JOB_PROFILES=0` (the test default) turns the reading off, and any failure to read leaves the suggested split. The same table feeds the gate claim (OPS.156) and the worker caps (OPS.157).

## What is not built

* Declarations for the registered commands and the use of `splitter.decide` in their `expand` functions: PERF.261, starting with the Stage 17 uncharted index (`uncharted_index.shell_runs`), then block lights, link, settle and the warm-up chunks; each behind PERF.260's byte-identical test.
* The overlap of declared independent jobs (OPS.131/OPS.132): the declaration carries the data (`reads`, `writes`, `after`); the executor does not use it yet.
* A learned figure that moves at most 25 percent per run (the design note's rule) applies to the stored profile, not to this engine.
