**Headline:** Two nebula field slabs drawn at once no longer deadlock: a slab takes its rows back only when it has rows to take (OPS.291).

### Fixed

- OPS.291 (bug): stage 17 of the OPS.247 full run (Nebula field, 184 units, 2 workers) failed 29 units for good with MySQL error 1213, Deadlock found. The OPS.264 line at the top of a slab unit, `DELETE FROM nebula_field WHERE cell_i = ?`, matches no row the first time a slab is drawn, and under REPEATABLE READ a DELETE that matches none still gap-locks the end of the key; two slabs in two workers each held that lock when the other came to its first INSERT. A slab now checks with a plain read (no lock) whether it has any row and runs the DELETE only then, so a slab drawn again after a try that wrote it still takes its rows back and a first draw takes no lock at all. One slab owns each `cell_i`, so the check cannot race another unit. A test draws two slabs in two connections with a barrier between the take-back and the first row; it fails with 1213 on the old code.

generation-output: unchanged
generation-affects: none
