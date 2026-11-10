**Headline:** Generating a system is about 9 percent faster, with every stored value identical.

### Changed

- `interior_heat.heat_flow` now keeps the result of each call it has answered (PERF.198): about two calls in three of a system ask the interior-heat box about the same body again (the settle-to-final-age step and the tide's decay-time loops), and the answer is the same object, so the stored numbers do not change. On 12 systems the run took 15.2 s before and 13.8 s after, with the same hash of every body's values.
- The item as first written (a faster `detmath.matmul` and a cache of equal `gas_equilibrium.solve` calls) found little left to take: `matmul` is 1 to 3 percent of a system since PERF.201 and PERF.246, and equal gas solves repeat in under 1 percent of calls.

generation-output: unchanged
generation-affects: none
