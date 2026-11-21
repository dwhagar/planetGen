**Headline:** The wiring burn-down starts: a trace-gas function nothing read is retired and two things the guard listed are marked as by design.

### Removed
- GEN.1559: `mixture_trace.draw_trace_gases`, the mixture-level draw of the trace gases, is gone. It duplicated what `physics.trace_air` does with `draw_trace_amounts`, which the generator reads; only one test called it, and that test now checks `draw_trace_amounts`.

### Changed
- GEN.1559: the wiring register marks `design_flags` (the check that the life-linked gases m3, m4 and m11 stay off unless a flag allows them, run by its test) and the Class M meter's row AT4 as by design, with the reason.

generation-output: unchanged
generation-affects: none
