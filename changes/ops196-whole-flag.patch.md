**Headline:** the `whole` flag of a command is gone: a command that is one unit by nature says so with the break declaration the queue already reads.

### Changed

- OPS.196: `CommandSpec.whole` is deleted. A command whose stage is one unit by nature gives `breaks=ONE_UNIT` (`commands.ONE_UNIT`, which names `breakpoints.whole_job`: one item, no cut), and `CommandSpec.one_unit` is what the unit-length ceiling (PERF.155) reads to leave it alone. The 33 commands that set `whole=True` are the same 33; the cut they get and the ceiling's exemption are unchanged.
