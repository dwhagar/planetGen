**Headline:** The stored mean molar mass, air density and scale height of a world with trace gases now include the trace gases (they left them out by parts in a hundred million).

### Changed

generation-output: changed

- `trace_air.apply` adds the trace gases (nitrous oxide, phosphine, nitrogen dioxide, methyl chloride) to the air after the air was solved and stored its mean molar mass, density and scale height, and did not make those three again, so they described the air without the traces while the partial pressures and the pressure described the air with them (TEST.616: the sweep's rules for them missed by 1e-9 to 4e-8 on 3 to 4 bodies of 361). It now makes the three again from the air it leaves, at the stored temperature and the stored pressure, which do not move.
- Stored values change by up to a part in 1e7 on the worlds that have trace gases, so this is a generator-output change: the epoch bump and golden rewrite are the Train Manager's.
