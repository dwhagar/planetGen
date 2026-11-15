**Headline:** An air that is still creeping toward an ice-albedo snowball when its 60 steps run out is given the steps to get there, instead of being stored as unsettled.

### Changed

- GEN.1718 part 2: `surface_air.solve` gives a fixed point that reaches `MAX_STEPS` with its last eight temperature residuals all of one sign (`_creeps`) one more `CREEP_EXTENSION` of 60 steps, with its damping restored to `DAMPING`. A swing early in the solve halves the damping to its floor of 0.05; an ice-albedo creep after it (the answers fall and the ice grows) then moved 5 percent of the residual a step and had not arrived in 60 (key 175 on the first Train 30 head: 236 K to 202 K, unsettled). On the Train 30 head ba4a3515a1, of 500 seeded keys (0 to 499) 14 were stored unsettled (ten of them no steady state, kept by the rule of GEN.1718 part 1, and unchanged); keys 81 and 276 are creeps and now settle (at 187 and 231 K). A solve that settles inside 60 steps, a giant, and an air whose root find found no steady state are not touched: 498 of 500 keys are identical. Keys 128 (a runaway at 950 K) and 194 (273 K) stay unsettled: they are not creeps.
- The new settings, `CREEP_STEPS` and `CREEP_EXTENSION`, are solver settings (GEN.1718), not physical values.

generation-output: changed
generation-affects: fill
