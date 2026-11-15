**Headline:** An air that is still creeping toward an ice-albedo snowball when its 60 steps run out is given the steps to get there, instead of being stored as unsettled.

### Changed

- GEN.1718 part 2: `surface_air.solve` gives a fixed point that reaches `MAX_STEPS` with its last eight temperature residuals all of one sign (`_creeps`) one more `CREEP_EXTENSION` of 60 steps, with its damping restored to `DAMPING`. A swing early in the solve halves the damping to its floor of 0.05; an ice-albedo creep after it (the answers fall and the ice grows) then moved 5 percent of the residual a step and had not arrived in 60 (key 175: 236 K to 202 K, unsettled). Of 400 seeded keys on Train 30, 12 were stored unsettled, now 9 (key 156 is not extended: it has no steady state and keeps the rule of GEN.1718 part 1); keys 175, 242 and 275 settle, at 190, 230 and 229 K. A solve that settles inside 60 steps, a giant, and an air whose root find found no steady state are not touched: 397 of 400 keys are identical.
- The new settings, `CREEP_STEPS` and `CREEP_EXTENSION`, are solver settings (GEN.1718), not physical values.

generation-output: changed
generation-affects: fill
