**Headline:** The gas equilibrium solver no longer gives up on near-stoichiometric water with a trace of argon.

### Fixed
- **GEN.1003: the equilibrium solver could stall with "no convergence" on a mixture of hydrogen and oxygen at the proportion of water with a trace noble gas.** The last row of its scaled Newton system was rounding noise, so the step to a trace element jittered for ever. A solve that does not converge now retries with that row scaled by the total; a solve that converges the first time is unchanged.
