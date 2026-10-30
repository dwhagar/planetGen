**Headline:** The planet architecture now keeps a record of every body it loses and says which ones end thrown out of the system, which collide and which stay bound.

### Added
- `planet_architecture.losses` and `ejection_yield` (GEN.1119, RP1): each body the architecture loses (planets a moved giant sweeps or a crowded giant throws out, chain bodies below the planet minimum, the mass a merge or the core cap does not keep) is shared among ejected, collided and bound by the Safronov number of the encounter (`scatter.loss_fates`). The ejected ones carry an ejection record (mass, core, orbit, perturber, Hill radius for the moons that could stay, time, formed beyond the snow line) for the rogue builder (GEN.1111). A star's planets are exactly what they were; `Architecture` gains a `lost` field. Registry rows EJ1 and EJ2 (section 6.102). On 3,000 sun-like seeded stars with a giant, about 2.5 planets per giant-hosting star are expected to end ejected, the scale Barclay et al. 2017 find; over all stars it is about 0.15 per star.

generation-output: unchanged
