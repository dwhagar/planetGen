**Headline:** A planet with a heavy moon no longer takes a minute and a half to generate: the moon orbit solver gives up on a stage that is not finishing and moves to a looser one.

### Changed
- **The moon orbit solver's ladder is capped (PERF.231).** A rocky planet with a moon of a few percent of its own mass (a Mars-sized planet with a moon a thirtieth of its mass, and a thermal tide) made the spin and orbit equations so stiff that the second integrator, at a tolerance of one part in a million, ran 400,000 evaluations (100 seconds) and then failed, so the class census timed out in 8 cells. Each stage now has a few tens of thousands of evaluations: the same two stages as before, then the first stage's own tolerance (one part in a hundred thousand), then one part in ten thousand. The worst draw in the census cells is under 13 seconds (it was past 100), a planet that was solved by the second integrator late now gets a slightly looser answer (a few parts in a hundred thousand), and one no stage finishes keeps its moons as born, as it did after the long wait.

generation-output: changed

TODO items resolved: PERF.231.
