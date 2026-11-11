**Headline:** The light elements in a planet's core now range from 0.02 to 0.25 of its metal, so Mars-like lighter cores appear.

### Changed
- **The core's light-element share is a lognormal, not a flat 0.03 to 0.12 (GEN.993).** Earth's core holds 0.05 to 0.10 of its metal as sulfur, silicon, oxygen, carbon and hydrogen, and Mars's holds 0.07 to 0.20; the flat draw never reached the high half of Mars. The share is now drawn lognormally with a median of 0.08 and a width of 0.25 dex, kept between 0.02 and 0.25. About 14 percent of rocky planets now have a core lighter than 0.15, which lowers their core mass, so a few percent of planets come out a little larger and less dense, and their core heat is a little lower. The stand-in is BC18 in the registry (tier T7, tweak L2, tuned to the Mars total), and a test checks that the draw covers the Earth's range and the Mars total.
generation-output: changed

generation-affects: fill
