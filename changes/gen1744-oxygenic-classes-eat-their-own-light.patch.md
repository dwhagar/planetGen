**Headline:** The blue and far-red oxygenic classes now make oxygen from the light of their own band, not from the 400 to 700 nm of class 1a.

### Changed

- GEN.1744 (every galaxy): `biogenic_oxygen.active_light_ratio`, `productivity` and `refuge_productivity` take the `Band` of the oxygenic class (`band_of`: its `band_nm` and `efficiency` in `life_classes.json`: 1a 400 to 700 nm at 1.0, 1b blue 300 to 500 nm at 1.0, 1c far red 700 to 800 nm at 0.8), where all three used class 1a's 400 to 700 nm. The source `S` of a world is the sum of its oxygenic classes' shares, each counted in its own band (`source_of`, `profile_of(..., producers=...)`), and the oxygen rise (`rise_age_gyr`) takes the same sum. The methane that a burial-fuelled class makes follows its own band's production. Class 1a is unchanged under any star (1.0 under the Sun), so worlds with only class 1a keep their air; a world whose oxygenic life is 1b or 1c now has less or more oxygen by how much of its star's light lies in that band: a hot star is rich in blue, a cool one in far red.

generation-output: changed
generation-affects: fill
