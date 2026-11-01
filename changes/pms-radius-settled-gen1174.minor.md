**Headline:** A star counts as pre-main-sequence until its radius has settled, as the BHAC15 tracks have it, not until MIST's later nuclear arrival.

### Changed
- **The pre-main-sequence clock now follows the radius (GEN.1174, Boss: both times, class by radius).** The age at which a star reaches the main sequence was MIST's nuclear arrival, 3 to 4 times longer than the age at which the radius settles for stars of 0.1 to 0.3 solar masses (0.2 solar masses: 727 Myr against 188). A star is now pre-main-sequence, and its contraction track runs, until the radius has settled to within 5 percent in the BHAC15 tracks (402 Myr at 0.1 solar masses, 130 at 0.3, 50 at 0.8; the model's own time scaled to the 0.8 value above), `phase_end_age_gy` of such a star is that age, and MIST's arrival is kept beside it as `pre_main_sequence.nuclear_arrival_yr` (and `nuclear_arrival_gy` of the state). The step of 718 to 1358 Myr between 0.105 and 0.11 solar masses (GEN.1173) goes with MIST's clock. The share of living stars that are pre-main-sequence falls to about a third (5.8 to 1.9 percent in the generator's own sampler). Registry row LM6.

generation-output: changed
