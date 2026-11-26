**Headline:** The zone mix test counts planets, not the belt members promoted to planets

### Fixed

- TEST.838: `test_zone_mix_is_spread` counted the ten heaviest round bodies of every belt that GEN.1510 (11-11) promotes to planets (they keep `belt_rank`, all cold): 553 of 837 counted bodies in a 60-system sample, and a hot share of 0.058 on main (218/3736). The test now skips them, as it skips asteroids. No product change, no band change. The hot share is still 0.181 (207/1143) against the floor of 0.20, which is a generation question (GEN.1858), not part of this fix.

generation-output: unchanged
