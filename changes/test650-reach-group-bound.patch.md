**Headline:** The reach-group density test no longer demands over 1000 expected stars in its smallest shell, and passes again.

### Changed

- `test_gen_reach_groups.py` (TEST.650): the test of the expected stars beyond the outline required more than 1,000 expected stars in the test galaxy's shell 3 and got 873.8 on main. The count is the test galaxy's, not a model fall: the thin disk's infall (GEN.1242, 3,058 down to 310 stars; GEN.1246 back up to 1,205) and then the American spelling sweep (DOC.103, which renamed the nuclear disc population and so moved the shape's random draw, 1,205 down to 874) changed it. The estimate matches the density to 0.04 percent (873.5 against 873.8), which is what the test is for. The bound now only guards against a shell too empty to compare (300).

generation-output: unchanged
