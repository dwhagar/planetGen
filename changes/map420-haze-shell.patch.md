**Headline:** A hazy body, Titan-like, now has a coloured haze shell past its limb, as deep as its haze.

### Added
- MAP.420: a body with a stored haze depth (`tau_haze`) of 0.3 or more has a haze shell in the haze's own colour. The look record carries `hazeReach` (0.04 to 0.14 radii past the limb) and `hazeStrength` (0.25 to 0.8), both growing with the haze depth, and the System Map marker, the 3D system view and the body picture each draw a second glow shell in the haze colour with the glow shader they already have. Nothing new is stored, and a body with a thin haze, dust only, or no air has no shell.

generation-output: unchanged
generation-affects: none
