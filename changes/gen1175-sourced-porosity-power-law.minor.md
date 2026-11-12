**Headline:** The roundness engine's off-by-default porosity rule is now the sourced carbonate-rock power law, not a recalled ceramics exponential; the rule in use is still "none", so nothing moves.

### Changed

- GEN.1175: `roundness.POROSITY_RULES` loses `exponential` and gains `power`: strength times `(1 - phi / 0.5)^1.67` ([R0083], 253 uniaxial tests of carbonate rock and gypsum, 0.4 to 50 percent porosity; relative strength 0.69, 0.43, 0.22, 0.07 at porosity 0.1 to 0.4), with the argument held at `POROSITY_FLOOR` 0.1 (stand-in RN12, tier T3, tweak L1; the floor T9, L3). Right form, wrong material (no ice-rich law exists), so the rule stays `none` in `roundness_strengths.json`; switching it on is a separate, lock decision after the law is judged on Proteus.

generation-output: unchanged
