**Headline:** The three spreads of a trojan swarm (how widely its bodies librate, how eccentric and how tilted they are) are now fitted to the Jupiter and Neptune Trojans, with the values proposed; no constant has changed yet.

### Changed

- GEN.1643 (slice 0): new `scripts/fit_trojan_scales.py` turns the published Trojan laws into the Rayleigh scale the generator draws from, for the stand-ins TJ9 (libration amplitude), TJ10 (eccentricity) and TJ11 (inclination) of section 6.121 of the gap-filling plan. The sources are the debiased Jupiter Trojan laws of the proper amplitudes and of the libration integral (the Jupiter Trojan paper, [R5245]) and the Rayleigh widths of the Neptune Trojans [R5244]. The proposal, recorded in the registry rows and in `reports/gen1643-trojan-scales-2026-11-12.md`: libration 13 degrees (now 10), eccentricity 0.06 and inclination 12 degrees stay (the fit gives 0.058 and 11.7). Slice 1 sets them and lowers the rows to L1 or L2.
- The same fit shows that `trojans.TADPOLE_LIMIT_DEG` (23.9) is the nearest a tadpole comes to the secondary, not a libration amplitude (that is 78 degrees at the edge of the zone), and that about a fifth of the real Jupiter Trojans librate wider than 23.9. The note is on the TJ3 row; the correction rides with slice 1.
- No constant of the generator is changed and nothing is stored.

generation-output: unchanged
