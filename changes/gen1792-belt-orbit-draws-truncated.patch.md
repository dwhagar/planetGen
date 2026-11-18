**Headline:** A belt's bodies are no longer piled on the eccentricity and inclination caps: their orbit draws are truncated, not clipped.

### Fixed
- **A third of a wide belt's bodies sat at exactly e = 0.5 (and the widest at exactly 60 degrees).** `belt_roster.entry_of` cut each Rayleigh draw with `min(draw, cap)`, which moved the whole tail onto the cap (37 percent of a belt three times as wide). The draw is now truncated by its inverse distribution (`trojan_swarm.truncated_rayleigh`, reused; one draw per body, no extra stream numbers), so no body is piled on a cap and the shape below it is the Rayleigh. Registry row 6.182 MR3 says truncated (GEN.1792).
- Every stored belt member's eccentricity and inclination change slightly for the same seed; the forced eccentricity sum of `belt.py` is unchanged.

generation-output: changed
generation-affects: fill
