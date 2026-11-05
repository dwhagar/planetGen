**Headline:** A census script counts how many keep-out spheres overlap in a generated galaxy.

### Added
- `scripts/census_keep_out_overlap.py` reads a galaxy's stored objects and counts overlapping keep-out spheres, nested pairs, and clusters. The answer decides whether the course planner needs its cluster fallback. An estimate from the model's own numbers says a quarter to a third of the star spheres near the Sun overlap another, so the fallback is needed.

generation-output: unchanged
