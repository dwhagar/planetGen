**Headline:** Stored values, API keys and names that were spelled the British way are now spelled the American way.

### Changed
- **Stored words and API keys (DOC.103, schema v210).** The values the site stores and the API returns use American spelling, with no alias for the old words: `center`, `center_pc` and `center_distance_pc` (were `centre...`), `vapor` (was `vapour`), `nuclear_disk` (was `nuclear_disc`), the states `canceling` and `canceled`, the column `ground_oxidized_fraction`, and the belt origin sentences (`asteroid_belts.origin_sentence`, `origin_history`: "the disk past ..."). A migration rewrites the stored rows.
- **Names in code, commands, tests and the docs.** They follow the same spelling: disk, flavor, ionizing, oxidized, link-neighbors. `scripts/american_spelling_sweep.py` (`apply`, `check`) redoes the rewrite on any tree, and the spelling test fails if a British word is left anywhere in it.
- **API version 13.** The keys and values named above changed.

api-version: changed
generation-output: changed
