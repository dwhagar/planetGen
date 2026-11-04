**Headline:** Stored values, API keys and names that were spelled the British way are now spelled the American way.

api-version: changed

### Changed

- The values and keys the site stores and the API returns use American spelling, with no alias for the old words: `center`, `center_pc` and `center_distance_pc` (were `centre...`), `vapor` (was `vapour`), `nuclear_disk` (was `nuclear_disc`), the states `canceling` and `canceled`, the column `ground_oxidized_fraction`, and the belt origin sentences (`asteroid_belts.origin_sentence`, `origin_history`: "the disk past ..."). A migration rewrites the stored rows, and the API version is 13 (DOC.103).
- Names in code, commands, tests and the docs follow the same spelling: disk, flavor, ionizing, oxidized, link-neighbors. The spelling test covers the new forms.
