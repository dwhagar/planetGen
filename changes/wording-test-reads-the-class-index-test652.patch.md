**Headline:** The wording test that checks a body only claims life when the index gives it life now reads the life-class index, so a world whose life lives in acid clouds is not reported as a lifeless body that claims life.

### Changed

- `tests/test_wording_agrees_with_the_index.py` (TEST.652): `test_a_generated_body_claims_life_only_when_the_index_gives_it_life` asked for `phi_bio > 0` wherever the page had a life or timeline sentence. `phi_bio` is the water-and-carbon index (Class 1a), 0 for a Class 5a world whose life lives in sulfuric-acid clouds aloft (`acid_medium` level 1, which needs no surface temperature), so a hot vapor world (surface 2595 K, seed 13 of the test) with Class 5a life, small mobile animals and a stored timeline failed. The body is not lifeless: it has a stored class that the class index supports (`life_class_scores`, `acid_medium` cloud layer). The test now asks that the text names a stored class and that the class index supports it.

generation-output: unchanged
