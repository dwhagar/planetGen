**Headline:** The count of Class M worlds that nothing read is gone from a system, so no count is taken by class letter.

### Changed

- GEN.1400 (the last leftover of the class-last audit, F7): `StarSystem.m_count` and the Class M tally in `count_habitable` are removed; `count_habitable` returns the habitable count alone and the loader, `from_dict` and the tests use that. Nothing generated from `m_count` and no page showed it, so nothing stored or drawn changes.

generation-output: unchanged
generation-affects: none
