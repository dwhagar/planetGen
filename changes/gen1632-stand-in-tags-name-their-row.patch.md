**Headline:** Every stand-in tag in the code now names its registry section and matches its row, and the registry lists every stand-in table of the plan.

### Fixed
- The registry now reads 21 tables of the gap-filling plan that it skipped (66 rows), and rows that sat after a blank line or a repeated block (LQ12 to LQ16, GN1 to GN4, OB9 to OB14, LM6) are in their tables; two stale copies (LQ11, CMX7 and CMX8) are gone, and the second LM6 row is the first LM5.
- A tag such as `[Stand-in 6.294 AD1, tier T7, tweak L1.]` carries the section, because the same id is used in several sections; its tier and tweak are the row's. Ids of the nightside collapse tags (NC1 to NC6) are the rows NEB1 to NEB6.
- Rows that named constants that moved (OE9, EC8, CC7, DE4, SC6, GM3) now name the real ones.
- `src/tests/test_stand_in_tags.py` fails when a tag has no section, names no row, states another tier or tweak than its row, or when a plan row falls outside the registry, a section repeats an id, or a row names a constant that does not exist.
- TEST.768: the 19 constants that `test_constants_completeness` found without a tag (the far-zoom lanes settings, the run and dealing settings, the air history steps, the starlight sub-steps) say what they are; `ice_growth` reads the Julian year and the water density from `constants` instead of its own copies.

generation-output: unchanged
