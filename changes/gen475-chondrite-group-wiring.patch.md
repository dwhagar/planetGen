**Headline:** Every asteroid belt and field now says which meteorite group its rock is nearest to (CI, CM, H, L and so on, or none for an icy mix).

### Changed

- Belts and asteroid fields gain a "Nearest Meteorite Group" fact, worked out from the mix of grains the body was drawn with (the iron in its olivine, and its metal, refractory and fine-grained shares) and the average of each chondrite group. A mix with 5 percent or more ice or water is called none, since the groups are rocks. The belt's own mix and density are unchanged. Stored as `chondrite_group` on `asteroid_belts` and `asteroid_fields` (schema migration, nothing backfilled: a regenerated galaxy fills it). generation-output: changed
