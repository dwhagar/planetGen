**Headline:** The rule that tells whether two stars went through a common envelope in the past now exists as a tested module, ready to replace the white-dwarf-only test.

### Changed
- **The history of a pair (GEN.1012, first part).** A new rule (`physics/pair_history`) finds the orbit a pair had when its members were at their heaviest, from the mass they lost, and sets each member's largest radius against its Roche lobe at that orbit's closest approach. A pair whose member filled its lobe or swallowed its partner can only be a post-common-envelope pair, for any member and any tree depth, not for white dwarfs only. Nothing reads it yet; it is wired into `classify_pair` on the generation train with the member record (GEN.1009).

generation-output: unchanged
