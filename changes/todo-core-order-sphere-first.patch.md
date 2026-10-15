**Headline:** The plan now finishes the spherical grid first, then the whole sphere with real stars for every nebula and Wolf-Rayet stars.

### Changed
- Core's order is rebuilt: A the spherical grid (GEN.232, then GEN.230 and GEN.223 in one migration, then GEN.233), B the whole sphere (GEN.234 to GEN.243 and the globular cluster items GEN.158 to GEN.164), C the spherical-region database query.
- GEN.242 (pre-main-sequence state) and GEN.243 (Wolf-Rayet stars and ring nebula class X) are filed; GEN.241 now covers every nebula's required stars, with MAP.182 folded in and retired; DOC.48 and TEST.159 follow steps 14 and 15.
- The early alpha rule: the database schema is not a constraint on any item (TODO.md rule and notes on GEN.223, GEN.230, GEN.233, GEN.235, GEN.237); the scatter from the origin keeps the banded layered approach (GEN.240).

TODO items resolved: MAP.182 (folded into GEN.241) / Added: GEN.242, GEN.243
