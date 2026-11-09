**Headline:** A cold giant layer whose last oxygen sits in lithium hydroxide no longer stops generation with "the ice caps did not settle".

### Changed

- `mixture_giant.cap_ices` cycled on a layer where water, methane and ammonia had all condensed and 1e-37 of oxygen per hydrogen was left, held mostly by lithium as hydroxide at the chemistry's temperature: water was over its cap by 40 orders, each round took a sliver of it and the equilibrium made it again, for 1000 rounds, then `MixtureError`. After the rounds, if every species still over its cap is one the budget cannot make more of than `MIN_CONDENSED` per hydrogen atom (the amount a condensate is dropped from the record at), the layer is returned as it stands; a cycle with more to condense still raises. One system in 700 of a default star mix hit it (seed 533 of a G2V system with PLANETS and MAX_PLANETS; the same error stopped a `test_db_column_round_trip` run). The rounds and arithmetic before the 1000th are unchanged, so every layer that settled before settles the same.

generation-output: unchanged
