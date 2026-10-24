**Headline:** A world's four life scores, Microbial, Flora, Fauna and Sentient, now have their own symbols with a number from 0 to 9.

### Added
- **The life score symbols (UX.101).** Four shapes in the icon sprite (a cluster of cells, a stem with two leaves, a paw print, a head and shoulders), each drawn with a small number from 0 (none) to 9 (abundant) and a label that says the same in words. `habitability_life(name, score)` draws one and `habitability_symbol_row(body)` ends a body's row with the cluster when it has `life_microbial`, `life_flora`, `life_fauna` and `life_sentient`. The habitability levels page explains them. The scores themselves come from the model (the habitability model lane).
