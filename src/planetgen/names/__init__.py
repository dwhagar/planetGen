"""Names: how generated things are named.

Stars and sectors get word-salad names (`wordsalad`, `wordlists`); their
uniqueness rules are `uniqueness`. Planets, moons and belts take their names
from their star (`bodies`). Objects with no star get a name from their ID:
`object_id` (the 76-bit position ID of an interstellar object), `naming_key`
(the galaxy's naming key) and `gated_phoneme_codec` (the codec that turns an ID
and the key into a pronounceable name)."""
