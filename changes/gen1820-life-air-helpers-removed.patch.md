**Headline:** The three life-ledger helpers that only the tests called are removed.

### Removed
- GEN.1820: `life_air.reactions_of`, `takes` and `gives`. The air and the air history read the class rows of the ledger, which already carry what each class draws down and makes. The 35 balanced reactions stay as the ledger `tests/test_life_air.py` holds the class rows to, and the tests read them directly.

generation-output: unchanged
