### Added
- The ID-to-words codec (`gatedPhonemeCodec.py`, which sat in the repo root and was used by nothing) is now `planetgen.names.gated_phoneme_codec`, with tests that pin its words. It needs only the standard library, and its decoding is exact when the identifier's length is given (GEN.120).
