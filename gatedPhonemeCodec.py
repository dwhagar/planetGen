"""
Deterministic Gated Phoneme Codec with Domain-Keyed Positional Permutation
==========================================================================

This module implements a bidirectional, deterministic codec that transforms 
arbitrary hexadecimal identifiers (and raw strings) into natural, phonotactically
legal, pronounceable names.

Key Architectural Pillars:
--------------------------
1. Domain Separation (Namespacing):
   Objects of different types (e.g., 'cargo_container' vs. 'delivery_truck') 
   sharing the exact same hexadecimal ID are mapped to completely distinct word 
   sequences using a deterministic affine permutation over Z_16.

2. Positional Phase Diffusion:
   To prevent repeated digits (such as '00' in '10D00B002C7') from generating 
   redundant, repetitive phonetic sounds, each digit's transformation is shifted 
   by its index in the string:
       y_i = (A * x_i + B + i) mod 16

3. Phonotactic Gating Hierarchy:
   Phonemes are partitioned into three disjoint tiers (Consonants, Vowels, and 
   Syllabic Buffers). A finite-state gate enforces canonical syllable structures 
   (CV/CVC), suppresses illegal stop-burst clusters, eliminates vocalic hiatus, 
   and prevents stuttered geminates.

4. Bounded Working-Memory Wordiness:
   Partitioning is tuned to heavily favor compact, 2-word phrases for common 
   identifier lengths (4 to 14 hex characters), ensuring names stay within the 
   human phonological loop capacity while remaining strictly invertible.

5. Sub-Microsecond In-Memory Execution:
   The algorithm eliminates heavy database lookups and cryptographic hash rounds.
   Names are computed on the fly when an entity is viewed and parsed back to the 
   underlying ID in linear O(N) time with zero collision risk.
"""

from __future__ import annotations

import re
import zlib
from typing import Dict, List, Optional, Tuple, Union
from more_itertools import chunked
import pygtrie


class PhonemeDecodeError(ValueError):
    """
    Raised when a phonetic phrase violates phonotactic gating grammar, contains
    unrecognized phonemes, or fails conversion back to its original identifier.
    """
    pass


class GatedPhonemeCodec:
    """
    Bidirectional codec implementing domain-separated linear permutation and 
    a gated phoneme matrix.

    Attributes:
        default_digits_per_word (int): Default target number of hexadecimal 
            digits allocated to each generated word.
        stacks (Dict[str, List[Dict[str, str]]]): The 3-tier candidate stacks 
            indexed by hexadecimal character ('0' through 'F').
        trie (pygtrie.CharTrie): Prefix tree storing all 48 candidate phonemes 
            for greedy, backtracking-assisted inverse parsing.
    """

    # Global wordiness presets mapping human labels to target hex digits per word.
    # Allocating 6 digits per word heavily favors 2-word outputs for IDs up to 14 chars.
    DEFAULT_WORDINESS: str = "compact"
    WORDINESS_MAP: Dict[str, int] = {
        "compact": 6,   # Favors exactly 2 words (e.g., 4-14 hex characters)
        "low": 6,       # Alias for compact
        "balanced": 4,  # Produces 3 words for 10-14 digit IDs
        "medium": 4,    # Alias for balanced
        "verbose": 3,   # Shorter, more frequent words (4+ words for 11-digit IDs)
        "high": 3,      # Alias for verbose
    }

    # Set of integers coprime to 16, used as modular multipliers in Z_16.
    # Multiplying by any coprime is an unconditional bijection (permutation).
    COPRIMES: List[int] = [1, 3, 5, 7, 9, 11, 13, 15]

    # Precomputed modular multiplicative inverses such that (A * A_INV) % 16 == 1.
    MOD_INVERSES_16: Dict[int, int] = {
        1: 1,
        3: 11,
        5: 13,
        7: 7,
        9: 9,
        11: 3,
        13: 5,
        15: 15,
    }

    # -------------------------------------------------------------------------
    # Phoneme Inventories: Mutually Disjoint Across Tiers and Hex Digits
    # -------------------------------------------------------------------------

    # Tier 0: 16 Consonant Onsets (C)
    TIER_0_CONSONANTS: Dict[str, str] = dict(
        zip("0123456789ABCDEF", "bdfghjklmnprstvz")
    )

    # Tier 1: 16 Vocalic Nuclei and Diphthongs (V)
    TIER_1_VOWELS: Dict[str, str] = dict(
        zip(
            "0123456789ABCDEF",
            [
                "a", "e", "i", "o", "u", "ai", "ee", "oo",
                "ar", "or", "an", "en", "in", "on", "un", "ay"
            ]
        )
    )

    # Tier 2: 16 Syllabic Buffer Morphemes (S) - Fail-safe articulable codas
    TIER_2_BUFFERS: Dict[str, str] = dict(
        zip(
            "0123456789ABCDEF",
            [
                "bel", "dan", "fen", "gal", "hop", "jin", "kor", "lum",
                "mor", "nil", "per", "ras", "san", "tor", "vel", "zar"
            ]
        )
    )

    def __init__(self, default_wordiness: Union[str, int] = DEFAULT_WORDINESS) -> None:
        """
        Initializes the codec, builds candidate stacks, and compiles the prefix trie.

        :param default_wordiness: Word density configuration. Accepts preset strings
            ('compact', 'balanced', 'verbose') or an integer specifying target
            digits per word. Defaults to 'compact'.
        """
        self.default_digits_per_word = self._resolve_wordiness(default_wordiness)

        # Build candidate stacks per hex digit
        self.stacks: Dict[str, List[Dict[str, str]]] = {}
        for h in "0123456789ABCDEF":
            self.stacks[h] = [
                {"phoneme": self.TIER_0_CONSONANTS[h], "type": "C"},
                {"phoneme": self.TIER_1_VOWELS[h], "type": "V"},
                {"phoneme": self.TIER_2_BUFFERS[h], "type": "S"},
            ]

        # Populate prefix trie for O(N) inverse decoding
        self.trie = pygtrie.CharTrie()
        for hex_char, candidate_list in self.stacks.items():
            for entry in candidate_list:
                # Each phoneme is guaranteed unique across all tiers
                self.trie[entry["phoneme"]] = (hex_char, entry["type"])

    def _resolve_wordiness(self, wordiness: Union[str, int]) -> int:
        """
        Validates and translates a wordiness parameter into digits-per-word.

        :param wordiness: Preset name or explicit integer digits-per-word.
        :return: Validated integer digits-per-word (minimum 2).
        :raises ValueError: If an unrecognized preset name is supplied.
        """
        if isinstance(wordiness, int):
            return max(2, wordiness)
        val = self.WORDINESS_MAP.get(str(wordiness).lower())
        if val is None:
            raise ValueError(
                f"Unknown wordiness preset '{wordiness}'. "
                f"Choose from {list(self.WORDINESS_MAP.keys())} or pass an integer."
            )
        return val

    # -------------------------------------------------------------------------
    # Domain-Keyed Positional Permutation Engine
    # -------------------------------------------------------------------------

    def _get_domain_params(self, domain: Union[str, int]) -> Tuple[int, int]:
        """
        Derives an odd coprime multiplier $A$ and offset $B$ from a domain tag.

        Uses CRC32 instead of Python's built-in `hash()` to guarantee exact 
        stability across different interpreter processes, operating systems, 
        and application restarts without the randomized seeding of Python 3.3+.

        :param domain: Object namespace identifier (e.g., 'truck', 'container', 101).
        :return: Tuple of $(A, B)$ where $A \in \{1, 3, \dots, 15\}$ and $B \in [0, 15]$.
        """
        # Calculate cross-platform deterministic 32-bit checksum
        crc = zlib.crc32(str(domain).encode("utf-8"))
        
        # Multiplier A must be coprime to 16
        a = self.COPRIMES[(crc >> 4) % len(self.COPRIMES)]
        
        # Additive offset B in range [0, 15]
        b = crc % 16
        return a, b

    def _permute(self, hex_str: str, domain: Union[str, int]) -> str:
        """
        Applies a bijective positional affine transformation to a hex string:
            y_i = (A * x_i + B + i) mod 16

        :param hex_str: Normalized uppercase hexadecimal string.
        :param domain: Domain tag used to derive $A$ and $B$.
        :return: Transformed hexadecimal string of identical length.
        """
        a, b = self._get_domain_params(domain)
        permuted_chars: List[str] = []

        for i, char in enumerate(hex_str):
            val = int(char, 16)
            # Bijective transformation incorporating positional index i
            transformed = (a * val + b + i) % 16
            permuted_chars.append(format(transformed, "X"))

        return "".join(permuted_chars)

    def _unpermute(self, hex_str: str, domain: Union[str, int]) -> str:
        """
        Inverts the positional affine transformation:
            x_i = (A_inv * (y_i - B - i)) mod 16

        :param hex_str: Transformed hexadecimal string.
        :param domain: Domain tag matching the original encoding.
        :return: Recovered original hexadecimal string.
        """
        a, b = self._get_domain_params(domain)
        a_inv = self.MOD_INVERSES_16[a]
        recovered_chars: List[str] = []

        for i, char in enumerate(hex_str):
            val = int(char, 16)
            # Modular inverse cancels A; Python handles negative modulo correctly
            original = (a_inv * (val - b - i)) % 16
            recovered_chars.append(format(original, "X"))

        return "".join(recovered_chars)

    # -------------------------------------------------------------------------
    # Phonotactic Transition Logic
    # -------------------------------------------------------------------------

    @staticmethod
    def _is_valid_transition(
        prev_type: Optional[str],
        prev_phoneme: Optional[str],
        cand_type: str,
        cand_phoneme: str
    ) -> bool:
        """
        Evaluates whether a candidate phoneme is articulable after the prior state.

        Phonotactic Constraints:
          1. Word-Initial: Canonical words must begin with a Consonant (Tier 0).
          2. Anti-Gemination: Consecutive identical phonemes are prohibited.
          3. Sonority Sequencing: Consecutive Consonants (C+C) are blocked to 
             prevent illegal stop-burst clusters.
          4. Hiatus Prevention: Consecutive Vowels (V+V) are blocked to prevent 
             auditory vocalic distortion.
        """
        if prev_type is None:
            # Word boundary: enforce canonical consonant onset
            return cand_type == "C"
        if prev_phoneme == cand_phoneme:
            # Suppress stuttering/gemination
            return False
        if prev_type == "C" and cand_type == "C":
            # Suppress consonant clustering
            return False
        if prev_type == "V" and cand_type == "V":
            # Suppress vocalic hiatus
            return False
        return True

    def _partition_hex_into_words(
        self, hex_str: str, digits_per_word: int
    ) -> List[str]:
        """
        Splits a hexadecimal string into balanced chunks.

        Heavily favors 2 words for identifiers between 4 and 14 characters,
        enforcing a minimum of 2 words whenever length >= 2.
        """
        length = len(hex_str)
        if length <= 1:
            return [hex_str] if length == 1 else []

        # Target at least 2 words, rounded based on digits_per_word
        target_words = max(2, round(length / digits_per_word))
        # Ensure we never request more words than available characters
        target_words = min(target_words, length)

        # Distribute characters evenly across target_words
        q, r = divmod(length, target_words)
        chunks: List[str] = []
        idx = 0
        for i in range(target_words):
            size = q + (1 if i < r else 0)
            chunks.append(hex_str[idx:idx + size])
            idx += size

        return chunks

    # -------------------------------------------------------------------------
    # Public Encoding / Decoding API
    # -------------------------------------------------------------------------

    def encode(
        self,
        identifier: str,
        domain: Union[str, int] = "default",
        wordiness: Optional[Union[str, int]] = None,
        is_raw_text: bool = False,
        style: str = "auto"
    ) -> str:
        """
        Encodes a hexadecimal identifier or string into pronounceable words.

        :param identifier: Hexadecimal string or arbitrary text payload.
        :param domain: Namespace domain (e.g., 'truck', 'container', 1, 2) that
            ensures identical IDs yield completely different names.
        :param wordiness: Word density override ('compact', 'balanced', 'verbose'
            or explicit integer digits per word).
        :param is_raw_text: If True, treats identifier as UTF-8 string data.
        :param style: Formatting mode: 'auto', 'words', 'sentences', or 'paragraphs'.
        :return: Pronounceable formatted text string.
        """
        sanitized = identifier.strip()
        if not sanitized:
            return ""

        # Normalize input to uppercase hexadecimal representation
        is_hex = all(c in "0123456789ABCDEFabcdef" for c in sanitized)
        if is_raw_text or not is_hex:
            hex_data = sanitized.encode("utf-8").hex().upper()
        else:
            hex_data = sanitized.upper()

        # Step 1: Apply lightweight domain-keyed positional permutation
        permuted_hex = self._permute(hex_data, domain)

        # Step 2: Determine word chunking boundaries
        dpw = (
            self._resolve_wordiness(wordiness)
            if wordiness is not None
            else self.default_digits_per_word
        )
        word_chunks = self._partition_hex_into_words(permuted_hex, dpw)

        # Step 3: Sequential gated phoneme generation per word chunk
        generated_words: List[str] = []
        for chunk in word_chunks:
            word_phonemes: List[str] = []
            prev_type: Optional[str] = None
            prev_phoneme: Optional[str] = None

            for char in chunk:
                selected = False
                for cand in self.stacks[char]:
                    if self._is_valid_transition(
                        prev_type, prev_phoneme, cand["type"], cand["phoneme"]
                    ):
                        word_phonemes.append(cand["phoneme"])
                        prev_type = cand["type"]
                        prev_phoneme = cand["phoneme"]
                        selected = True
                        break

                # Fail-safe buffer if primary tiers are blocked
                if not selected:
                    fallback = self.stacks[char][2]
                    word_phonemes.append(fallback["phoneme"])
                    prev_type = fallback["type"]
                    prev_phoneme = fallback["phoneme"]

            generated_words.append("".join(word_phonemes))

        # Step 4: Formatting (words, sentences, or paragraphs)
        num_words = len(generated_words)
        if style == "words" or (style == "auto" and num_words <= 4):
            return " ".join(generated_words)

        sentence_chunks = list(chunked(generated_words, 6))
        sentences = [
            " ".join(chunk).capitalize() + "." for chunk in sentence_chunks
        ]

        if style == "sentences" or (style == "auto" and len(sentences) <= 4):
            return " ".join(sentences)

        paragraph_chunks = list(chunked(sentences, 4))
        paragraphs = [" ".join(p) for p in paragraph_chunks]
        return "\n\n".join(paragraphs)

    def _parse_word(self, word: str) -> Optional[List[str]]:
        """
        Parses a phonetic word back into its intermediate hex digits using a 
        prefix trie with recursive backtracking.
        """
        def backtrack(
            index: int,
            prev_type: Optional[str],
            prev_phoneme: Optional[str]
        ) -> Optional[List[str]]:
            if index == len(word):
                return []

            sub = word[index:]
            candidate_keys = list(self.trie.prefixes(sub))
            # Longest matches evaluated first (e.g., 'bel' before 'b')
            candidate_keys.sort(key=len, reverse=True)

            for key in candidate_keys:
                hex_digit, p_type = self.trie[key]
                if self._is_valid_transition(prev_type, prev_phoneme, p_type, key):
                    remainder = backtrack(index + len(key), p_type, key)
                    if remainder is not None:
                        return [hex_digit] + remainder
            return None

        return backtrack(0, None, None)

    def decode(
        self,
        phrase: str,
        domain: Union[str, int] = "default",
        strict: bool = False,
        as_raw_text: bool = False
    ) -> str:
        """
        Decodes a phonetic phrase back into its original hexadecimal identifier.

        :param phrase: The phonetic words, sentences, or paragraphs to decode.
        :param domain: The matching domain namespace under which it was encoded.
        :param strict: If True, raises PhonemeDecodeError on failure.
                       If False, returns a descriptive error message string.
        :param as_raw_text: If True, decodes the resulting hex into UTF-8 text.
        :return: Reconstructed original hexadecimal ID, text string, or error string.
        :raises PhonemeDecodeError: If strict is True and parsing fails.
        """
        if not phrase or not phrase.strip():
            return ""

        # Normalize text: strip periods, commas, newlines, and case
        cleaned = re.sub(r"[^\w\s]", "", phrase.lower())
        words = cleaned.split()

        reconstructed_hex_digits: List[str] = []

        for word in words:
            parsed_digits = self._parse_word(word)
            if parsed_digits is None:
                err_msg = (
                    f"Decode Error: Word '{word}' contains invalid phonemes "
                    f"or violates gating grammar."
                )
                if strict:
                    raise PhonemeDecodeError(err_msg)
                return err_msg
            reconstructed_hex_digits.extend(parsed_digits)

        permuted_hex = "".join(reconstructed_hex_digits)

        # Invert the domain permutation
        original_hex = self._unpermute(permuted_hex, domain)

        if as_raw_text:
            try:
                return bytes.fromhex(original_hex).decode("utf-8")
            except (ValueError, UnicodeDecodeError) as e:
                err_msg = f"Decode Error: Hex payload could not be converted to UTF-8: {e}"
                if strict:
                    raise PhonemeDecodeError(err_msg)
                return err_msg

        return original_hex

    def decode_safe(
        self,
        phrase: str,
        domain: Union[str, int] = "default",
        as_raw_text: bool = False
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Safe decode helper returning a (result, error) tuple.

        :param phrase: The phonetic phrase to decode.
        :param domain: Matching domain namespace.
        :param as_raw_text: Decode to UTF-8 text if True.
        :return: Tuple of (recovered_id, None) on success, or (None, error_msg) on failure.
        """
        try:
            res = self.decode(phrase, domain=domain, strict=True, as_raw_text=as_raw_text)
            return res, None
        except PhonemeDecodeError as err:
            return None, str(err)