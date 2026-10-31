# planetgen/util/wordlist.py

"""
A list of things written as words (DOC.86, rule 15 of `docs/design/interface/text-style.md`)
===========================================================================================

`join_words(["red", "green", "blue"])` is "red, green, and blue": every list of three or more puts a comma before its "and" or "or"
(Boss, 2:50 PM EDT 2026-10-31), and a pair has none ("red and blue"). Text built from a list of names goes through here, so a sentence
the code writes keeps the rule the text files keep. Standard library only.
"""


def join_words(items, word="and"):
    """
    `items` as words: "a", "a and b", "a, b, and c" (or with `word="or"`). An empty list is an empty string.
    """
    items = [str(item) for item in items]
    if len(items) <= 2:
        return f" {word} ".join(items)
    return f"{', '.join(items[:-1])}, {word} {items[-1]}"
