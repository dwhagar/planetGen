# The plain-words style guide (DOC.83)

Informs: DOC.83, UX.252 (planet, moon and belt text), UX.253 (star, system and star-class text), DOC.84 (the life help pages), GEN.1208 (class descriptions). Related: API.33 (units), DOC.15 (glossary), DOC.57 and DOC.58 (the "?" help links).
Status: written. The checkable rules are `planetgen/util/plainwords.py`; `src/tests/test_plain_words.py` reads the generated text against them.

## What it is for

A visitor reads a planet's page to learn what it is like, not to decode it. Every sentence the site writes about a planet, a moon, a belt, a star or a system is meant to be read once, at a glance, by someone who has never seen the generator. Today the sentences come from many builders (`describe/build.py`, `describe/sheet.py`, `web/lib/classref.py`, `web/lib/planet_layout.py`, the templates) and each was written on its own day, so the same fact is worded three ways and the same class is named two. This note fixes the rules once. The two rewrite items apply it; a test keeps new text from breaking it.

## The rules

### Checked by the test

1. **Short sentences.** At most 30 words (`MAX_WORDS`); aim at 20. A sentence over 30 words is two ideas. Cut it at the colon, the semicolon or the "which".
2. **One idea each.** At most one semicolon (`MAX_SEMICOLONS`). A second semicolon is a second idea: give it its own sentence. A list inside a sentence ("1: 10 to 100 AU; 2: ...") is a table, not a sentence.
3. **Plain words.** The words and phrases in `BANNED` are not used; each has the plain word beside it (`approximately` is "about", `in order to` is "to", `e.g.` is "for example"). The list also keeps out the words of the program's own workshop (`stand-in`, `tweak`, `epoch`, `placeholder`): a reader does not need to know how a number was set.
4. **One name for a class.** A life class is named by its whole id: "Class 3a", "Class 2", never "Class 3" for a class that has letters (the classes are `physics/life_classes.py`'s ids). When a sentence means the whole family, say "the Class 3 family" or name the classes: "Class 3a and 3b".
5. **Labels in sentence case.** A table label starts with a capital and the rest is lower case except names and abbreviations: "Surface temperature", "Distance from the star", "Urey ratio". Not "Surface Temperature". Where a label needs a possessive, write it plain ("Planetary limit of the primary").
6. **American spelling.** The site is in American English (Boss, 10:33 AM EDT 2026-10-31): "color", "gray", "center", "meter", "kilometer", "liter", "neighbor", "vapor", "sulfur", "analyze", "organize", "modeling". Not mixed. `BRITISH` and `ISE_STEMS` in `plainwords.py` list the spellings the checker flags and the American ones; text in quotation marks and the names in `BRITISH_NAMES` (an organization's own name) are left alone.

### Kept by the writer (not machine-checkable)

7. **Say the answer first.** The first sentence of a paragraph says the thing; the reason or the exception follows.
8. **Never repeat a fact on a page.** A fact stated in the summary is not stated again in the description, the notes or a table caption. A table row and a sentence may both give a number only when the sentence is the only place that explains it.
9. **One unit for one quantity, spelled as the page's unit setting says.** A quantity is a typed quantity in the Document (API.33; `tests/test_counts_have_units.py` checks counts), so the unit follows the visitor's setting and no sentence hard-codes `AU` or `km`. A number in a sentence always has its noun or unit beside it.
10. **Active, concrete verbs.** "The star warms the planet to 280 K", not "The planet is characterised by a surface temperature of 280 K". No stacked "of the ... of the ...".
11. **No stage or age in a class description.** A class description says what the class is, never how old its life is or which stage it has reached; that is the world's timeline (GEN.1208).
12. **No jargon without a help link.** A term a visitor may not know ("Roche limit", "hycean", "suitability") is linked to its help page (DOC.57's "?" link, the DOC.15 glossary) the first time it appears on a page, and its explanation lives on the help page, not inline.
13. **No process words.** Text never mentions the generator, its version, a seed, a stand-in or a tier. The Debug view and the Properties sheet are where those live.
14. **No codes.** No item code (`GEN.866`) reaches the page; `tests/test_no_todo_codes_shown.py` checks this.
15. **Oxford comma (the serial comma).** Every list of three or more puts a comma before its "and" or "or": "red, green, and blue", never "red, green and blue" (Boss, 2:50 PM EDT 2026-10-31: always, in every text page, help page and generated sentence). Kept by the writer until DOC.86 adds the check.
16. **No out-of-universe words outside References.** The words in the banned list (stand-in, tier, fitted, speculative, generator and the like, set by DOC.87) appear only on the References pages (`/references`), which carry the Out of universe banner; everywhere else the site speaks in universe and links to the References page instead. Checked by the page crawl once DOC.87 and DOC.93 land.

## How the test uses the checker

`tests/test_plain_words.py` builds the Document of ten seeded systems and reads every prose sentence and table label of it, then the class pages (`web/lib/classref.py`) and the label of every row of the property sheet (`describe/sheet.py`). Each break of rules 1 to 6 must be in `KNOWN`, which maps the rule and the start of the sentence (or the label) to the item that owns the fix. A break that is not listed fails the test, so new text keeps the guide while old text is rewritten. `KNOWN` is the list of sentences and labels the two rewrite items must fix: UX.252 for planets, moons, belts, comets and their class pages, UX.253 for stars, systems and star-class pages. Each item removes its entries as it rewrites, and is done when none carry its code.

The life help pages (`web/lib/life_help.py`, DOC.84) are held to the same rules by `tests/test_life_help.py`, with no `KNOWN` list: they were written to the guide. The rules themselves are tested on made-up sentences (`test_the_checker_flags_each_rule`), so a rule cannot go quiet without a test failing.

## Adding a rule

A rule goes into `plainwords.py` with a test sentence for it in the checker's test, then a pass over the corpus to see what it flags. A rule that flags current text adds those entries to `KNOWN` under the item that will fix them; it never loosens the rule to make the current text pass.
