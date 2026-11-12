**Headline:** A wide pair loaded from the database keeps its two star words, so editing or renaming it names its planets as before.

### Changed
- **`star_words` of a wide pair survive a save and load (DB.125).** `store.load_star_system` set `star_words` to nothing, so the next `assign_names` on a loaded wide pair drew two new words and renamed its stars and the planets of the secondary. The two words are not a column: `assign_names` names the stars `<first word> <primary word>` and `<first word> <secondary word>` (GEN.62), so the load reads them back from the two stored star names with the new `names.bodies.wide_pair_star_words` (which reuses `wide_pair_first_word`). Names of another shape (a pair renamed by hand) give nothing, as before. No migration, no column; stored data is unchanged.

generation-output: unchanged
generation-affects: none
