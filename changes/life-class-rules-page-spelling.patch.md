**Headline:** The life class rules page is rebuilt with the American spelling DOC.103 gave the rule file, and its test passes again.

### Changed

- `docs/design/habitability/life-class-rules.md` is rebuilt from the rule file (`scripts/build_life_class_rules.py`): one word, "center" for "centre", had changed in the file and not in the page, so `test_life_class_rules_doc::test_the_page_is_the_rule_file_as_built` failed on main.

generation-output: unchanged
