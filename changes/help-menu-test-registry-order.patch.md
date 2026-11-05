**Headline:** The Help menu test no longer pins the exact list of help pages, so a new help page such as the sky page does not turn it red.

### Changed

- `test_the_help_menu_holds_classes_and_habitability_and_the_old_links_are_gone` compared the whole menu with a fixed list. The menu is now built from the help registry (`web/lib/help_registry.py`), and VIEW.38 registered "The sky behind a system" in it, so the test failed on main. The test now demands All help first, References last, and Classes, Habitability and Life in that order, and allows other registered pages between them. Tests only.

generation-output: unchanged
