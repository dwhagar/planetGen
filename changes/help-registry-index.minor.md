**Headline:** The Help menu now has an index page that lists every help page, the class types, and the life topics and classes, all read from the program's own tables.

### Changed
- **DOC.6.** A new help registry (`web/lib/help_registry.py`) holds no list of its own: each part of the program registers its own help page, category and link terms (classes, life, properties and habitability do so now). The Help menu, the new `/help` index and the sitemap read what is registered. The pages under a topic (the class types, the life topics and life classes) come from `classref.catalog()` and `life_classes.CLASSES`, so a new class appears with no edit.
- The habitability page no longer carries typed-in numbers: the gear limits, the score edges, the life class count on the map row and the temperature bands with their shares of the surface are read from `protection_limits`, `habitability`, `life_classes` and `regional_temperature`. The text the page shows is unchanged.
- `tests/test_help_registry.py` fails when a `/help` route is missing from the registry, the menu or the index differs from it, or a number on the habitability page drifts from the code. Display only, no stored value changes.
generation-output: unchanged
