**Headline:** Every kind of object the generator makes is now one row in a single registry, and the lists that used to repeat the kinds read it.

### Changed
- A new phenomena registry (`planetgen/generation/phenomena/registry.py`) holds one row per kind: family, rate, gate, mode, builder, system model, what it carries, object ID code, class page and map marker. The sector roll list, the scatter list, the `--type` choices, the object ID codes, the naming key kinds and the rate scale table now read from it; nothing the generator makes changes.
- A completeness test checks every row has every column, every builder resolves, no ID code repeats, and that the database tables, insert functions, class pages and the scatter table's constraint name exactly the registry's kinds.

TODO items resolved: GEN.250 / Added: none
