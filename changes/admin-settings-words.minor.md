**Headline:** Every setting now has a plain-language name, a one-line description and a purpose-based group, ready for the redesigned Settings page.

### Changed
- **Settings are named in words.** Each of the 89 options in `config.json` carries a label, a one-sentence summary, a group, an optional control kind and an Advanced mark in the settings model (`settings.opt`), instead of a title built from its key. The Settings page shows the label as each setting's heading, `docs/config.md` lists it beside the key, and a test holds every option to the same wording rules (short label without key punctuation, one sentence, a known group).
- **Groups follow purpose, not file layout.** The page's sections are now Site identity, Look, Search engines, Link previews, Speed and disk, Abuse protection and Housekeeping for what an admin can change, then the server-only groups (database connection, logging, reverse proxy, security, folders and addresses, wiki upload, settings file).
