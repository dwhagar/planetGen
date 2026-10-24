**Headline:** The wiki upload is gone: planetGen no longer talks to Wiki.js or MediaWiki, and a system's write-up is exported as Markdown (or plain text).

### Removed
- **The wiki upload and download code (OPS.69).** The `planetgen.wiki` package (the MediaWiki and Wiki.js clients), the wikitext renderer, `POST /api/systems/<id>/wiki`, `POST /api/sectors/<id>/wiki`, `GET /api/wiki-config`, the "Upload to Wiki" forms on the system and sector pages, the sector "Set wiki link" form, the `wiki` settings section (`PLANETGEN_WIKIJS_*`, `PLANETGEN_MEDIAWIKI_*`), and the wiki integration tests. The old addresses answer 404.
- **The wikitext format.** `GET /api/systems/<id>/text`, `planetgen system --format` and the one-off system page offer Markdown and plain text only; the `?code=wikitext` view and the `.wiki` download are gone.
- **The wiki link columns (schema v149, migration 0149).** `sectors.wiki_url`, `star_systems.wikijs_url` and `star_systems.mediawiki_url` are dropped.

### Changed
- `GET /api/systems/<id>/text` answers Markdown when no `format` is given (it was wikitext). The Markdown export, the Markdown code view with its Copy button, the `.md` download and `planetgen system --format markdown --output FILE` work as before.
