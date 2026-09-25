### Changed
- **System and phenomenon pages moved to the new site.** A star system is
  now at `/system/<id>`, the phenomena list at `/phenomena` (paged with
  `?page=N`), and one phenomenon at `/phenomenon/<type>/<id>`. These are
  plain, bookmarkable URLs with no database name in them, and every link
  on them (neighbouring systems, sector, Navigate from/to here, the
  Wikitext/Markdown views with `?code=...`) is an ordinary link, so Back,
  reload and open-in-new-tab work. The pages use the new header and
  breadcrumbs, and the system map, body list, code views with the Copy
  button and phenomenon diagram work as before. Old `system.py`,
  `phenomena.py` and `phenomenon.py` links redirect to the new addresses.
- The admin "Upload to Wiki" form on a system page is now CSRF-protected
  and redirects back to the page afterwards with a fixed status message,
  so reloading never uploads twice. It also says so when the default
  admin credentials must be changed first.
