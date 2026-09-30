### Removed
- **The CGI pages are gone.** Every page is served by the Flask app, so the
  `src/html/*.py` redirect shims, the CGI page shell (`lib/page.py`),
  `fmt.post_link`/`data_nav_params`, `static/navform.js`, the POST mode of
  the pager and the helpers' `LEGACY_PAGES` are removed, along with the
  CSS for the old side rail and link buttons.
- **The example Apache vhost has no CGI rules.** `ScriptAliasMatch`,
  `ExecCGI`, `AddHandler cgi-script .py`, the `<Files "wsgi.py">` override
  and the CGI-only `SetEnv` lines are gone; `Alias /static/` and
  `WSGIScriptAlias /` remain. `install.sh` enables `wsgi` instead of
  `cgid`, and prints what to remove when an existing site file still has
  the CGI rules (see `docs/apache-deployment.md`, "Updating an existing
  server").

### Changed
- **Old `/<name>.py` URLs redirect in the app.** `web/old_urls.py` answers
  `/index.py`, `/sector.py?id=5`, `/search.py?...` and the rest with a 301
  to the page that replaced them, keeping their parameters, so old
  bookmarks still work once the vhost's `ScriptAliasMatch` is removed. Any
  other `.py` name is a 404.
