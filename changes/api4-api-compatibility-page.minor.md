**Headline:** The docs now say which API version each release speaks, what each version changed, and which clients work with it.

### Added
- **API compatibility data (API.4).** `docs/reference/api-versions.json` holds one row per API version (the release that first carried it, the routes and payload formats it changed, which clients speak it); `scripts/build_api_compatibility.py` builds `docs/reference/api-compatibility.md` from it (`--check` fails when the page is stale). `docs/reference/api.md` links the page.
- **Release tooling keeps the history honest.** A change that raises `API_VERSION` must add its row to the history file and carry an `api-version: changed` line in its release note (`bump_version.py --check-pr` fails if only some go together); the stamp fills the row's first release.
