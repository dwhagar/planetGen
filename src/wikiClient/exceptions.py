# wikiClient/exceptions.py

"""Exception hierarchy shared by every `wikiClient` backend.

Every backend (`wikijs.WikiJsBackend`, `mediawiki.MediaWikiBackend`) raises
these same four types rather than its own -- calling code that talks to a
`client.WikiClient` never needs to know which wiki software is actually on
the other end to catch "bad credentials" vs. "page already exists" vs.
"everything else"."""


class WikiClientError(Exception):
    """Base class for every error this package raises."""


class WikiClientAuthError(WikiClientError):
    """The configured credentials were rejected -- an HTTP 401/403 response,
    a failed login, or an API-level error naming an authentication/
    authorization failure. Distinct from `WikiClientRequestError` so a
    caller can show "check your wiki credentials" specifically, rather than
    a generic failure."""


class WikiClientPageExistsError(WikiClientError):
    """`WikiClient.create_page` was called for a page that already exists
    on the target wiki. Every backend here is create-only -- none of them
    check for an existing page before creating, so this is how that
    conflict surfaces; the caller decides whether to pick a different
    path/title or show the conflict to whoever triggered the upload."""


class WikiClientRequestError(WikiClientError):
    """Everything else: the wiki instance couldn't be reached at all, it
    returned a non-2xx/unparseable response, or the create call itself
    reported failure for a reason other than a duplicate page."""
