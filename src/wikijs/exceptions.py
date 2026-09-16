# wikijs/exceptions.py

"""Exception hierarchy for the `wikijs` package."""


class WikiJsError(Exception):
    """Base class for every error this package raises."""


class WikiJsAuthError(WikiJsError):
    """The configured API token was rejected -- an HTTP 401/403 response,
    or a GraphQL-level error naming an authentication/authorization
    failure. Distinct from `WikiJsRequestError` so a caller can show "check
    your Wiki.js API token" specifically, rather than a generic failure."""


class WikiJsPageExistsError(WikiJsError):
    """`WikiJsClient.create_page` was called with a `path` that already has
    a page on the target Wiki.js instance. This client is create-only (see
    `client.py`'s module docstring) -- it never checks for an existing page
    before creating, so this is how that conflict surfaces; the caller
    decides whether to pick a different path, or show the conflict to
    whoever triggered the upload."""


class WikiJsRequestError(WikiJsError):
    """Everything else: the Wiki.js instance couldn't be reached at all, it
    returned a non-2xx/unparseable response, or the `pages.create` mutation
    itself reported failure (`responseResult.succeeded == False`) for a
    reason other than a duplicate path."""
