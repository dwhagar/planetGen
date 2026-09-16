# wikijs/client.py

"""
`WikiJsClient`: a small, standalone GraphQL client for publishing pages to
a [Wiki.js](https://js.wiki/) instance.

Stdlib-only (`urllib.request`/`urllib.error`/`json`) -- this package has no
dependency on the rest of this project (no `stellarObjects`/`html`
imports, no `config.json`); `base_url` and `api_token` are passed straight
to the constructor by whatever caller wires this up later. Authenticates
with a Wiki.js Personal API Token (Admin -> API Access -> "New API Key"),
sent as `Authorization: Bearer <token>` on every request -- Wiki.js's
GraphQL API accepts this the same way it accepts a logged-in user's own
session JWT, with no separate login call needed.

Create-only: `create_page` never checks whether a page already exists at
the target path first: it always sends `pages.create` and lets Wiki.js's
own duplicate-path rejection surface as `WikiJsPageExistsError`. This
avoids a check-then-create race (a page created between an existence check
and the create call), and matches this package's own scope: there is
deliberately no `update_page`/`get_page`/`page_exists` here yet -- add
those as a separate increment if upsert behavior is ever wanted.
"""

import json
import urllib.error
import urllib.request
from dataclasses import dataclass

from .exceptions import WikiJsAuthError, WikiJsPageExistsError, WikiJsRequestError

_DEFAULT_TIMEOUT_SECONDS = 15

_CREATE_PAGE_MUTATION = """
mutation (
  $content: String!
  $description: String!
  $editor: String!
  $isPublished: Boolean!
  $isPrivate: Boolean!
  $locale: String!
  $path: String!
  $tags: [String]!
  $title: String!
) {
  pages {
    create(
      content: $content
      description: $description
      editor: $editor
      isPublished: $isPublished
      isPrivate: $isPrivate
      locale: $locale
      path: $path
      tags: $tags
      title: $title
    ) {
      responseResult {
        succeeded
        errorCode
        slug
        message
      }
      page {
        id
        path
        title
        updatedAt
      }
    }
  }
}
"""

# Substrings Wiki.js's own `pages.create` responseResult.message uses for a
# duplicate-path rejection (observed across Wiki.js 2.x releases -- there is
# no single documented, version-stable errorCode to key off instead, so this
# matches on the message text; confirm this still matches the exact wording
# of the target deployment's Wiki.js version, e.g. by triggering a real
# duplicate create against a dev instance, since a wording change there
# would otherwise silently fall through to WikiJsRequestError instead).
_DUPLICATE_PATH_MESSAGE_MARKERS = ("already exists", "already in use", "duplicate")

# Substrings indicating an authentication/authorization failure inside a
# GraphQL response's top-level "errors" array (as opposed to an HTTP-level
# 401/403, which is checked separately) -- Wiki.js returns these for an
# invalid/expired/insufficiently-privileged API token.
_AUTH_ERROR_MESSAGE_MARKERS = ("unauthorized", "access denied", "not authorized", "invalid authentication")


@dataclass(frozen=True)
class WikiPage:
    """The result of a successful `create_page` call.

    Attributes:
        id (int): The new page's Wiki.js internal id.
        path (str): The page's path, as Wiki.js stored it (normalized --
            may differ slightly from the requested `path`, e.g. a leading
            slash stripped).
        title (str): The page's title, as stored.
        url (str): `base_url` + `path` -- not returned by Wiki.js itself,
            built here so a caller can link straight to the new page
            without a second lookup.
    """

    id: int
    path: str
    title: str
    url: str


class WikiJsClient:
    """A client for one Wiki.js instance's GraphQL API.

    Args:
        base_url (str): The instance's root URL, e.g.
            `"https://wiki.example.com"` -- a trailing slash is stripped;
            `/graphql` is appended for every request.
        api_token (str): A Wiki.js Personal API Token (Admin -> API Access),
            sent as `Authorization: Bearer <api_token>`.
        timeout (float): Per-request socket timeout, in seconds.
    """

    def __init__(self, base_url, api_token, timeout=_DEFAULT_TIMEOUT_SECONDS):
        self._base_url = base_url.rstrip("/")
        self._api_token = api_token
        self._timeout = timeout

    def create_page(
        self,
        path,
        title,
        content,
        description="",
        editor="markdown",
        locale="en",
        tags=None,
        is_published=True,
        is_private=False,
    ):
        """
        Creates a new page on the target Wiki.js instance.

        Args:
            path (str): The page's path (no leading slash, e.g.
                `"systems/kepler-442"`) -- must not already exist on this
                instance (see the module docstring's "Create-only" note).
            title (str): The page's title.
            content (str): The page's body, in `editor`'s format (Markdown
                by default).
            description (str): Short summary shown in Wiki.js's page
                listings/search results.
            editor (str): Wiki.js editor type -- `"markdown"` matches this
                project's own generated `markdown_content`.
            locale (str): Wiki.js locale code for the page.
            tags (Sequence[str], optional): Page tags.
            is_published (bool): Whether the page is immediately visible
                (vs. saved as a draft).
            is_private (bool): Whether the page is restricted to specific
                users/groups (Wiki.js's own access rules decide who, this
                only flags the page as private).

        Returns:
            WikiPage: The newly created page.

        Raises:
            WikiJsPageExistsError: A page already exists at `path`.
            WikiJsAuthError: `api_token` was rejected.
            WikiJsRequestError: The instance couldn't be reached, returned
                an unparseable response, or rejected the create for any
                other reason.
        """
        variables = {
            "content": content,
            "description": description,
            "editor": editor,
            "isPublished": bool(is_published),
            "isPrivate": bool(is_private),
            "locale": locale,
            "path": path,
            "tags": list(tags) if tags else [],
            "title": title,
        }
        data = self._graphql(_CREATE_PAGE_MUTATION, variables)

        try:
            result = data["pages"]["create"]
            response_result = result["responseResult"]
            succeeded = response_result["succeeded"]
        except (KeyError, TypeError) as exc:
            raise WikiJsRequestError(f"Unexpected response shape from Wiki.js pages.create: {exc}")

        if not succeeded:
            message = response_result.get("message") or "Wiki.js rejected the page create request."
            lowered = message.lower()
            if any(marker in lowered for marker in _DUPLICATE_PATH_MESSAGE_MARKERS):
                raise WikiJsPageExistsError(message)
            raise WikiJsRequestError(message)

        page = result.get("page") or {}
        try:
            page_path = page["path"]
        except KeyError:
            raise WikiJsRequestError("Wiki.js reported success but returned no page in its response.")

        return WikiPage(
            id=page["id"],
            path=page_path,
            title=page["title"],
            url=f"{self._base_url}/{page_path.lstrip('/')}",
        )

    def _graphql(self, query, variables):
        """
        Runs one GraphQL request and returns its `data` object.

        Raises:
            WikiJsAuthError: An HTTP 401/403, or a top-level GraphQL
                `errors` entry naming an auth failure.
            WikiJsRequestError: Any other non-2xx/unreachable/unparseable
                response, or a non-auth top-level GraphQL error.
        """
        url = f"{self._base_url}/graphql"
        body = json.dumps({"query": query, "variables": variables}).encode("utf-8")
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json",
            "Authorization": f"Bearer {self._api_token}",
        }
        request = urllib.request.Request(url, data=body, headers=headers, method="POST")

        try:
            with urllib.request.urlopen(request, timeout=self._timeout) as response:
                raw_body = response.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            detail = _read_error_detail(exc)
            if exc.code in (401, 403):
                raise WikiJsAuthError(f"Wiki.js rejected the API token ({exc.code}): {detail}")
            raise WikiJsRequestError(f"Wiki.js returned HTTP {exc.code}: {detail}")
        except urllib.error.URLError as exc:
            raise WikiJsRequestError(f"Could not reach Wiki.js at {self._base_url}: {exc.reason}")

        try:
            parsed = json.loads(raw_body)
        except ValueError as exc:
            raise WikiJsRequestError(f"Wiki.js returned an unparseable response: {exc}")

        errors = parsed.get("errors")
        if errors:
            combined = "; ".join(error.get("message", str(error)) for error in errors)
            lowered = combined.lower()
            if any(marker in lowered for marker in _AUTH_ERROR_MESSAGE_MARKERS):
                raise WikiJsAuthError(f"Wiki.js rejected the request: {combined}")
            raise WikiJsRequestError(f"Wiki.js GraphQL error: {combined}")

        return parsed.get("data") or {}


def _read_error_detail(http_error):
    try:
        return http_error.read().decode("utf-8", errors="replace") or http_error.reason
    except Exception:
        return http_error.reason
