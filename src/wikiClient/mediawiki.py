# wikiClient/mediawiki.py

"""
`MediaWikiBackend`: a small client for publishing pages to a
[MediaWiki](https://www.mediawiki.org/) instance's Action API (`api.php`)
-- the other `wikiClient` backend (see `wikijs.py` for the Wiki.js one),
dispatched to by `client.WikiClient`.

Stdlib-only (`urllib.request`/`urllib.error`/`http.cookiejar`/`json`).
Authenticates with a [Bot Password](https://www.mediawiki.org/wiki/Special:BotPasswords)
(`username` like `"SomeBot@planetgen"`, plus the password `Special:BotPasswords`
generated for it) -- MediaWiki's standard machine-account login, needing no
OAuth consumer registration. Unlike `WikiJsBackend`'s bearer token, this API
is session-based: login and every later request share cookies via an
internal `http.cookiejar.CookieJar`, so one `MediaWikiBackend` instance logs
in at most once (lazily, on its first `create_page` call) and reuses that
session for every call after.

Create-only: `create_page` always sends `action=edit` with `createonly=1`,
which MediaWiki itself rejects (rather than overwriting) if the page
already exists -- surfaced here as `WikiClientPageExistsError`, the same
check-then-create race this avoids as `WikiJsBackend`'s own "Create-only"
note explains.

MediaWiki has no separate "path" distinct from a page's title the way
Wiki.js does -- the title itself is the address (subpages nest with `/`
the same way a Wiki.js `path` does). So `create_page`'s `path` argument is
used as the actual MediaWiki page title here; `title` is accepted only for
interface parity with `WikiBackend.create_page` and is not sent to
MediaWiki at all (there is no way to give a MediaWiki page an address
distinct from its title, short of the `{{DISPLAYTITLE}}` magic word, which
is out of scope here).
"""

import http.cookiejar
import json
import urllib.error
import urllib.parse
import urllib.request

from .base import WikiBackend, WikiPage
from .exceptions import WikiClientAuthError, WikiClientPageExistsError, WikiClientRequestError

_DEFAULT_TIMEOUT_SECONDS = 15

# `action=edit`'s error `code` for "a page already exists at this title and
# createonly=1 was set" -- the one, stable, documented value MediaWiki uses
# for this (unlike Wiki.js, which has no errorCode stable enough to key off
# instead of message text -- see wikijs.py's own note).
_PAGE_EXISTS_ERROR_CODE = "articleexists"

# `action=edit`/`action=query` error `code`s that mean the session's
# credentials aren't good enough for this request -- an expired/invalid
# login, a blocked account, or a permission the bot password's grants don't
# include.
_AUTH_ERROR_CODES = ("permissiondenied", "readapidenied", "mustbeloggedin", "blocked")


class MediaWikiBackend(WikiBackend):
    """A client for one MediaWiki instance's Action API.

    Args:
        base_url (str): The instance's API entry point's directory, e.g.
            `"https://wiki.example.com/w"` for an instance whose API lives
            at `https://wiki.example.com/w/api.php` -- a trailing slash is
            stripped; `/api.php` is appended for every request. Pass the
            full path up to (not including) `api.php` if it isn't at the
            instance root.
        username (str): A Bot Password username, in `"User@BotName"` form
            (`Special:BotPasswords` on the target instance).
        password (str): That bot password.
        timeout (float): Per-request socket timeout, in seconds.
    """

    def __init__(self, base_url, username, password, timeout=_DEFAULT_TIMEOUT_SECONDS):
        self._api_url = f"{base_url.rstrip('/')}/api.php"
        self._username = username
        self._password = password
        self._timeout = timeout
        self._cookie_jar = http.cookiejar.CookieJar()
        self._opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(self._cookie_jar))
        self._logged_in = False

    def create_page(self, path, title, content, summary="", **kwargs):
        """
        Creates a new page on the target MediaWiki instance.

        Args:
            path (str): The page's title (see the module docstring --
                MediaWiki addresses a page by its title alone, so this is
                sent as `action=edit`'s `title`).
            title (str): Accepted for interface parity with
                `WikiBackend.create_page`; not sent to MediaWiki (see the
                module docstring).
            content (str): The page's body, in wikitext -- pass a caller's
                `wikitext_content`, not `markdown_content` (see `wikijs.py`
                for that one; MediaWiki has no Markdown editor of its own
                to hand Markdown to the way Wiki.js does).
            summary (str): Edit summary shown in the page's history.

        Returns:
            WikiPage: The newly created page.

        Raises:
            WikiClientPageExistsError: A page already exists at this title.
            WikiClientAuthError: `username`/`password` were rejected, or
                the account can't create pages.
            WikiClientRequestError: The instance couldn't be reached,
                returned an unparseable response, or rejected the create
                for any other reason.
        """
        self._ensure_logged_in()
        csrf_token = self._fetch_token("csrf")

        response = self._call(
            {
                "action": "edit",
                "title": path,
                "text": content,
                "summary": summary,
                "createonly": "1",
                "token": csrf_token,
                "format": "json",
            }
        )

        error = response.get("error")
        if error is not None:
            self._raise_for_error(error)

        try:
            edit = response["edit"]
            if edit.get("result") != "Success":
                raise KeyError("result")
            page_id = edit["pageid"]
            page_title = edit["title"]
        except KeyError as exc:
            raise WikiClientRequestError(f"Unexpected response shape from MediaWiki action=edit: {exc}")

        base_url = self._api_url.rsplit("/api.php", 1)[0]
        return WikiPage(
            id=page_id,
            path=page_title,
            title=page_title,
            url=f"{base_url}/index.php?title={urllib.parse.quote(page_title.replace(' ', '_'))}",
        )

    def _ensure_logged_in(self):
        """Logs in at most once per instance, on the first call that needs
        it -- every later call reuses the same cookie-jar session."""
        if self._logged_in:
            return

        login_token = self._fetch_token("login")
        response = self._call(
            {
                "action": "login",
                "lgname": self._username,
                "lgpassword": self._password,
                "lgtoken": login_token,
                "format": "json",
            }
        )

        login = response.get("login") or {}
        if login.get("result") != "Success":
            reason = login.get("reason") or login.get("result") or "MediaWiki rejected the login."
            raise WikiClientAuthError(f"MediaWiki login failed: {reason}")

        self._logged_in = True

    def _fetch_token(self, token_type):
        response = self._call({"action": "query", "meta": "tokens", "type": token_type, "format": "json"})
        try:
            return response["query"]["tokens"][f"{token_type}token"]
        except KeyError as exc:
            raise WikiClientRequestError(f"MediaWiki returned no {token_type} token: {exc}")

    def _raise_for_error(self, error):
        code = error.get("code", "")
        info = error.get("info") or "MediaWiki rejected the request."
        if code == _PAGE_EXISTS_ERROR_CODE:
            raise WikiClientPageExistsError(info)
        if code in _AUTH_ERROR_CODES:
            raise WikiClientAuthError(info)
        raise WikiClientRequestError(f"MediaWiki API error ({code}): {info}")

    def _call(self, params):
        """
        Runs one `api.php` POST and returns its parsed JSON body. Cookies
        set by any previous call on this instance (the login session) are
        sent automatically, via `self._opener`'s cookie jar.

        Raises:
            WikiClientAuthError: An HTTP 401/403.
            WikiClientRequestError: Any other non-2xx/unreachable/
                unparseable response.
        """
        body = urllib.parse.urlencode(params).encode("utf-8")
        request = urllib.request.Request(
            self._api_url,
            data=body,
            headers={"Content-Type": "application/x-www-form-urlencoded", "Accept": "application/json"},
            method="POST",
        )

        try:
            with self._opener.open(request, timeout=self._timeout) as response:
                raw_body = response.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            detail = _read_error_detail(exc)
            if exc.code in (401, 403):
                raise WikiClientAuthError(f"MediaWiki rejected the request ({exc.code}): {detail}")
            raise WikiClientRequestError(f"MediaWiki returned HTTP {exc.code}: {detail}")
        except urllib.error.URLError as exc:
            raise WikiClientRequestError(f"Could not reach MediaWiki at {self._api_url}: {exc.reason}")

        try:
            return json.loads(raw_body)
        except ValueError as exc:
            raise WikiClientRequestError(f"MediaWiki returned an unparseable response: {exc}")


def _read_error_detail(http_error):
    try:
        return http_error.read().decode("utf-8", errors="replace") or http_error.reason
    except Exception:
        return http_error.reason
