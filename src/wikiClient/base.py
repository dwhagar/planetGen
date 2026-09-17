# wikiClient/base.py

"""Shared result type and backend interface every `wikiClient` backend
implements -- `client.WikiClient` is a thin dispatcher over whichever one
of these it was constructed with."""

from dataclasses import dataclass


@dataclass(frozen=True)
class WikiPage:
    """The result of a successful `create_page` call, from any backend.

    Attributes:
        id (int | str | None): The new page's id on the target wiki, if
            that wiki assigns one visible in the create response (a Wiki.js
            page id is always present; a MediaWiki `pageid` is too, but a
            backend that has nothing to report here uses `None` rather than
            inventing a value).
        path (str): The page's path/title on the target wiki, as it stored
            it (normalized -- may differ slightly from what was requested,
            e.g. a leading slash stripped, or MediaWiki's own title
            capitalization/whitespace normalization).
        title (str): The page's display title, as stored.
        url (str): A URL a caller can follow straight to the new page,
            without a second lookup.
    """

    id: object
    path: str
    title: str
    url: str


class WikiBackend:
    """Abstract interface a `wikiClient` backend implements.

    A concrete backend (`wikijs.WikiJsBackend`, `mediawiki.MediaWikiBackend`)
    takes whatever connection/credential arguments its own wiki software
    needs in `__init__`, but exposes this same `create_page` contract so
    `client.WikiClient` can call either one interchangeably.
    """

    def create_page(self, path, title, content, **kwargs):
        """
        Creates a new page on the target wiki. Every backend here is
        create-only (see `exceptions.WikiClientPageExistsError`).

        Args:
            path (str): The page's path/slug -- how a backend that
                separates "where the page lives" from "what it's called"
                (Wiki.js) addresses it. A backend whose wiki software has
                no separate concept of this (MediaWiki, where the title
                itself is the address) uses `path` as that address instead
                -- see that backend's own docstring.
            title (str): The page's display title.
            content (str): The page's body, in whatever markup the target
                wiki expects (backend-specific -- see each backend's own
                docstring for which of a caller's pre-rendered
                Markdown/wikitext copies it wants).
            **kwargs: Backend-specific optional arguments (e.g. Wiki.js's
                `tags`/`is_published`; MediaWiki's `summary`).

        Returns:
            WikiPage: The newly created page.

        Raises:
            WikiClientPageExistsError: A page already exists at this
                path/title.
            WikiClientAuthError: The configured credentials were rejected.
            WikiClientRequestError: The instance couldn't be reached,
                returned an unparseable response, or rejected the create
                for any other reason.
        """
        raise NotImplementedError
