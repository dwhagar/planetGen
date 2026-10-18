# planetgen/util/user_agent.py

"""
The `User-Agent` of every request the software makes to another site (API.31).

One function builds it, so this is the only place that spells the string: the release from `planetgen/_version.py`,
the project URL in a comment and the API version from `planetgen/api/version.py`, in the RFC 9110 form of product
tokens with a comment, for example `PlanetGen/8.0.1721 (+https://github.com/dwhagar/planetGen) PlanetGen-API/4`.
"""

from planetgen._version import REPO_URL, __version__
from planetgen.api.version import API_VERSION

PRODUCT = "PlanetGen"
"""str: The software's product token."""


def user_agent():
    """The `User-Agent` for an outgoing request: release, project URL and API version."""
    return f"{PRODUCT}/{__version__} (+{REPO_URL}) {PRODUCT}-API/{API_VERSION}"
