# planetgen/web/preferences_page.py

"""
The Preferences page (UX.236): `/preferences` lists each per-browser setting with its plain-words description and its default, and a
Reset to defaults button. The page is drawn in the browser by `static/preferencespage.js` from `static/mappreferences.js`, the one place
that knows the keys, their ranges and their words, and the settings are kept in the browser's `localStorage` (the same keys move into
the account when accounts land, USR.7). The server holds no setting, so the page is the same for everyone.
"""

from . import bp
from .helpers import crumb, render_page


@bp.route("/preferences")
def preferences():
    """The settings a visitor can change for this browser, starting with the Galaxy Map's."""
    return render_page(
        "preferences.html",
        title="Preferences",
        section="preferences",
        breadcrumbs=[crumb("Preferences")],
        description="The settings this browser keeps for the Galaxy Map, each with what it does and its default.",
    )
