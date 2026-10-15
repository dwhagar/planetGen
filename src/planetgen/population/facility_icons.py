# planetgen/population/facility_icons.py

"""
The icons a facility type may use (POP.10).

A type's `icon` is the name of one Bootstrap Icons glyph (MIT licence, the
set Shoelace ships). The pages are served under a `default-src 'self'`
policy, so each icon is a vendored file; this is the list of names the admin
may pick from and the only names a type may store.
"""

ICON_NAMES = (
    "activity", "award", "bank", "book", "box-seam", "bricks", "broadcast", "broadcast-pin", "building",
    "buildings", "bullseye", "cart", "cash-stack", "check-circle-fill", "cloud", "compass", "cpu", "crosshair",
    "dash-circle-fill", "diagram-3", "droplet", "exclamation-triangle-fill", "flag-fill", "flask", "flower1",
    "fuel-pump", "gear", "gem", "geo-alt", "globe", "globe2", "hammer", "hdd-network", "heart-pulse", "hexagon",
    "hospital", "house-door", "leaf", "lightbulb", "lightning-charge", "moon", "moon-stars", "mortarboard",
    "people-fill", "pin-map", "radioactive", "rocket", "rocket-takeoff", "shield", "shield-check",
    "shield-exclamation", "shop", "star", "star-fill", "stars", "sun", "thermometer-half", "tools", "tree",
    "trophy", "truck", "wifi", "wrench-adjustable", "x-octagon-fill",
)
"""tuple: Every icon name a facility type may store."""
