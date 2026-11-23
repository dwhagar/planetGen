# tests/trojan_support.py

"""The check the order and spacing tests of a star's bodies share for a trojan (TEST.818, TEST.821)."""


def assert_trojans_on_their_giants(planets):
    """The trojans of `planets` (a region or a planet with a `lagrange_point`, GEN.1150, GEN.1151), each checked against its giant, and the list of them. A trojan shares its giant's orbit by design, so it is no slot of its own: it stands at its
    giant's distance, 60 degrees ahead of or behind it (L4 or L5), and its giant is a gas giant of the same list. A test of the order or the spacing of the bodies takes the returned trojans out of its list and keeps every other check as it is."""
    trojans = [obj for obj in planets if getattr(obj, "lagrange_point", None) is not None]
    for trojan in trojans:
        host = getattr(trojan, "lagrange_host", None) or getattr(trojan, "host", None)
        assert trojan.lagrange_point in ("L4", "L5"), f"trojan at {trojan.distance} AU stands at {trojan.lagrange_point!r}"
        assert host is not None and any(host is obj for obj in planets) and host.body_type == "g", f"trojan at {trojan.distance} AU has no giant of its star's list for its host"
        assert trojan.distance == host.distance, f"trojan at {trojan.distance} AU is not on its giant's orbit at {host.distance} AU"
    return trojans
