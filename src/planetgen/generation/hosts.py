"""One star (or pair of stars) that holds planets, and the planets and belts that go round it (GEN.1664).

`StarSystem.hosts` is the list of these in the one order every pass over a system's bodies keeps: the system's own star (a
single star, or a close pair's merged proxy), a wide pair's secondary, each further star of the tree, then the proxy of each
inner pair node that holds circumbinary planets. Replaces the three lists `planets`, `secondary_planets` and `member_planets`
(and `node_planets`) the system kept before, so a pass over the bodies is written once.
"""

PRIMARY = "primary"      # the system's star, or a close pair's proxy (the `stars` row role `primary`)
SECONDARY = "secondary"  # the wide pair's second star
MEMBER = "member"        # a star beyond the pair (`StarSystem.member_stars`)
NODE = "node"            # the proxy of an inner pair node of the tree (`StarSystem.node_hosts`)


class Host:
    """The star `star` (or the proxy of a pair), its `role`, and `orbiters`, the planets and belts that go round it."""

    __slots__ = ("star", "role", "orbiters")

    def __init__(self, star, role, orbiters=None):
        self.star = star
        self.role = role
        self.orbiters = [] if orbiters is None else orbiters

    def __iter__(self):
        """Unpacks as `star, orbiters`, the `(star, planets)` pair a pass over a system has always worked on."""
        yield self.star
        yield self.orbiters

    def __repr__(self):
        return f"Host({self.role}, {getattr(self.star, 'name', self.star)!r}, {len(self.orbiters)} bodies)"
