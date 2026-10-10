# planetgen/galaxy/uid.py

"""
A sector's ID (GEN.69): its designation (`geometry.provisional_sector_designation`),
ring, biased layer and slot packed into one integer. It needs no row, so a
sector nobody has generated has its ID already; `sector_uid` is the integer,
`format_sector_uid` the hex the site prints. Every other object's ID is the
80-bit birth-location ID of `galaxy/object_uid.py` (GEN.170).
"""

from planetgen.galaxy.geometry import parse_sector_designation, provisional_sector_designation


def sector_uid(ring_index, layer_index, slot_index):
    """A sector's ID: its designation as an integer. Raises `ValueError`
    like `provisional_sector_designation` for an address outside the grid."""
    return int(provisional_sector_designation(ring_index, layer_index, slot_index), 16)


def sector_address(uid):
    """The `(ring_index, layer_index, ring_slot_index)` a sector ID names."""
    return parse_sector_designation(f"{uid:X}")


def format_sector_uid(uid):
    """A sector ID as the uppercase hex its designation is printed in."""
    return f"{uid:X}"
