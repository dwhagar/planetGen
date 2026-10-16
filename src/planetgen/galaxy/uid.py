# planetgen/galaxy/uid.py

"""
A sector's ID (GEN.69): its designation (`spherical.provisional_sector_designation`),
shell, biased band and wedge packed into one integer. It needs no row, so a
sector nobody has generated has its ID already; `sector_uid` is the integer,
`format_sector_uid` the hex the site prints. Every other object's ID is the
96-bit birth-location ID of `galaxy/object_uid.py` (GEN.170).
"""

from planetgen.galaxy.spherical import parse_sector_designation, provisional_sector_designation


UNPLACED_SECTOR_BASE = 1 << 56
"""int: A sector with no grid address has no designation. Its ID is this plus its row id, above any real
designation (which fits 56 bits), written once when the sector is saved and never changed."""


def unplaced_sector_uid(sector_row_id):
    """The ID of a sector that has no grid address."""
    return UNPLACED_SECTOR_BASE + int(sector_row_id)


def sector_uid(shell_index, band_index, wedge_index):
    """A sector's ID: its designation as an integer. Raises `ValueError`
    like `provisional_sector_designation` for an address outside the grid."""
    return int(provisional_sector_designation(shell_index, band_index, wedge_index), 16)


def sector_address(uid):
    """The `(shell_index, band_index, wedge_index)` a sector ID names."""
    return parse_sector_designation(f"{uid:X}")


def format_sector_uid(uid):
    """A sector ID as the uppercase hex its designation is printed in."""
    return f"{uid:X}"
