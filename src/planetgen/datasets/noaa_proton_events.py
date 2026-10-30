# planetgen/datasets/noaa_proton_events.py

"""
The NOAA solar proton event list, 1976 to 2017 (GEN.747): every event whose >10 MeV proton flux reached 10 particle flux units
(pfu, protons per cm^2 per second per steradian) at the GOES spacecraft, with the start, the time of maximum and the peak flux,
and the flare that went with it.

The file is `seps.html`, NOAA Space Environment Services Center's "Solar Proton Events Affecting the Earth Environment"
preliminary listing as NASA's Solar Data Analysis Center serves it (https://umbra.nascom.nasa.gov/SEP/seps.html). It is an HTML
table with a year row before each year's events. The reader keeps the year, start, maximum, peak flux and the X-ray/optical
importance of the flare. Start and maximum are datetimes: the listing gives day and hour without a year, so the year of the
heading above the row is used, and a maximum earlier in the calendar than its start belongs to the next year. Ten of the newest
rows (2012 to 2017) write the date as `Aug/10 0435` or `June 19/0450` or have a non-breaking space; the reader reads those too, so it
finds 265 events where the 255-event table of the stand-in report of 10-24 (`reports/best-guess-standins-2026-10-24-files`) left
those ten out, among them the 1,490 pfu event of September 2017.

Until now the solar particle event rate of the radiation model was one number (`SEP_FREE_MSV_YR`) and an exponential in
rigidity.

The reader raises `fetch.DataMissing` when the file has not been downloaded.
"""

import datetime
import functools
import html
import re
from dataclasses import dataclass

from planetgen.datasets import fetch

FILE = "seps.html"

DATASET = fetch.Dataset(
    key="noaa_proton_events",
    title="NOAA solar proton events above 10 MeV, 1976 to 2017",
    source="NOAA Space Environment Services Center, Solar Proton Events Affecting the Earth Environment, via NASA SDAC",
    pieces=(fetch.Piece(FILE, "https://umbra.nascom.nasa.gov/SEP/seps.html",
                        "8bf7e4ed2972c333d43373a9eee4455d3362ce859e2fe20704544da18ea0d321", 45002),),
    replaces="the one solar particle dose and the exponential in rigidity of the radiation model",
    fill="SC5",
    fallback="the stand-in dose and rigidity law of physics/radiation.py",
)

_YEAR = re.compile(r"<strong>(\d{4})</strong>")
_MOMENT = re.compile(r"^([A-Za-z]{3})[A-Za-z]* (\d{1,2}) (\d{4})$")
_MONTHS = {name: number for number, name in enumerate(("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 1)}


@dataclass(frozen=True)
class Event:
    """One proton event: when it began and peaked (UTC), the peak >10 MeV flux in pfu and the flare's X-ray/optical class."""

    start: datetime.datetime
    maximum: datetime.datetime
    peak_flux_pfu: int
    flare: str


def _normalise(text):
    """`Month day hhmm` with single spaces, from the three spellings the listing uses."""
    return re.sub(r"\s+", " ", html.unescape(text).replace("\xa0", " ")).strip().replace("/", " ")


def _moment(text, year):
    """The datetime of `Apr 30/2120`, `Aug/10 0435` or `June 19/0450` in `year`."""
    match = _MOMENT.match(_normalise(text))
    if not match:
        raise ValueError(f"not a date and hour of a proton event: {text!r}")
    hhmm = match.group(3)
    return datetime.datetime(year, _MONTHS[match.group(1).title()], int(match.group(2)), int(hhmm[:2]), int(hhmm[2:]))


def parse(text):
    """The events of a `seps.html`, in the listing's order."""
    events, year = [], None
    for row in re.split(r"<tr>", text):
        heading = _YEAR.search(row)
        if heading:
            year = int(heading.group(1))
            continue
        cells = [re.sub(r"<[^>]+>", " ", cell) for cell in re.split(r"<td[^>]*>", row)[1:]]
        if len(cells) < 3 or not _MOMENT.match(_normalise(cells[0])):
            continue
        if year is None:
            raise ValueError("a proton event before any year heading")
        start, maximum = _moment(cells[0], year), _moment(cells[1], year)
        if maximum < start:
            maximum = maximum.replace(year=year + 1)
        flux = re.sub(r"[^0-9]", "", cells[2])
        if not flux:
            raise ValueError(f"a proton event without a peak flux: {cells[:3]!r}")
        flare = re.sub(r"\s+", " ", cells[6]).strip() if len(cells) > 6 else ""
        events.append(Event(start, maximum, int(flux), flare))
    if not events:
        raise ValueError("a proton event listing with no events")
    return tuple(events)


@functools.lru_cache(maxsize=2)
def _load(path):
    with open(path, encoding="utf-8", errors="replace") as handle:
        return parse(handle.read())


def events():
    """The 265 events of the listing, 1976 to 2017, oldest first."""
    return _load(fetch.locate(DATASET).require().path(FILE))
