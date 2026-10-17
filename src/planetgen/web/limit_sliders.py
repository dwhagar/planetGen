# planetgen/web/limit_sliders.py

"""
The Generate page's limit sliders (UX.116)
==========================================

Every scatter and fill limit on the Generate page is one control: a slider with a number box at its end that shows the
value and accepts a typed one (`templates/generate.html`'s `limit_slider` macro, `static/limitsliders.js`). A
`LimitSlider` is what the macro and the form parser share: the control's form field, label, the stops its slider moves
between and the default. The box is the form field; the slider only moves the box, so the page still works typed into
without JavaScript.

Luminosities and system counts run over decades, so their sliders are logarithmic: the stops are a 1-2-5 ladder (the
baseline luminosity floor keeps its own presets, `generation.luminosity_floor.PRESETS`) and a typed value between two
stops is kept as typed. A mass limit is a linear slider (the baseline's presets, the core's in steps of 0.2 solar
masses).
"""

import bisect
import json
import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.generation import luminosity_floor


def ladder(low, high):
    """The 1-2-5 stops from `low` to `high`, both ends included."""
    stops = {float(low), float(high)}
    exponent = math.floor(math.log10(low))
    while 10.0 ** exponent <= high:
        for multiple in (1, 2, 5):
            value = multiple * 10.0 ** exponent
            if low <= value <= high:
                stops.add(value)
        exponent += 1
    return tuple(sorted(stops))


def steps(low, high, step):
    """The evenly spaced stops from `low` to `high`."""
    count = round((high - low) / step)
    return tuple(round(low + index * step, 10) for index in range(count + 1))


def plain(value):
    """A number as the form writes it: no exponent, no trailing zeros."""
    return format(float(value), ".10g")


@dataclass(frozen=True)
class LimitSlider:
    """
    One limit's control.

    Attributes:
        name (str): The form field (the number box's `name`).
        label (str): The label on the control and in form errors.
        scale (str): "log" or "linear": how the slider's stops are spaced.
        stops (tuple): The values the slider moves between, ascending; the first and last are the control's range.
        default (float or None): The default; `None` for a limit that follows another until it is set.
        snap (bool): Whether a typed value is moved to the nearest stop (a limit that takes only presets).
    """
    name: str
    label: str
    scale: str
    stops: tuple
    default: float | None = None
    snap: bool = False

    @property
    def low(self):
        return self.stops[0]

    @property
    def high(self):
        return self.stops[-1]

    def position(self, value):
        """The index of the stop nearest `value`, the slider's position."""
        value = min(max(float(value), self.low), self.high)
        index = bisect.bisect_left(self.stops, value)
        if index == 0:
            return 0
        below, above = self.stops[index - 1], self.stops[index]
        if self.scale == "log":
            return index if math.log(value / below) > math.log(above / value) else index - 1
        return index if value - below > above - value else index - 1

    @property
    def stops_json(self):
        return json.dumps([float(stop) for stop in self.stops])

    @property
    def tick_positions(self):
        """The stops the slider marks: every stop of a short ladder, each decade of a long one."""
        if len(self.stops) <= 12:
            return list(range(len(self.stops)))
        return [index for index, stop in enumerate(self.stops)
                if self.scale != "log" or abs(math.log10(stop) - round(math.log10(stop))) < 1e-9]


MASS_LIMIT = LimitSlider("phenomenon_min_mass", "Mass limit (solar masses)", "linear",
                         tuple(float(preset) for preset in tuning.PHENOMENON_MIN_MASS_PRESETS),
                         tuning.PHENOMENON_MIN_MASS_SOLAR, snap=True)
BRIGHT_FLOOR = LimitSlider("bright_min_luminosity", "Bright stars from (solar luminosities)", "log",
                           tuple(float(preset) for preset in luminosity_floor.PRESETS),
                           tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL)
CORE_FLOOR = LimitSlider("core_min_luminosity", "Core: bright stars from (solar luminosities)", "log",
                         ladder(tuning.WD_LUMINOSITY_RANGE_SOL[1], tuning.BRIGHT_STAR_FLOOR_MAX_SOL))
CORE_MASS = LimitSlider("core_mass_limit", "Core: mass limit (solar masses)", "linear", steps(8.0, 24.0, 0.2))
HALO_FLOOR = LimitSlider("halo_min_luminosity", "Halo: bright stars from (solar luminosities)", "log",
                         ladder(*tuning.HALO_MIN_LUMINOSITY_RANGE_SOL), tuning.HALO_MIN_LUMINOSITY_SOL)
CELL_CEILING = LimitSlider("cell_star_ceiling", "Cell star ceiling (systems)", "log", ladder(1.0, 1.0e8),
                           tuning.CELL_STAR_CEILING)
FILL_CAP = LimitSlider("fill_luminosity_cap", "Fill luminosity cap (solar luminosities)", "log",
                       ladder(tuning.WD_LUMINOSITY_RANGE_SOL[1], tuning.BRIGHT_STAR_FLOOR_MAX_SOL),
                       tuning.FILL_LUMINOSITY_CAP_SOL)

SLIDERS = {slider.name: slider for slider in
           (MASS_LIMIT, BRIGHT_FLOOR, CORE_FLOOR, CORE_MASS, HALO_FLOOR, CELL_CEILING, FILL_CAP)}
"""dict: Every limit slider by form field."""
