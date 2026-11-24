# planetgen/physics/comet_coma.py

"""
The coma of an interstellar comet, read from the ices it was drawn with (GEN.1227)
=================================================================================

A comet's colours come from the species of its coma: the carbon-chain radicals CN and C2 give the green of most comets, and a comet with CO shows the blue ion tail of CO+. The volatiles that set them
are drawn when the nucleus is made (`comet_population.draw_volatiles`, in percent of the water by number, and the carbon-chain-depleted flag) and sit in the `detail` of its makeup, but only the three
largest components of the nucleus were stored. This module reads three things off the drawn volatiles, with no draw of its own and no threshold of its own, so the numbers the generator draws do not move:

* the **class**: which of the two carbon oxides is the more abundant, `CO2-dominated` or `CO-dominated` (the two larger classes of the coma recipe, rows k1 and k2 of the atmosphere recipes; the third, comparable
  amounts, is the measure-zero case of two equal drawn numbers and is read as CO2-dominated);
* the **dominant volatiles**: the species after water, from the most abundant, as many as `DOMINANT_COUNT` (the comet's three largest ices are named the same way by `small_body_composition.component_names`);
* the **carbon-chain-depleted** flag as it was drawn (a depleted comet has little CN and C2 and so no green coma).

Use:

    coma = coma_of(comet.makeup)     # a Coma: coma_class, volatiles, carbon_chain_depleted
    coma.text()                      # "CO2, CO, CH3OH"
"""

from dataclasses import dataclass

CO2_DOMINATED = "CO2-dominated"
"""str: The coma class of a comet whose CO2 is at least its CO."""

CO_DOMINATED = "CO-dominated"
"""str: The coma class of a comet whose CO is more than its CO2."""

COMA_CLASSES = (CO2_DOMINATED, CO_DOMINATED)
"""tuple: The coma classes."""

DOMINANT_COUNT = 3
"""int: How many volatiles after water are named as the dominant ones (GEN.1227): a count of what is listed, not a physical quantity, the same as the three largest components a nucleus names (`small_body_composition.component_names`)."""

SEPARATOR = ", "
"""str: What the stored list of species is joined with."""


class ComaError(ValueError):
    """The makeup carries no drawn volatiles (it is not a comet's)."""


@dataclass(frozen=True)
class Coma:
    """
    What a comet's drawn volatiles say about its coma.

    Attributes:
        coma_class (str): `CO2_DOMINATED` or `CO_DOMINATED`.
        volatiles (tuple): The formulas of the dominant volatiles after water, the most abundant first.
        carbon_chain_depleted (bool): The comet's flag: little CN and C2.
    """

    coma_class: str
    volatiles: tuple
    carbon_chain_depleted: bool

    def text(self):
        """The dominant volatiles as the one line the comet's row keeps."""
        return SEPARATOR.join(self.volatiles)


def coma_of(makeup):
    """
    The coma of a comet from its drawn makeup.

    Args:
        makeup (SmallBodyMakeup): A comet's, whose `detail` holds `volatiles` (percent of water by number, by formula) and `carbon_chain_depleted`.

    Returns:
        Coma: Its class, its dominant volatiles and the carbon-chain flag.

    Raises:
        ComaError: The makeup has no `volatiles` or no carbon oxides in them.
    """
    detail = makeup.detail
    percent = detail.get("volatiles")
    if not percent or "CO" not in percent or "CO2" not in percent:
        raise ComaError("the makeup holds no drawn volatiles with the carbon oxides")
    coma_class = CO2_DOMINATED if percent["CO2"] >= percent["CO"] else CO_DOMINATED
    ranked = sorted(percent, key=lambda formula: (-percent[formula], formula))
    return Coma(coma_class, tuple(ranked[:DOMINANT_COUNT]), bool(detail.get("carbon_chain_depleted", False)))
