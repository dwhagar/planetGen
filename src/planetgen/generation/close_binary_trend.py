# generation/close_binary_trend.py

"""
The metallicity trend of the close companions (GEN.103, Boss 2026-10-16): metal-poor stars have more companions inside
10 AU than metal-rich ones (Moe, Kratter and Badenes 2019). The factor is a leaf of its own, so the tree draw
(`generation.multistar`) and the place's terms (`generation.place_terms`) can both read it.
"""

from planetgen import tuning


def close_companion_factor(fe_h):
    """
    `F(fe_h) / F(0)`, the factor on the odds of a star's companion inside 10 AU (`scaled_close_chance`): `F` is the
    straight line through `tuning.CLOSE_BINARY_FRACTION_BY_FE_H`, held flat past its ends. 1 for a star with no [Fe/H].
    """
    if fe_h is None:
        return 1.0
    return _fraction_at(fe_h) / _fraction_at(0.0)


def scaled_close_chance(chance, fe_h):
    """
    A close-companion chance moved by the star's [Fe/H] (Boss, 2026-10-16, 7:19 PM EDT): the factor
    `close_companion_factor` multiplies the *odds* `p / (1 - p)`, not `p`, for every primary mass. The result never passes 1,
    and a star already likely to have a companion rises less in absolute terms: a G dwarf (0.22) goes 0.43, 0.36, 0.25
    and 0.12 at [Fe/H] -3.0, -1.0, -0.2 and +0.5, inside Moe, Kratter and Badenes 2019's ranges, and an O star (0.45)
    goes to 0.58 at -0.7, inside the "no change within about 30%" the O and early-B surveys find.
    """
    if chance <= 0.0 or chance >= 1.0:
        return min(max(chance, 0.0), 1.0)
    odds = chance / (1.0 - chance) * close_companion_factor(fe_h)
    return odds / (1.0 + odds)


def _fraction_at(fe_h):
    points = sorted(tuning.CLOSE_BINARY_FRACTION_BY_FE_H.items())
    if fe_h <= points[0][0]:
        return points[0][1]
    for (low_x, low_y), (high_x, high_y) in zip(points, points[1:]):
        if fe_h <= high_x:
            return low_y + (high_y - low_y) * (fe_h - low_x) / (high_x - low_x)
    return points[-1][1]
