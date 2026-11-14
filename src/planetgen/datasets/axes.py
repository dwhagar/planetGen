# planetgen/datasets/axes.py
"""The arithmetic the dataset readers share on a tabulated axis: where a value sits between two nodes of a grid, and the
integral and cell means of a curve given at nodes (a spectrum, a cross section, an absorption coefficient)."""

import bisect

import numpy as np


def bracket(values, value):
    """`(lower, upper, weight)`: the indices of the two ascending `values` around `value` (held to the ends of the list) and the
    share of the upper one; `lower == upper` and a weight of 0 when `value` is on a node or past an end."""
    value = min(max(value, values[0]), values[-1])
    upper = min(bisect.bisect_left(values, value), len(values) - 1)
    lower = upper if values[upper] == value else upper - 1
    if lower == upper:
        return lower, upper, 0.0
    return lower, upper, (value - values[lower]) / (values[upper] - values[lower])


def cumulative_integral(x, y):
    """The trapezoid integral of `y` over `x` from `x[0]` to each node: `len(x)` values, the first 0."""
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


def cell_means(x, y, edges):
    """The mean of the piecewise-linear curve `y(x)` over each cell of `edges`; NaN for a cell the table does not cover."""
    cumulative = cumulative_integral(x, y)
    low, high = edges[:-1], edges[1:]
    inside = (low >= x[0]) & (high <= x[-1])
    area = np.interp(high, x, cumulative) - np.interp(low, x, cumulative)
    return np.where(inside, area / (high - low), np.nan)
