# planetgen/util/checks.py

"""
Checks
======

`finite_domain`, the decorator that gives a numeric physics helper one
contract: it returns finite numbers or raises `ValueError`.
"""

import functools
import inspect
import math



def _numeric_leaves(value):
    if isinstance(value, (tuple, list)):
        for item in value:
            yield from _numeric_leaves(item)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield value


def finite_domain(*, allow_inf=(), clamped=()):
    """
    Decorator giving a numeric physics helper one explicit contract: it
    returns finite numbers or raises `ValueError` -- never a
    `ZeroDivisionError`/`OverflowError` from an out-of-range input, and
    never a silent NaN/inf that would otherwise surface far away (in a
    stored row, a rendered page).

    * Every float argument must be finite, except those named in
      `allow_inf` (which may be +/-inf, e.g. `vis_viva_speed_kms`'s
      parabolic `semi_major_axis_au = math.inf`) and those named in
      `clamped` (which the helper clamps into a valid range itself, so any
      non-NaN value is fine; NaN still surfaces through the result check).
    * A `ZeroDivisionError`/`OverflowError` raised inside (a zero or
      subnormal divisor, a power overflowing) becomes a `ValueError`.
    * A non-finite number anywhere in the result becomes a `ValueError`.

    Plain unit conversions and formatting helpers deliberately don't use
    this: they pass non-finite values through unchanged.
    """
    allow_inf = frozenset(allow_inf)
    clamped = frozenset(clamped)

    def decorate(func):
        params = list(inspect.signature(func).parameters)
        where = func.__name__

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for name, value in (*zip(params, args), *kwargs.items()):
                if isinstance(value, float) and not math.isfinite(value):
                    if name in clamped or (name in allow_inf and math.isinf(value)):
                        continue
                    raise ValueError(f"{where}: {name} must be a finite number, got {value!r}")
            try:
                result = func(*args, **kwargs)
            except (ZeroDivisionError, OverflowError) as exc:
                raise ValueError(f"{where}: input out of range ({exc})") from exc
            for leaf in _numeric_leaves(result):
                if not math.isfinite(leaf):
                    raise ValueError(f"{where}: result out of range ({leaf!r}) for {args!r}")
            return result

        return wrapper

    return decorate
