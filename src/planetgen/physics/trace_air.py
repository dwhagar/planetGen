# planetgen/physics/trace_air.py

"""
The trace gases of a finished air (GEN.524 part 1, table m of the possibility tables)
====================================================================================

`mixture_trace` draws which of the trace gases of table m an air holds, given the life stage; until GEN.524 nothing called it. This module builds
the quantities its rows read from the finished body (the oxygen the life has built, the air's character, the life stage), draws the rows on a
stream of the body's own and gives the air the gases they draw present: nitrous oxide, phosphine, nitrogen dioxide and methyl chloride, each
at a mole fraction of its design bracket (`mixture_trace.ABUNDANCES`). The other gases keep their ratios and the total pressure stays; the
trace gases are at most one part in a million, so the temperature, the greenhouse and the colour do not move and nothing that reads the air is
made again.

Rows m3, m4 and m11 are design flags (GEN.525, `physics.design_flags`): off unless life or lightning or stellar energetic particles allow them,
and never tuned. No code models lightning or the stellar energetic particle flux of an air yet, so both are false here and the rows run only
through the life stage; a sourced rate would replace them. HCl and HF are not rows of the table (GEN.524 part 2): the halogen budget of the outgassing gives them.

The words of the rows are read from the air without a number of their own: an air is *CO2-rich* when carbon dioxide is its largest gas and
*H2-poor* when hydrogen is below the share the air's text calls a trace (`atmosphere.TRACE_FRACTION`).

A body that is drawn again (an edit, `validation.reapply_life`) loses the trace gases it had before the new ones are drawn, so applying twice gives the same air.
"""

from planetgen.generation.possibility_table import Subject
from planetgen.physics import atmosphere, mixture_subject, mixture_trace, ozone
from planetgen.util import draw

TRACE_GASES = tuple(sorted({formula for row in mixture_trace.WIRED_ROWS for formula in mixture_trace.ABUNDANCES[row]}))
"""tuple: The formulas this module adds to an air (the gases of the rows `mixture_trace.WIRED_ROWS`)."""

LIGHTNING = False
SEP_REDUCING = False
"""bool: The two design flags no model defines (see the module text); both false until a lightning rate and a stellar energetic particle flux are sourced."""


def _strip(partials):
    """`partials` (formula to Pa) without the trace gases, renormalised to the same total pressure."""
    total = sum(partials.values())
    kept = {formula: pascal for formula, pascal in partials.items() if formula not in TRACE_GASES}
    left = sum(kept.values())
    if not kept or left == total:
        return dict(kept)
    return {formula: pascal * total / left for formula, pascal in kept.items()}


def subject_of(body, partials):
    """The `Subject` of table m for `body` and its air `partials` (formula to Pa, without the trace gases)."""
    pressure_pa = sum(partials.values())
    quantities = mixture_subject.quantities_of_body(
        body, ts_k=body.surface_temperature, pressure_bar=pressure_pa / 1.0e5, star_distance_au=getattr(body, "_star_distance_au", None))
    total = pressure_pa or 1.0
    largest = max(partials, key=partials.get)
    quantities.update({
        "o2_pal": ozone.o2_level_pal(partials.get("O2", 0.0), body.gravity),
        "co2_rich": largest == "CO2",
        "h2_poor": partials.get("H2", 0.0) / total < atmosphere.TRACE_FRACTION,
        "has_solid_ground": getattr(body, "body_type", "t") == "t",
        "gas_has_co_h2s": partials.get("CO", 0.0) > 0.0 and partials.get("H2S", 0.0) > 0.0,
        "life_chemical": getattr(body, "life_stage", None) is not None,
        "lightning": LIGHTNING,
        "sep_reducing": SEP_REDUCING,
    })
    quantities.setdefault("host_type", None)                     # a rogue body has no host
    quantities.setdefault("t_layer_k", body.surface_temperature)
    needed = mixture_trace.required_trace_inputs("terrestrial")
    return Subject("terrestrial", {name: quantities[name] for name in needed})


def apply(body):
    """
    Gives the air of a rocky planet or moon the trace gases table m draws for it (and takes away the ones an earlier draw gave it). Run after the
    life stage, because the rows read it. Draws on a stream of the body's own, so the run's draws do not move. Returns the formulas added.
    """
    if getattr(body, "body_type", "t") == "g":
        return ()
    partials = _strip(atmosphere.gas_partials_pa(body))
    if not partials:
        return ()
    subject = subject_of(body, partials)
    with draw.bound(f"trace-gases:{body.is_moon}:{body.mass!r}:{body.radius!r}:{body.distance!r}"):
        amounts, _, _ = mixture_trace.draw_trace_amounts(subject)
    if amounts:
        total = sum(partials.values())
        share = sum(amounts.values())
        partials = {formula: pascal * (1.0 - share) for formula, pascal in partials.items()}
        partials.update({formula: fraction * total for formula, fraction in amounts.items()})
    atmosphere.store_air(body, partials, body.mantle_delta_iw)
    return tuple(sorted(amounts))
