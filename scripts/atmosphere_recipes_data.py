# scripts/atmosphere_recipes_data.py

"""
The machine columns of the atmosphere recipe tables (GEN.478)
=============================================================

Part 09 of the atmosphere research (docs/design/composition/atmosphere-composition-parts/09-weighted-composition-tables.md) gives the
104 rows in prose: the option, the gate, the weight, its source, its status and the notes. `build_atmosphere_recipes.py`
reads those words from the part itself and takes from here what the part leaves as text: for each row its kind, the classes
and conditions of its gate, the weight as an expression the loader evaluates, the distribution of a continuous quantity,
the constants a row draws within a bracket, and the corrections the later research made.

Expressions are Python arithmetic over the stored quantities of the subject (the names of `INPUTS` in
`planetgen.physics.atmosphere_recipes`), the row's constants and these functions: `phi` (normal CDF), `window` (a soft
window, lo and hi may be None), `exp log log10 sqrt min max abs clamp loglog`, `w('id')` (the gated weight of a sibling
row), `p('id')` (the probability of a row in its own table at this subject) and `pdist('id', lo, hi)` (the probability a
continuous row's value lies in [lo, hi]).
"""

import json
import os

ROCKY = ["terrestrial"]   # a rocky or icy body that is not a giant: the kind (`body_type` "t") the cut-over's classifier reads (GEN.840)
GIANT = ["giant"]
BD = ["brown-dwarf"]
AIR = ROCKY               # a row about a body's air also lists the surface-pressure condition `AIR_WHEN`
AIR_MIN_BAR = 1e-4        # 10 Pa: below it a body is airless (the classifier's D, S and C rows use the same line)
AIR_WHEN = [{"quantity": "pressure_bar", "minimum": AIR_MIN_BAR}]
ICE_GIANT_MAX_EARTH = 35  # the classifier's ice-giant row (`planet_classes.json`, class I) ends here
BELOW_ICE_GIANT = [{"quantity": "mass_earth", "maximum": ICE_GIANT_MAX_EARTH}]   # a giant of this mass or less has a rock and water core (table b)
ANY = ["*"]


def c(quantity, minimum=None, maximum=None, equals=None, one_of=None):
    """A condition of a gate."""
    out = {"quantity": quantity}
    if minimum is not None:
        out["minimum"] = minimum
    if maximum is not None:
        out["maximum"] = maximum
    if equals is not None:
        out["equals"] = equals
    if one_of is not None:
        out["one_of"] = list(one_of)
    return out


def drawn(group, *ids):
    """The gate of a row that is drawn only under a parent row of an earlier group."""
    return c(f"drawn.{group}", one_of=ids) if len(ids) > 1 else c(f"drawn.{group}", equals=ids[0])


_WATER_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "planetgen", "physics", "data", "ngpps_formation_water.json")
with open(_WATER_FILE, encoding="utf-8") as _handle:
    _WATER = json.load(_handle)


def _water_logit(k):
    """The logit of the k-th nested event of the formation-distance law (GEN.802), as an expression of lx, lm and hinge."""
    a, b, c, h = _WATER["model"]["coefficients"][k]
    return f"{a} + {b}*lx + {c}*lm + {h}*hinge"


_WATER_X, _WATER_M = _WATER["model"]["x_range"], _WATER["model"]["mass_range_earth"]
WATER_LETS = {
    "lx": f"log(clamp(snow_ratio, {_WATER_X[0]}, {_WATER_X[1]}))",
    "lm": f"log10(clamp(mass_earth, {_WATER_M[0]}, {_WATER_M[1]}))",
    "hinge": f"max(0, lx - log({_WATER['model']['knot']}))",
    "nd": f"1/(1 + exp(-({_water_logit(0)})))",
    "hi": f"min(nd, 1/(1 + exp(-({_water_logit(1)}))))",
    "ri": f"min(hi, 1/(1 + exp(-({_water_logit(2)}))))",
}
"""dict: The lets of rows b1 to b4: `nd` is the probability that the planet is not dry, `hi` that it holds at least 0.16 percent water, `ri` that it
holds at least 10 percent (the nested events of `ngpps_formation_water.json`, held in order). The four weights are the class shares."""

WATER_LIMIT = ("a MODEL label: the shares are the Bern population synthesis against the distance the planet formed at (a_form / a_snow, with the orbit "
               "standing for it because the generator does not migrate planets) and its mass, not an observation; the class thresholds mapping to Ortenzi's "
               "cases are an assumption; the plateau c_b4 (1 is the law as fitted; the water-rich planets it removes become 'high') is bracketed 0 to 1 because "
               "published water-rich shares run from at most 3 percent (Rogers 2025) to 87 percent; the law is not tuned to the Luque and Palle densities (GEN.510)")
WATER_CORRECTION = ("Replaced 2026-10-25 (GEN.802): the weights (1-s)^2, 2s(1-s), s^2 with s = 1/(1 + exp(-ln x/0.4)) and the water-rich mass logistic of GEN.510 "
                    "(centre 6.1 Earth masses) with the gate x >= 1 are one law now: the planets of the population that formed inside the snow line are 26 to 60 percent "
                    "water-rich at x = 0.55 to 0.75, and the mass dependence at fixed x is weak (the mass law was a migration effect, counted once here).")


# The tables that compete for 100 percent: a group is drawn once, from the rows of it that pass the subject.
# `fallback` is the row that wins when every weight is zero ("none" is the explicit outcome that no row of the group
# applies). `constants` are the design constants the group's rows draw once inside a bracket, [low, high, shape].
GROUPS = {
    "a": {"rows": ["a5", "a6", "a7"], "fallback": None},
    "b": {"rows": ["b1", "b2", "b3", "b4"], "fallback": None,
          "constants": {"c_b4": [0.0, 1.0, "uniform"]}},
    "c": {"rows": ["c1", "c2", "c3", "c4"], "fallback": None},
    "c-after": {"rows": ["c5", "c6"], "fallback": "c6"},
    "d": {"rows": ["d1", "d2", "d3", "d4", "d5"], "fallback": "d1"},
    "e": {"rows": ["e1", "e2", "e3", "e4"], "fallback": "e1",
          "constants": {"f_env0": [0.05, 0.2, "log"], "eta": [0.1, 0.5, "log"]}},
    "f": {"rows": ["f1", "f2", "f3", "f4"], "fallback": None},
    "g": {"rows": ["g1", "g2", "g3", "g4", "g5"], "fallback": None},
    "h": {"rows": ["h1", "h2", "h3"], "fallback": None},
    "i": {"rows": ["i1", "i2", "i3", "i4", "i5"], "fallback": None},
    "j": {"rows": ["j1", "j2", "j3", "j4", "j5"], "fallback": "none",
          "constants": {"c_j1": [0.0, 1.0, "uniform"], "c_j5": [0.05, 0.5, "log"]}},
    "k": {"rows": ["k1", "k2", "k3"], "fallback": "none"},
    "l": {"rows": ["l1", "l2", "l3"], "fallback": "l2", "constants": {"c_l3": [0.0, 0.1, "uniform"]}},
    "p": {"rows": ["p1", "p2", "p3", "p4", "p5"], "fallback": None},
}

# The Berta-Thompson, Wachiraphan and Murray 2025 shoreline (arXiv:2507.02136), posteriors of the main sample, read in
# docs/design/composition/tunable-gaps-sources/group2-delivery-retention.md: log10 f0 2.62, p 5.9, q 1.17, w 0.247 dex.
_DELTA = "log10(insolation) - (2.62 + 5.9*log10(v_esc_kms/11.186) + 1.17*log10(luminosity_solar))"
_P_ATM = "1/(1 + 10**(delta/0.247))"
_SHORE = {"delta": _DELTA, "p_atm": _P_ATM}

RETIRED = frozenset({"m5", "m6"})
"""frozenset: The rows of the part-09 table that the generator no longer has (GEN.524 part 2: HCl and HF come from the halogen budget of the outgassing, not
from a presence draw); the part keeps them as the record, the builder leaves them out of the data file."""

MACHINE = {
    # ---- a. mantle redox ----
    "a1": dict(state="partly verified", kind="param", classes=ROCKY, sets="redox_diw",
               limit="the surface redox of the magma ocean from physics.magma_ocean_redox (stand-in rows MR1 to MR7): the base pressure is a power of the mass, held at the 140 GPa where the "
                     "first-principles calculations end (5.4 Earth masses; no source above), the ferric-share power (0.5 to 1.0) and the offset are read from two papers whose equations were not "
                     "available when it was built, there is no temperature or composition term, and the graphite-flotation alternative of Gaillard 2022 is not drawn (GEN.482, GEN.1005)",
               source="Hirschmann 2022 (GCA 328, section 6.3 and Table 1): Fe3+/FeT of the ocean surface 0.034 to 0.103 and the dIW it gives; Deng 2020: the ferric share of the melt against the "
                      "pressure at the base of the ocean (factor of two from 55 to 25 GPa) and the base pressures; Sossi 2020 0.037 Fe3+/sumFe (READ, COMPUTED); the form, the pressure "
                      "scaling and the width are DESIGN brackets of docs/design/composition/magma-ocean-redox-model.md",
               correction="Corrected 2026-10-23 (GEN.479): the gate M >= 0.01 left a body under 0.01 Earth masses with no redox row at all. GEN.1005 (2026-10-30): the mean and width are the "
                          "model of physics.magma_ocean_redox in place of the line through three Solar System bodies, which are now tests of the model.",
               dist={"kind": "normal", "mu": "magma_ocean_diw(mass_earth)", "sigma": "magma_ocean_diw_sigma(mass_earth)"}),
    "a2": dict(kind="param", classes=AIR, when=[c("long_lived_geology", equals=True)], sets="redox_diw_later",
               state="partly verified",
               limit="the magma-ocean surface of the body (physics.magma_ocean_redox, at most that of one Earth mass) plus an oxidation increment of 1.2 to 2.4 that grows with the mass: a "
                     "line between two numbers (stand-in row MR7), not a model of the oxidation (hydrogen loss, lower-mantle disproportionation, recycling); the sources give a gap of 0 to 3 "
                     "that grows with mass; the later mantle is never below the magma-ocean surface of the same body (GEN.482, GEN.482b, GEN.1005)",
               source="the later upper mantle of the rocky bodies measured: dIW +3.5 for the one with plate tectonics (Deng 2020 with its refs 37, 6, 7; Hirschmann 2022 section 4: above IW+3), "
                      "IW to IW+1 for the primitive mantle of the next lighter one and IW to IW-2 for the lightest (Hirschmann 2022 introduction and section 4; Deng 2020 refs 8 to 10); the "
                      "increment, its growth with mass and the width 1.0 are DESIGN (READ, COMPUTED; docs/design/composition/magma-ocean-redox-model.md)",
               dist={"kind": "normal", "mu": "later_mantle_diw(mass_earth)", "sigma": 1.0, "range": [-6, 6]},
               correction="GEN.1005 (2026-10-30): the mean is the magma-ocean model's surface value plus the increment, in place of the straight line through three Solar System bodies; GEN.482b takes the larger of a2 and a1 in `draw_redox`."),
    "a3": dict(state="partly verified", limit="22 polluted white dwarf values (READ, COMPUTED), a selection-biased sample of core-formation redox; the scale differs from a1; the width 0.22 is tuned with the detection model (GEN.513)", kind="param", classes=["asteroid"], sets="redox_diw", tunable=True,
               dist={"kind": "normal", "mu": -1.25, "sigma": 0.22},
               source="22 polluted white dwarf values: median -1.26, mean -1.25, sample sd 0.45 of the reported medians (Doyle 2020 Table "
                      "1, Doyle 2023 Table 4, via docs/design/composition/tunable-gaps-sources/group3-redox-wd-radio.md; READ, COMPUTED); "
                      "the width 0.22 is the intrinsic spread of the 13 well-constrained stars (COMPUTED) and lies inside the "
                      "0 to 0.38 the detection model accepts (GEN.513, docs/design/composition/wd-redox-tuning.md); the meteoritic bracket "
                      "-1.4 to 0 (Schaefer and Fegley 2017 Table 2, READ) is kept as a check",
               correction="Corrected 2026-10-18 (tunable-gaps 3.2): the uniform U[-1.4, 0] has its mean 0.55 dex above the "
                          "white dwarf ensemble; the core-formation scale differs from the surface magma-ocean scale of a1."),
    "a4": dict(state="partly verified", limit="22 polluted white dwarf values (READ, COMPUTED), a selection-biased sample; the width 0.22 is tuned with the detection model, the settling is not modelled (GEN.513)", kind="param", classes=ROCKY, when=[c("host_type", equals="WD")], sets="redox_diw", tunable=True,
               dist={"kind": "normal", "mu": -1.25, "sigma": 0.22},
               source="22 polluted white dwarf values: mean -1.25 (standard error 0.10), sample sd 0.45 of the reported medians "
                      "(Doyle 2020, Doyle 2023; READ, COMPUTED); the width 0.22 is the intrinsic spread of the 13 "
                      "well-constrained stars and lies inside the 0 to 0.38 the detection model accepts (GEN.513)",
               correction="Corrected 2026-10-18 (tunable-gaps 3.2): the centre -1.0 came from the first six stars; the "
                          "22-value mean is -1.25 +/- 0.10 and the bracket [-1.37, -0.6] of six stars is dropped."),
    "a5": dict(kind="option", group="a", classes=ROCKY, weight="phi((0 - redox_diw)/0.5)"),
    "a6": dict(kind="option", group="a", classes=ROCKY, weight="1 - w('a5') - w('a7')"),
    "a7": dict(kind="option", group="a", classes=ROCKY, weight="phi((redox_diw - 2.5)/0.5)"),

    # ---- b. delivered volatile inventory ----
    "b1": dict(state="partly verified", limit=WATER_LIMIT, kind="option", group="b", classes=ROCKY + GIANT, when=BELOW_ICE_GIANT, weight="1 - nd", lets=WATER_LETS,
               correction=WATER_CORRECTION),
    "b2": dict(state="partly verified", limit=WATER_LIMIT, kind="option", group="b", classes=ROCKY + GIANT, when=BELOW_ICE_GIANT, weight="nd - hi", lets=WATER_LETS),
    "b3": dict(state="partly verified", limit=WATER_LIMIT, kind="option", group="b", classes=ROCKY + GIANT, when=BELOW_ICE_GIANT, weight="hi - ri + (1 - c_b4)*ri", lets=WATER_LETS),
    "b4": dict(state="partly verified", limit=WATER_LIMIT, kind="option", group="b", classes=ROCKY + GIANT, when=BELOW_ICE_GIANT, weight="c_b4*ri", lets=WATER_LETS),
    "b5": dict(state="stand-in", replacement="GEN.510", kind="param", classes=ROCKY + GIANT, when=BELOW_ICE_GIANT + [drawn("b", "b2", "b3")],
               dist={"kind": "loguniform", "lo": 5e-4, "hi": 5e-3}),
    "b6": dict(state="stand-in", replacement="GEN.510", kind="binary", classes=ROCKY + GIANT, when=[c("zone", equals="r"), c("mass_earth", maximum=ICE_GIANT_MAX_EARTH)],
               weight="0.8 if body_class == 'giant' else 0.35",
               dist={"kind": "loguniform", "lo": 1e-4, "hi": 0.1}),

    # ---- c. outgassed regime ----
    "c1": dict(kind="option", group="c", classes=ROCKY, when=[drawn("a", "a5")], weight="phi((ts_k - 950)/50)"),
    "c2": dict(kind="option", group="c", classes=ROCKY, when=[drawn("a", "a5"), c("pressure_bar", minimum=10)],
               weight="phi((1100 - ts_k)/50)"),
    "c3": dict(kind="option", group="c", classes=ROCKY, when=[drawn("a", "a6")], weight=1),
    "c4": dict(kind="option", group="c", classes=ROCKY, when=[drawn("a", "a7")], weight=1),
    "c5": dict(kind="option", group="c-after", classes=ROCKY,
               when=[drawn("c", "c1", "c2", "c3", "c4"), c("ts_k", maximum=647), c("ocean_present", equals=True)],
               weight=1),
    "c6": dict(kind="option", group="c-after", classes=ROCKY, weight="1 - w('c5')",
               correction="The gate Ts >= 647 K is the complement of c5's gate, so c6 is gated by its weight 1 - w[c5] and "
                          "wins wherever c5 does not: it is the table's fallback row."),

    # ---- d. retention regime ----
    "d1": dict(kind="option", group="d", classes=ROCKY + ["asteroid"], weight="1 - p_atm", lets=_SHORE,
               source="Berta-Thompson, Wachiraphan and Murray 2025 (arXiv:2507.02136): the logistic shoreline "
                      "P(atmosphere) = 1/(1 + 10^(delta/w)), log10 f0 2.62, p 5.9, q 1.17, w 0.247 dex, 917 labelled "
                      "planets (READ); the Meni-Gallardo and Palle 2025 line, slope 5.77, agrees (READ)",
               correction="Corrected 2026-10-18 (tunable-gaps 3.3): the K window LU[3.0e-4, 6.6e-4] with slope 4 rests on 13 "
                          "Solar System bodies and is contradicted by 55 Cnc e and TOI-561 b; the two-class probability of an "
                          "atmosphere now follows the published logistic, and the K window is no longer used."),
    "d2": dict(kind="option", group="d", classes=ROCKY + ["asteroid"], lets=_SHORE,
               weight="p_atm*window(delta, -0.3, None, 0.247)",
               source="share of P(atmosphere) (Berta-Thompson 2025, READ) by distance below the shoreline; the edges "
                      "-0.3, -0.8, -1.4 dex and the width 0.247 are DESIGN",
               correction="The five-way split of the atmosphere cases is constrained by no measurement (tunable-gaps 3.3); "
                          "the four windows add up to 1, so together they carry P(atmosphere)."),
    "d3": dict(kind="option", group="d", classes=ROCKY + ["asteroid"], lets=_SHORE,
               weight="p_atm*window(delta, -0.8, -0.3, 0.247)",
               source="same as d2"),
    "d4": dict(kind="option", group="d", classes=ROCKY + ["asteroid"], lets=_SHORE,
               weight="p_atm*window(delta, -1.4, -0.8, 0.247)",
               source="same as d2"),
    "d5": dict(kind="option", group="d", classes=ROCKY + ["asteroid"], lets=_SHORE,
               weight="p_atm*window(delta, None, -1.4, 0.247)",
               source="same as d2"),
    "d6": dict(state="stand-in", replacement="GEN.511", kind="derived", classes=ROCKY + ["asteroid"], table_of="d",
               table={"d1": 0.0, "d2": 1e-5, "d3": 0.01, "d4": 0.1, "d5": None}),
    "d7": dict(kind="flag", classes=GIANT, when=[c("mass_earth", minimum=1, maximum=30)],
               weight="1 if 2.016*AMU_KG*(v_esc_kms*1000)**2/(2*K_B*(teq_k + min(1500, 750*sqrt(xuv_index)))) >= 30 else 0"),

    # ---- e. sub-Neptune and super-Earth envelope ----
    "e1": dict(state="partly verified", limit="f_env0 and the onset age of the loss are tuned to the radius valley (slope -0.11 and 1.84 Earth radii at 10 d); the loss itself is the hydro-based rate of Kubyshkina et al. 2018 over the star's XUV history (GEN.822); the core-mass width, the envelope width and the host sample are design choices; the mechanism (photoevaporation or core-powered loss) is not chosen (GEN.509, GEN.722)", tunable=True, kind="option", group="e", classes=GIANT, when=[c("mass_earth", minimum=1, maximum=30)],
               weight="(1 - p('b4'))*f_lost", lets={"f_lost": "envelope_lost_fraction(mass_earth, f_env0, insolation, host_mass_solar, luminosity_solar, age_gy, teq_k)"}),
    "e2": dict(kind="option", group="e", classes=GIANT, when=[c("mass_earth", minimum=1, maximum=30)],
               weight="(1 - p('b4'))*(1 - f_lost)",
               lets={"f_lost": "envelope_lost_fraction(mass_earth, f_env0, insolation, host_mass_solar, luminosity_solar, age_gy, teq_k)"}),
    "e3": dict(state="partly verified", limit="the gate is Madhusudhan 2021 (READ); the share is the b4 stand-in (GEN.510)", kind="option", group="e", classes=GIANT, when=[c("mass_earth", minimum=1, maximum=30)],
               weight="p('b4')*(1 - eta)"),
    "e4": dict(state="partly verified", limit="the gate is Madhusudhan 2021 (READ); the share is the b4 stand-in (GEN.510)", kind="option", group="e", classes=GIANT, when=[c("mass_earth", minimum=1, maximum=30)],
               weight="p('b4')*eta"),
    "e5": dict(state="partly verified", limit="bracket edges from about 8 JWST retrievals and one 15-planet fit; no population mean exists (GEN.509)", kind="param", classes=GIANT, when=[drawn("e", "e2")],
               dist={"kind": "loguniform", "lo": 0.01, "hi": 0.6}, bins=[0.01, 0.03, 0.1, 0.3, 0.6],
               values={"mean_molecular_weight": "2.4 + 3.1*log10(value/0.01)/log10(58)"}),
    "e6": dict(kind="binary", classes=GIANT, when=[drawn("e", "e2")], weight="phi((500 - teq_k)/100)"),

    # ---- f. giant planets ----
    "f1": dict(kind="option", group="f", classes=GIANT, weight="window(teq_k, None, 200, 100)"),
    "f2": dict(kind="option", group="f", classes=GIANT, weight="window(teq_k, 200, 700, 100)"),
    "f3": dict(kind="option", group="f", classes=GIANT + BD, weight="window(teq_k, 700, 1300, 100)"),
    "f4": dict(kind="option", group="f", classes=GIANT, when=[c("zone", equals="h")],
               weight="window(teq_k, 1300, None, 100)"),
    "f5": dict(kind="param", classes=GIANT, sets="metallicity", tunable=True,
               dist={"kind": "normal", "log10": True, "mu": "log10(2.6/mass_jupiter + 3.3)", "sigma": 0.26,
                     "range": [-0.5, 2.0]},
               source="147 warm giants: Zp/Z* = 2.6/M + 3.3 (M in Jupiter masses; a +/- 0.3, b +/- 0.5), 0.03 to 12 Jupiter masses "
                      "(Chachan 2025, arXiv 2509.20428; READ); intrinsic scatter of the heavy-element mass 0.26 dex, band 0.25 to "
                      "0.27 (Thorngren 2016, N=47; READ, COMPUTED)",
               correction="Retuned 2026-10-23 (GEN.516): the centre 9.7 M^-0.45 (Thorngren 2016) is 1.6 times the 147-giant fit at one "
                          "Jupiter mass (9.7 against 5.9) and the two meet at 0.1 and 10 Jupiter masses; the later, larger sample "
                          "wins. The width 0.3 was a recalled design value, the fit gives 0.26 dex."),
    "f6": dict(kind="bins", classes=GIANT, of="f5", bins=["10**-0.5", 1, "10**0.5", 10, "10**1.5", 100]),
    "f7": dict(state="partly verified", limit="stellar C/O is measured (N=849); the planet tail is 6 of 69 (LExACoM) with a 0.2 to 0.3 method systematic (GEN.516)", kind="param", classes=GIANT + BD, sets="c_to_o", constants={"f_tail": [0.043, 0.087, "log"]}, tunable=True,
               dist={"kind": "mixture", "weight": "f_tail",
                     "core": {"kind": "normal", "log10": True, "mu": "log10(0.57)", "sigma": 0.21,
                              "range": ["log10(0.2)", "log10(0.8)"]},
                     "tail": {"kind": "uniform", "lo": 0.8, "hi": 1.5}},
               bins=[0.2, 0.35, 0.7, 0.8, 1.0, 1.5],
               source="stellar C/O median 0.47 with no star above 0.8 (N=849); planets (LExACoM, 69): median 0.57, tail "
                      "0.043 to 0.087 above 0.8 (tunable-gaps 3.4, COMPUTED from READ tables); core width 0.21 dex "
                      "reproduces the 12 of 69 planets below 0.35 (GEN.516; the single planet WASP-77Ab 0.59 +/- 0.08, Line 2021, "
                      "READ, gave 0.15)",
               correction="Corrected 2026-10-18 (tunable-gaps 3.4): the tail share LU[0.0013, 0.34] is outside the data "
                          "(planet tail 0.043 to 0.087) and the core centre is 0.57. Part 05's recipe 0.55*10^N(0, 0.15) "
                          "truncated to [0.2, 1.5] put 13.9 percent of giants above C/O 0.8, where about 2 to 5 percent "
                          "was intended: the core stops at 0.8 and the tail is explicit."),
    "f8": dict(kind="binary", classes=GIANT + BD,
               weight="phi((1100 - teff_k)/100) if body_class == 'brown-dwarf' else phi((500 - teq_k)/100)"),

    # ---- g. brown dwarfs ----
    "g1": dict(state="partly verified", limit="late M is the extension of the map above L0 (2238 K); the map has no M-type fit (Table 13) and the 2950 K cap is the BHAC15 temperature of the hydrogen-burning limit at its peak (GEN.494)", kind="option", group="g", classes=BD, weight="phi((2950 - teff_k)/100) - phi((2238 - teff_k)/134)"),
    "g2": dict(kind="option", group="g", classes=BD, weight="phi((2238 - teff_k)/134) - phi((1273 - teff_k)/79)"),
    "g3": dict(kind="option", group="g", classes=BD, weight="phi((1273 - teff_k)/79) - phi((686 - teff_k)/79)"),
    "g4": dict(state="partly verified", limit="late T (T8 to T9) is a subclass of T, not a type of its own; its edges are the map's T8 and Y0 values (Kirkpatrick 2021 Table 13, GEN.494)", kind="option", group="g", classes=BD, weight="phi((686 - teff_k)/79) - phi((460 - teff_k)/79)"),
    "g5": dict(kind="option", group="g", classes=BD, weight="phi((460 - teff_k)/79)"),
    "g6": dict(kind="binary", classes=BD, weight="phi((1100 - teff_k)/100)"),

    # ---- h. hot Jupiters ----
    "h1": dict(state="partly verified", limit="the edge is the equilibrium line of methane and carbon monoxide at the photospheric pressure, computed with the NASA polynomials at solar metals and C/O 0.55; the softness of 100 K is DESIGN (GEN.518)", kind="option", group="h", classes=GIANT, when=[c("zone", equals="h")],
               lets={"edge": "loglog(photosphere_bar, 0.001, 726, 0.01, 830, 0.1, 966, 1, 1157, 10, 1439)"},
               weight="window(tday_k, None, edge, 100)"),
    "h2": dict(state="partly verified", limit="the edges are computed equilibrium lines at the photospheric pressure (methane equal to carbon monoxide; water 10 percent destroyed) with the NASA polynomials at solar metals and C/O 0.55, gas only and no ions; the softness of 100 K is DESIGN (GEN.518)", kind="option", group="h", classes=GIANT, when=[c("zone", equals="h")],
               lets={"low_edge": "loglog(photosphere_bar, 0.001, 726, 0.01, 830, 0.1, 966, 1, 1157, 10, 1439)",
                     "high_edge": "loglog(photosphere_bar, 0.001, 2331, 0.01, 2542, 0.1, 2793, 1, 3096, 10, 3471)"},
               weight="window(tday_k, low_edge, high_edge, 100)"),
    "h3": dict(state="partly verified", limit="the edge is the computed temperature at which a tenth of the water is destroyed, at the photospheric pressure, gas only and no ions; the softness of 100 K is DESIGN (GEN.518)", kind="option", group="h", classes=GIANT, when=[c("zone", equals="h")],
               lets={"edge": "loglog(photosphere_bar, 0.001, 2331, 0.01, 2542, 0.1, 2793, 1, 3096, 10, 3471)"},
               weight="window(tday_k, edge, None, 100)"),
    "h4": dict(state="partly verified", limit="a Ti threshold between 2350 and 2641 K from N=5 planets (TiO seen at 2641 K, not at 2210, 2250, 2260 and 2350 K); no population measures Kzz, so none is drawn (GEN.519)", kind="binary", classes=GIANT, when=[c("zone", equals="h"), drawn("h", "h2", "h3")],
               constants={"c_h4": [2350.0, 2641.0, "uniform"]}, weight="1 if teq_k >= c_h4 else 0"),
    "h5": dict(kind="binary", classes=GIANT, when=[c("zone", equals="h")], weight="pdist('f7', 1.0, None)"),

    # ---- i. lava worlds ----
    "i1": dict(kind="option", group="i", classes=ROCKY, when=[c("zone", equals="h")], weight="window(tp_k, None, 2033, 60)"),
    "i2": dict(kind="option", group="i", classes=ROCKY, when=[c("zone", equals="h")], weight="window(tp_k, 2033, 2588, 60)"),
    "i3": dict(kind="option", group="i", classes=ROCKY, when=[c("zone", equals="h")], weight="window(tp_k, 2588, 2890, 60)"),
    "i4": dict(kind="option", group="i", classes=ROCKY, when=[c("zone", equals="h")], weight="window(tp_k, 2890, 3168, 60)"),
    "i5": dict(kind="option", group="i", classes=ROCKY, when=[c("zone", equals="h")], weight="window(tp_k, 3168, None, 60)"),
    "i6": dict(state="partly verified", limit="2 of 2 molten-dayside ultra-short-period planets keep volatiles (N=2): the bracket is the 95 percent interval 0.29 to 0.99 of Beta(3, 1) (GEN.511)", kind="binary", classes=ROCKY, when=[c("zone", equals="h"), c("mass_earth", minimum=1, maximum=10)],
               constants={"c_i6": [0.29, 0.99, "log"]}, weight="c_i6",
               values={"p_eq_bar": "10**(-34.7*exp(-0.00112*tp_k))"}),

    # ---- j. frost-equilibrium and evolved atmospheres ----
    "j1": dict(state="partly verified", limit="the runaway gate is Kopparapu 2013 (READ); the retained fraction is a design residual that must allow 0, since 0 of 2 bare M-dwarf planets (TRAPPIST-1 b, c) cannot be forbidden (GEN.511)", kind="option", group="j", classes=ROCKY,
               when=[c("redox_diw", minimum=1), c("water_lost", equals=True), c("runaway_ratio", minimum=1.0)],
               weight="c_j1"),
    "j2": dict(kind="option", group="j", classes=ROCKY,
               when=[c("mass_earth", maximum=0.5), drawn("d", "d3"), c("ts_k", minimum=190, maximum=270)], weight=1,
               state="derived",
               correction="Decided (Boss, 2026-10-18, fork 7): the Mars-, Titan- and Pluto-like occurrence is an output of "
                          "the escape and retention physics, not a design flag; the gate carries it and the weight is 1."),
    "j3": dict(kind="option", group="j", classes=ROCKY,
               when=[c("is_moon", equals=True), c("radius_km", minimum=1500, maximum=4000), c("insolation", maximum=0.02),
                     drawn("d", "d5")],
               weight="window(ts_k, 57, 125, 5)", state="derived",
               correction="Decided (Boss, 2026-10-18, fork 7): an output of the retention physics; the gate carries it."),
    "j4": dict(kind="option", group="j", classes=ROCKY, when=[drawn("d", "d4", "d5")],
               weight="window(ts_k, 30, 67, 5)", state="derived",
               correction="Decided (Boss, 2026-10-18, fork 7): an output of the retention physics; the gate carries it."),
    "j5": dict(state="stand-in", replacement="GEN.511", kind="option", group="j", classes=ROCKY,
               when=[c("host_type", equals="M"), c("age_gy", minimum=1), c("mass_earth", minimum=0.3, maximum=5),
                     c("zone", equals="e")],
               weight="c_j5"),

    # ---- k. comets ----
    "k1": dict(state="partly verified", limit="13 of 25 comets (HP22 Table 2); the class depends on the heliocentric distance (every coma beyond about 3.5 AU is CO-dominated) and 9 of the 13 rest on an upper limit on CO, so the share is inflated by an amount that cannot be corrected (GEN.520)", kind="option", group="k", classes=["comet", "asteroid"], tunable=True, weight=
               "13/25 if (body_class == 'comet' and perihelion_au < 3) or (body_class == 'asteroid' and family == 'icy') else 0"),
    "k2": dict(state="partly verified", limit="9 of 25 comets (HP22 Table 2); the heliocentric distance and the CO admission rule bias the share (GEN.520)", kind="option", group="k", classes=["comet", "asteroid"], tunable=True, weight=
               "9/25 if (body_class == 'comet' and perihelion_au < 3) or (body_class == 'asteroid' and family == 'icy') else 0"),
    "k3": dict(state="partly verified", limit="3 of 25 comets (HP22 Table 2); the heliocentric distance and the CO admission rule bias the share (GEN.520)", kind="option", group="k", classes=["comet", "asteroid"], tunable=True, weight=
               "3/25 if (body_class == 'comet' and perihelion_au < 3) or (body_class == 'asteroid' and family == 'icy') else 0"),
    "k4": dict(state="partly verified", limit="geometric mean 10.9 percent and 0.51 dex of 21 comets at r_h up to 2.5 AU (Robinson database); one common offset per system is not drawn (N=2 interstellar objects) (GEN.520)", kind="param", classes=["comet", "asteroid"], bins=[0.03, 0.1, 0.3], tunable=True,
               dist={"kind": "normal", "log10": True, "mu": "log10(0.109)", "sigma": 0.51}),
    "k5": dict(state="partly verified", limit="geometric mean 2.5 percent and 0.65 dex of 16 comets at r_h up to 2.5 AU (Robinson database); one common offset per system is not drawn (N=2 interstellar objects) (GEN.520)", kind="param", classes=["comet", "asteroid"], tunable=True,
               dist={"kind": "normal", "log10": True, "mu": "log10(0.025)", "sigma": 0.65}),
    "k6": dict(kind="param", classes=["comet", "asteroid"],
               dists={"CH3OH": {"kind": "loguniform", "lo": 0.002, "hi": 0.026},
                      "NH3": {"kind": "loguniform", "lo": 0.001, "hi": 0.0073},
                      "HCN": {"kind": "loguniform", "lo": 0.0009, "hi": 0.0029}}),
    "k7": dict(kind="param", classes=["comet", "asteroid"], state="stand-in", replacement="GEN.1021", tunable=True,
               limit="CH4, C2H6, H2CO and H2S take the geometric mean and spread of 9 to 18 comets at r_h up to 2.5 AU (Robinson database); N2, Ar, OCS and SO2 are the host's own nitrogen, argon "
                     "and sulfur per carbon over the nebula's carbon in CO, times a trapping share and the median CO/H2O, and O2, Kr and Xe keep one value (physics.comet_species, stand-in "
                     "rows CK1 to CK6): the trapping shares and the three values are tuned to one measured comet (tier T7, tweak L3), so one comet still fixes the seven (GEN.520, GEN.1005)",
               dists={**{name: {"kind": "normal", "log10": True, "mu": f"log10({v})", "sigma": sd}
                         for name, v, sd in (("CH4", 6.8e-3, 0.33), ("C2H6", 4.1e-3, 0.33), ("H2CO", 2.5e-3, 0.35), ("H2S", 7.1e-3, 0.39))},
                      **{name: {"kind": "normal", "log10": True,
                                "mu": f"log10(comet_species_median('{name}', n_over_c, ar_over_c, s_over_c))", "sigma": "comet_species_sigma()"}
                         for name in ("N2", "O2", "OCS", "SO2", "Ar", "Kr", "Xe")}},
               correction="GEN.1005 (2026-10-30): N2, Ar, OCS and SO2 follow the host's ratios to carbon and the width is 0.4 dex (Bar-Nun's 5 K in the formation temperature) in place of the design 0.6; "
                          "O2 presence cannot be set from one comet. Ne is only an upper limit (below 5e-8) and is not drawn."),

    # ---- l. white-dwarf and pulsar planets ----
    "l1": dict(state="stand-in", replacement="GEN.513", kind="option", group="l", classes=ROCKY, when=[c("host_type", one_of=["WD", "PSR"])],
               weight="phi((-34.7*exp(-0.00112*ts_k) + 6)/0.5)",
               values={"p_eq_bar": "10**(-34.7*exp(-0.00112*ts_k))"}),
    "l2": dict(state="stand-in", replacement="GEN.513", kind="option", group="l", classes=ROCKY, when=[c("host_type", one_of=["WD", "PSR"])],
               weight="1 - w('l1')"),
    "l3": dict(state="stand-in", replacement="GEN.513", kind="option", group="l", classes=GIANT, when=[c("host_type", equals="WD")], weight="c_l3"),

    # ---- m. trace gases ----
    "m1": dict(kind="binary", classes=AIR, when=AIR_WHEN + [c("o2_pal", minimum=1e-3)], weight=1,
               dist={"kind": "loguniform", "lo": 192, "hi": 279}),
    "m2": dict(state="stand-in", replacement="GEN.524", kind="binary", classes=AIR, when=AIR_WHEN + [c("o2_pal", maximum=1e-3)], constants={"c_m2": [0.1, 0.6, "log"]},
               weight="c_m2 if (co2_rich and h2_poor and host_type in ('M', 'K')) else 0"),
    "m3": dict(state="partly verified", limit="design flag (decision D4): no exoplanet measures N2O, so the bracket is a labelled design range, never tuned (GEN.525)", kind="binary", classes=AIR, when=AIR_WHEN, constants={"c_m3": [0.05, 0.5, "log"]},
               weight="c_m3 if (life_chemical or sep_reducing) else 0"),
    "m4": dict(state="partly verified", limit="design flag (decision D4): no exoplanet measures PH3, and the Venus claim is contested, so the bracket is a labelled design range, never tuned (GEN.525)", kind="binary", classes=AIR, when=AIR_WHEN, constants={"c_m4": [0.0, 0.05, "uniform"]},
               weight="c_m4 if life_chemical else 0"),
    "m5": dict(state="retired", kind="binary", classes=ROCKY),         # GEN.524 part 2: HCl comes from the halogen budget (see RETIRED)
    "m6": dict(state="retired", kind="binary", classes=ROCKY),         # the same for HF
    "m7": dict(kind="binary", classes=GIANT + BD + ROCKY, weight="phi((t_layer_k - 962)/33)"),
    "m8": dict(kind="binary", classes=GIANT + BD + ROCKY, weight="phi((t_layer_k - 775)/26)"),
    "m9": dict(kind="derived", classes=AIR + GIANT + BD, when=[c("gas_has_co_h2s", equals=True)],
               values={"log10_k_eq": "-1.477 + 4/t_layer_k", "k_eq": "10**(-1.477 + 4/t_layer_k)"}),
    "m10": dict(kind="binary", classes=AIR, when=AIR_WHEN + [c("has_solid_ground", equals=True)], weight=1),
    "m11": dict(state="partly verified", limit="design flag (decision D4): no exoplanet measures NO2, so the bracket is a labelled design range, never tuned (GEN.525)", kind="binary", classes=AIR, when=AIR_WHEN, constants={"c_m11": [0.0, 0.3, "uniform"]},
                weight="c_m11 if (life_chemical or lightning) else 0"),

    # ---- n. clouds, hazes and dust ----
    "n1": dict(kind="windows", classes=GIANT + BD + ROCKY,
               windows=[["Al-Ti oxides", 1850, 2000], ["Fe", 1690, 1840], ["forsterite", 1597, 1698],
                        ["enstatite", 1513, 1597], ["Cr", 1416, 1521], ["MnS", 1271, 1343]]),
    "n2": dict(kind="windows", classes=GIANT + BD + ROCKY, windows=[["Na2S", 929, 996], ["KCl/ZnS", 749, 801]]),
    "n3": dict(kind="windows", classes=GIANT + BD, windows=[["NH4SH", 189, 223], ["H2S ice", 105, 128]]),
    "n4": dict(kind="windows", classes=ANY, when=[c("h2o_present", equals=True)], windows=[["water cloud or ice", 161, 253]]),
    "n5": dict(kind="windows", classes=GIANT + ROCKY + ["comet"],
               windows=[["CO2 ice", 114, 195], ["NH3 ice", 135, 179], ["CH4 ice", 53, 90], ["N2 ice", 31, 62],
                        ["CO ice", 34, 66]]),
    "n6": dict(state="partly verified", limit="A_H against Teq is measured (Fu 2017, N=34); the crossover Teq is a model output, not a measurement (GEN.519)", kind="binary", classes=GIANT + BD + ROCKY, weight="phi((teq_k - 950)/100)"),
    "n7": dict(state="stand-in", replacement="GEN.519", kind="binary", classes=AIR + GIANT + BD, when=[c("carbon_bearing", equals=True)],
               constants={"y": [0.01, 1.0, "log"]}, weight="(1 - w('n6'))*y"),
    "n8": dict(kind="binary", classes=AIR, when=AIR_WHEN, weight="1 if u_star > u_star_threshold else 0"),
    "n9": dict(kind="binary", classes=AIR, when=AIR_WHEN, weight="1 if sulfur_flux_earth >= (1 if h2_air else 100) else 0"),

    # ---- o. ground radioisotopes ----
    "o1": dict(kind="param", classes=ROCKY, when=[c("thin_disc_fgk", equals=True)], bins=[0.7, 1.4],
               dist={"kind": "normal", "log10": True, "mu": 0.08, "sigma": 0.058}),
    "o2": dict(state="partly verified", limit="thick-disc [Th/Mg] from Hypatia (N=32, moved onto the solar-twin scale of o1 by the Hypatia thick-thin offset of -0.019 dex, GEN.664); halo, young and M hosts have no Th data (GEN.514)", kind="param", classes=ROCKY, when=[c("thin_disc_fgk", equals=False)], tunable=True,
               dist={"kind": "normal", "log10": True, "mu": 0.061, "sigma": 0.19, "range": [-0.56, 0.60]},
               source="Hypatia (open API, retrieved 2026-10): [Th/Mg] of thick-disc stars mean -0.063, sd 0.193 (N=32); thin disc -0.044, sd 0.153 "
                      "(N=215): the offset inside the one catalogue is -0.019 +/- 0.036 dex, put on the solar-twin scale of o1 (+0.08, Botelho 2019) "
                      "as 0.08 - 0.019 = 0.061 (DATA, COMPUTED; tunable-gaps 3.2; docs/design/composition/radioisotope-tuning.md section 1)",
               correction="Retuned 2026-10-23 (GEN.514): the log-uniform bracket 59 to 251 percent was the range of 14 solar twins; "
                          "the 32 thick-disc stars give a centre of 86 percent of the solar ratio and a spread of 0.19 dex, the widest "
                          "of the intrinsic 0.056 (Botelho 2019) and the catalogue's 0.15 to 0.19."),
    "o3": dict(state="stand-in", replacement="GEN.514", kind="param", classes=ROCKY, dist={"kind": "loguniform", "lo": 8.9e3, "hi": 1.4e4}),
    "o4": dict(state="partly verified", limit="the production ratio 1.667 to 2.0 and the decay constants are published; the enrichment factor F comes from the disc's star formation history and a 1/tau merger delay law whose 30 Myr floor is a RECALLED stand-in; the gas before the disc began is left out, so F is a lower bound for old stars (GEN.663)", kind="param", classes=ROCKY,
               dist={"kind": "uniform", "lo": "0.97483*1.667*th_u_f*exp(0.105802*age_gy)", "hi": "0.97483*2.0*th_u_f*exp(0.105802*age_gy)"},
               source="Th/U (mass) = 0.97483 PR F(T) exp[(l238 - l232) age], PR = 1.667 to 2.0 (Frebel 2007 Table 1, READ), "
                      "F = th_u_f, the thorium-232 over uranium-238 left of the disc's r-process input at the star's birth: the disc's star formation "
                      "(galaxy.disk_history) convolved with a 1/tau merger delay law, each ejecta decayed with half-lives 14.05 and 4.468 Gyr "
                      "(physics/rprocess_enrichment.py, GEN.663; COMPUTED, tunable-gaps 2.6). A constant input of T Gyr gives F(3) = 1.16254 and F(13) = 1.71724",
               correction="Replaced 2026-10-23 (GEN.514): U[3.8, 4.5] was the Solar System value at 4.5 Gyr and did not grow with the planet's "
                          "age (x1.17 over 1.5 Gyr). Replaced again 2026-10-24 (GEN.663): the birth epoch T drawn between 3 and 13 Gyr is gone, F is the quantity "
                          "`th_u_f` that the star's galaxy supplies."),
    "o5": dict(kind="derived", classes=ROCKY,
               values={"heating_ratio": "0.384*exp(0.693147/4.468*(4.5 - age_gy)) + 0.016*exp(0.693147/0.7038*(4.5 - age_gy))"
                                        " + 0.402*exp(0.693147/14.05*(4.5 - age_gy)) + 0.193*exp(0.693147/1.248*(4.5 - age_gy))"}),
    "o6": dict(state="partly verified", limit="the enriched share spans 0.04 to 0.66 (clusters, white dwarf pollution); no single distribution, so the bracket is the widest and no target exists (GEN.514)", kind="binary", classes=ROCKY, when=[c("planetesimals_early", equals=True)],
               constants={"c_o6": [0.04, 0.66, "log"]}, weight="c_o6"),

    # ---- p. host-star UV regime ----
    "p1": dict(kind="option", group="p", classes=ANY, when=[c("host_type", one_of=["F", "G"]), c("age_gy", minimum=0.1)],
               weight=1),
    "p2": dict(kind="option", group="p", classes=ANY, when=[c("host_type", equals="K"), c("age_gy", minimum=0.1)], weight=1),
    "p3": dict(kind="option", group="p", classes=ANY, when=[c("host_type", equals="M"), c("age_gy", minimum=0.1)],
               weight="1 - phi((log10(loglog(host_mass_solar, 0.55, 0.162, 0.25, 0.6, 0.15, 2.3)) - log10(age_gy))/0.25)"),
    "p4": dict(kind="option", group="p", classes=ANY, when=[c("host_type", equals="M")],
               weight="phi((log10(loglog(host_mass_solar, 0.55, 0.162, 0.25, 0.6, 0.15, 2.3)) - log10(age_gy))/0.25)",
               values={"tau_sat_gy": "loglog(host_mass_solar, 0.55, 0.162, 0.25, 0.6, 0.15, 2.3)",
                       "lx_over_lbol": "10**-3.14*(age_gy/loglog(host_mass_solar, 0.55, 0.162, 0.25, 0.6, 0.15, 2.3))**-1.18"},
               source="saturation 10^-3.14 and time 0.162 Gyr for early M (Jackson 2012 Table 2, READ); 0.6 +/- 0.3 Gyr "
                      "for 0.2 to 0.3 and 2.3 +/- 1.3 Gyr for 0.1 to 0.2 solar masses (Medina 2022, READ); the width of "
                      "0.25 dex on log tau is DESIGN",
               correction="Corrected 2026-10-18 (tunable-gaps 3.3): tau_sat depends on stellar mass (0.16 Gyr at 0.5 to 0.6 "
                          "solar masses, 0.6 at 0.2 to 0.3, 2.3 at 0.1 to 0.2), interpolated in log mass, and replaces "
                          "LU[0.162, 2]."),
    "p5": dict(kind="option", group="p", classes=ANY, when=[c("age_gy", maximum=0.1)], weight=1),
    "p6": dict(state="partly verified", limit="the percentile of an observed rotation distribution is uniform by construction (drawn 0.02 to 0.98, the range Johnstone 2021 offers); the saturation time against mass is read off Johnstone 2021 Fig. 10 (`saturation_time.py`, registry RS1, GEN.493) and is held to 0.1 to 1.2 solar masses; the generator's own star draws its percentile in `activity.generate_activity`, and this row does not yet read that star's value", kind="param", classes=ANY, when=[drawn("p", "p4", "p5")],
               dist={"kind": "uniform", "lo": 0.02, "hi": 0.98}, sets="rotation_percentile",
               values={"t_sat_myr": "saturation_time_myr(host_mass_solar, 100*value)"},
               source="Johnstone, Bartel and Guedel 2021 (arXiv:2009.07695, READ): the saturation age against mass for the 5th, 50th and 95th "
                      "percentiles of the 150 Myr rotation distribution of the clusters (Fig. 10, read off by pixel position), which is the "
                      "observed distribution, so the percentile is uniform; at one solar mass the curves agree with Tu et al. 2015 "
                      "(arXiv:1504.04546, READ: 5.7, 23 and 226 Myr at the 10th, 50th and 90th percentiles)"),
}

# The group a row belongs to carries its fallback; these rows sit in no group (independent draws).

# The constants that are not row constants but part of the physics the expressions use.
CONSTANTS = {
    "pi": 3.141592653589793,
    "G": 6.6743e-11,
    "K_B": 1.380649e-23,
    "AMU_KG": 1.66053906660e-27,
    "R_EARTH_M": 6.371e6,
    "M_EARTH_KG": 5.972e24,
}

# The research items that replace the DESIGN rows (the part's section 4, R1 to R15 = GEN.480 to GEN.494).
RESEARCH = {f"R{n}": f"GEN.{479 + n}" for n in range(1, 16)}
