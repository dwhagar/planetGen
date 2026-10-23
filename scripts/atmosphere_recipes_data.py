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

ROCKY = list("ABCDEFGHKLMNOPQSV")
AIR = list("ABEFGHKLMNOPQV")
GIANT = list("IJT")
HAB = list("EFGHKLMOPQV")
LAVA = list("ABCS")
BD = ["brown-dwarf"]
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


# The tables that compete for 100 percent: a group is drawn once, from the rows of it that pass the subject.
# `fallback` is the row that wins when every weight is zero ("none" is the explicit outcome that no row of the group
# applies). `constants` are the design constants the group's rows draw once inside a bracket, [low, high, shape].
GROUPS = {
    "a": {"rows": ["a5", "a6", "a7"], "fallback": None},
    "b": {"rows": ["b1", "b2", "b3", "b4"], "fallback": None,
          "constants": {"c_b4": [0.05, 0.4, "log"]}},
    "c": {"rows": ["c1", "c2", "c3", "c4"], "fallback": None},
    "c-after": {"rows": ["c5", "c6"], "fallback": "c6"},
    "d": {"rows": ["d1", "d2", "d3", "d4", "d5"], "fallback": "d1"},
    "e": {"rows": ["e1", "e2", "e3", "e4"], "fallback": "e1",
          "constants": {"f_env0": [0.005, 0.10, "log"], "eta": [0.1, 0.5, "log"]}},
    "f": {"rows": ["f1", "f2", "f3", "f4"], "fallback": None},
    "g": {"rows": ["g1", "g2", "g3", "g4", "g5"], "fallback": None},
    "h": {"rows": ["h1", "h2", "h3"], "fallback": None},
    "i": {"rows": ["i1", "i2", "i3", "i4", "i5"], "fallback": None},
    "j": {"rows": ["j1", "j2", "j3", "j4", "j5"], "fallback": "none",
          "constants": {"c_j1": [0.1, 1.0, "log"], "c_j5": [0.05, 0.5, "log"]}},
    "k": {"rows": ["k1", "k2", "k3"], "fallback": "none"},
    "l": {"rows": ["l1", "l2", "l3"], "fallback": "l2", "constants": {"c_l3": [0.0, 0.1, "uniform"]}},
    "p": {"rows": ["p1", "p2", "p3", "p4", "p5"], "fallback": None},
}

# The Berta-Thompson, Wachiraphan and Murray 2025 shoreline (arXiv:2507.02136), posteriors of the main sample, read in
# docs/design/composition/tunable-gaps-sources/group2-delivery-retention.md: log10 f0 2.62, p 5.9, q 1.17, w 0.247 dex.
_DELTA = "log10(insolation) - (2.62 + 5.9*log10(v_esc_kms/11.186) + 1.17*log10(luminosity_solar))"
_P_ATM = "1/(1 + 10**(delta/0.247))"
_SHORE = {"delta": _DELTA, "p_atm": _P_ATM}

MACHINE = {
    # ---- a. mantle redox ----
    "a1": dict(state="stand-in", replacement="GEN.513", kind="param", classes=ROCKY, when=[c("mass_earth", minimum=0.01)], sets="redox_diw",
               dist={"kind": "normal", "mu": "clamp(-2 + 1.25*(log10(mass_earth) + 2), -2, 0.5)", "sigma": 0.5}),
    "a2": dict(kind="param", classes=AIR, when=[c("long_lived_geology", equals=True)], sets="redox_diw_later",
               state="stand-in", replacement="GEN.513",
               dist={"kind": "normal", "mu": "3.5 + 2*log10(mass_earth)", "sigma": 1.0, "range": [-6, 6]},
               correction="A law fitted to the Moon, Mars and Earth (tunable-gaps 3.2): a flagged stand-in until GEN.482."),
    "a3": dict(state="partly verified", limit="22 polluted white dwarf values (READ, COMPUTED), a selection-biased sample of core-formation redox; the scale differs from a1; the width 0.22 is tuned with the detection model (GEN.513)", kind="param", classes=["asteroid", "C", "D"], sets="redox_diw", tunable=True,
               dist={"kind": "normal", "mu": -1.25, "sigma": 0.22},
               source="22 polluted white dwarf values: median -1.26, mean -1.25, sample sd 0.45 of the reported medians (Doyle 2020 Table "
                      "1, Doyle 2023 Table 4, via docs/design/composition/tunable-gaps-sources/group3-redox-wd-radio.md; READ, COMPUTED); "
                      "the width 0.22 is the intrinsic spread of the 13 well-constrained stars (COMPUTED) and lies inside the "
                      "0 to 0.38 the detection model accepts (GEN.513, docs/design/wd-redox-tuning.md); the meteoritic bracket "
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
    "b1": dict(state="stand-in", replacement="GEN.510", kind="option", group="b", classes=ROCKY, weight="(1 - s)**2",
               lets={"s": "1/(1 + exp(-log(snow_ratio)/0.4))"}),
    "b2": dict(state="stand-in", replacement="GEN.510", kind="option", group="b", classes=ROCKY, weight="2*s*(1 - s)",
               lets={"s": "1/(1 + exp(-log(snow_ratio)/0.4))"}),
    "b3": dict(state="stand-in", replacement="GEN.510", kind="option", group="b", classes=ROCKY, weight="s**2",
               lets={"s": "1/(1 + exp(-log(snow_ratio)/0.4))"}),
    "b4": dict(state="stand-in", replacement="GEN.510", kind="option", group="b", classes=ROCKY + ["I"], when=[c("snow_ratio", minimum=1.0)],
               weight="c_b4*s**3", lets={"s": "1/(1 + exp(-log(snow_ratio)/0.4))"}),
    "b5": dict(state="stand-in", replacement="GEN.510", kind="param", classes=ROCKY, when=[drawn("b", "b2", "b3")],
               dist={"kind": "loguniform", "lo": 5e-4, "hi": 5e-3}),
    "b6": dict(state="stand-in", replacement="GEN.510", kind="binary", classes=["C", "D", "S", "I"], when=[c("zone", equals="r")],
               weight="0.8 if body_class == 'I' else 0.35",
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
    "d7": dict(kind="flag", classes=["S", "V", "I", "T"], when=[c("mass_earth", minimum=1, maximum=30)],
               weight="1 if 2.016*AMU_KG*(v_esc_kms*1000)**2/(2*K_B*(teq_k + min(1500, 750*sqrt(xuv_index)))) >= 30 else 0"),

    # ---- e. sub-Neptune and super-Earth envelope ----
    "e1": dict(kind="option", group="e", classes=["S", "V", "I", "T"], when=[c("mass_earth", minimum=1, maximum=30)],
               weight="(1 - p('b4'))*min(1, f_lost/f_env0)", lets={"f_lost": "0.1*pi*(core_radius_earth*R_EARTH_M)**3*3.1e16*insolation/(G*(mass_earth*M_EARTH_KG)**2)"}),
    "e2": dict(kind="option", group="e", classes=["S", "V", "I", "T"], when=[c("mass_earth", minimum=1, maximum=30)],
               weight="(1 - p('b4'))*(1 - min(1, f_lost/f_env0))",
               lets={"f_lost": "0.1*pi*(core_radius_earth*R_EARTH_M)**3*3.1e16*insolation/(G*(mass_earth*M_EARTH_KG)**2)"}),
    "e3": dict(state="partly verified", limit="the gate is Madhusudhan 2021 (READ); the share is the b4 stand-in (GEN.510)", kind="option", group="e", classes=["S", "V", "I", "T"], when=[c("mass_earth", minimum=1, maximum=30)],
               weight="p('b4')*(1 - eta)"),
    "e4": dict(state="partly verified", limit="the gate is Madhusudhan 2021 (READ); the share is the b4 stand-in (GEN.510)", kind="option", group="e", classes=["S", "V", "I", "T"], when=[c("mass_earth", minimum=1, maximum=30)],
               weight="p('b4')*eta"),
    "e5": dict(state="partly verified", limit="bracket edges from about 8 JWST retrievals and one 15-planet fit; no population mean exists (GEN.509)", kind="param", classes=["S", "V", "I", "T"], when=[drawn("e", "e2")],
               dist={"kind": "loguniform", "lo": 0.01, "hi": 0.6}, bins=[0.01, 0.03, 0.1, 0.3, 0.6],
               values={"mean_molecular_weight": "2.4 + 3.1*log10(value/0.01)/log10(58)"}),
    "e6": dict(kind="binary", classes=["S", "V", "I", "T"], when=[drawn("e", "e2")], weight="phi((500 - teq_k)/100)"),

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
    "g1": dict(state="partly verified", limit="late M is the extension of the map above L0 (2238 K); the map has no M-type fit (Table 13) and the 2800 K cap is the repo's rogue maximum (GEN.517)", kind="option", group="g", classes=BD, weight="phi((2800 - teff_k)/100) - phi((2238 - teff_k)/134)"),
    "g2": dict(kind="option", group="g", classes=BD, weight="phi((2238 - teff_k)/134) - phi((1273 - teff_k)/79)"),
    "g3": dict(kind="option", group="g", classes=BD, weight="phi((1273 - teff_k)/79) - phi((686 - teff_k)/79)"),
    "g4": dict(state="partly verified", limit="late T (T8 and later) is a subclass of T, not a type of its own; its edge at T8 (686 K) is a design choice on the map (GEN.517)", kind="option", group="g", classes=BD, weight="phi((686 - teff_k)/79) - phi((460 - teff_k)/79)"),
    "g5": dict(kind="option", group="g", classes=BD, weight="phi((460 - teff_k)/79)"),
    "g6": dict(kind="binary", classes=BD, weight="phi((1100 - teff_k)/100)"),

    # ---- h. hot Jupiters ----
    "h1": dict(kind="option", group="h", classes=GIANT, when=[c("zone", equals="h")],
               weight="window(tday_k, None, 1500, 100)"),
    "h2": dict(state="stand-in", replacement="GEN.518", kind="option", group="h", classes=GIANT, when=[c("zone", equals="h")],
               weight="window(tday_k, 1500, 2200, 100)"),
    "h3": dict(state="stand-in", replacement="GEN.518", kind="option", group="h", classes=GIANT, when=[c("zone", equals="h")],
               weight="window(tday_k, 2200, None, 100)"),
    "h4": dict(state="partly verified", limit="a Ti threshold between 2350 and 2640 K from N=3 to 6 planets; no population measures Kzz (GEN.519)", kind="binary", classes=GIANT, when=[c("zone", equals="h"), drawn("h", "h2", "h3")],
               weight="clamp((11 - (7.48 + 3.52*(log10(grain_um) + 1)/2))/4, 0, 1)"),
    "h5": dict(kind="binary", classes=GIANT, when=[c("zone", equals="h")], weight="pdist('f7', 1.0, None)"),

    # ---- i. lava worlds ----
    "i1": dict(kind="option", group="i", classes=LAVA, when=[c("zone", equals="h")], weight="window(tp_k, None, 2033, 60)"),
    "i2": dict(kind="option", group="i", classes=LAVA, when=[c("zone", equals="h")], weight="window(tp_k, 2033, 2588, 60)"),
    "i3": dict(kind="option", group="i", classes=LAVA, when=[c("zone", equals="h")], weight="window(tp_k, 2588, 2890, 60)"),
    "i4": dict(kind="option", group="i", classes=LAVA, when=[c("zone", equals="h")], weight="window(tp_k, 2890, 3168, 60)"),
    "i5": dict(kind="option", group="i", classes=LAVA, when=[c("zone", equals="h")], weight="window(tp_k, 3168, None, 60)"),
    "i6": dict(state="partly verified", limit="2 of 2 molten-dayside ultra-short-period planets keep volatiles (N=2, Beta interval 0.29 to 0.99; GEN.511)", kind="binary", classes=LAVA, when=[c("zone", equals="h"), c("mass_earth", minimum=1, maximum=10)],
               constants={"c_i6": [0.1, 0.9, "log"]}, weight="c_i6",
               values={"p_eq_bar": "10**(-34.7*exp(-0.00112*tp_k))"}),

    # ---- j. frost-equilibrium and evolved atmospheres ----
    "j1": dict(state="partly verified", limit="the runaway gate is Kopparapu 2013 (READ); the retained fraction is a design residual in [0, 1] (GEN.511)", kind="option", group="j", classes=["N", "V", "A", "S"],
               when=[c("redox_diw", minimum=1), c("water_lost", equals=True), c("runaway_ratio", minimum=1.0)],
               weight="c_j1"),
    "j2": dict(kind="option", group="j", classes=["K", "L", "U"],
               when=[c("mass_earth", maximum=0.5), drawn("d", "d3"), c("ts_k", minimum=190, maximum=270)], weight=1,
               state="derived",
               correction="Decided (Boss, 2026-10-18, fork 7): the Mars-, Titan- and Pluto-like occurrence is an output of "
                          "the escape and retention physics, not a design flag; the gate carries it and the weight is 1."),
    "j3": dict(kind="option", group="j", classes=["Y"],
               when=[c("is_moon", equals=True), c("radius_km", minimum=1500, maximum=4000), c("insolation", maximum=0.02),
                     drawn("d", "d5")],
               weight="window(ts_k, 57, 125, 5)", state="derived",
               correction="Decided (Boss, 2026-10-18, fork 7): an output of the retention physics; the gate carries it."),
    "j4": dict(kind="option", group="j", classes=["C", "D", "U"], when=[drawn("d", "d4", "d5")],
               weight="window(ts_k, 30, 67, 5)", state="derived",
               correction="Decided (Boss, 2026-10-18, fork 7): an output of the retention physics; the gate carries it."),
    "j5": dict(state="stand-in", replacement="GEN.511", kind="option", group="j", classes=["Z"],
               when=[c("host_type", equals="M"), c("age_gy", minimum=1), c("mass_earth", minimum=0.3, maximum=5),
                     c("zone", equals="e")],
               weight="c_j5"),

    # ---- k. comets ----
    "k1": dict(state="stand-in", replacement="GEN.520", kind="option", group="k", classes=["comet", "asteroid"], weight=
               "0.5 if (body_class == 'comet' and perihelion_au < 3) or (body_class == 'asteroid' and family == 'icy') else 0"),
    "k2": dict(state="stand-in", replacement="GEN.520", kind="option", group="k", classes=["comet", "asteroid"], weight=
               "1/3 if (body_class == 'comet' and perihelion_au < 3) or (body_class == 'asteroid' and family == 'icy') else 0"),
    "k3": dict(state="stand-in", replacement="GEN.520", kind="option", group="k", classes=["comet", "asteroid"], weight=
               "0.1 if (body_class == 'comet' and perihelion_au < 3) or (body_class == 'asteroid' and family == 'icy') else 0"),
    "k4": dict(state="stand-in", replacement="GEN.520", kind="param", classes=["comet", "asteroid"], bins=[0.03, 0.1, 0.3],
               dist={"kind": "normal", "log10": True, "mu": "log10(0.12)", "sigma": 0.5}),
    "k5": dict(state="stand-in", replacement="GEN.520", kind="param", classes=["comet", "asteroid"],
               dist={"kind": "normal", "log10": True, "mu": "log10(0.03)", "sigma": 0.6, "range": ["log10(0.005)", "log10(0.23)"]}),
    "k6": dict(kind="param", classes=["comet", "asteroid"],
               dists={"CH3OH": {"kind": "loguniform", "lo": 0.002, "hi": 0.026},
                      "NH3": {"kind": "loguniform", "lo": 0.001, "hi": 0.0073},
                      "HCN": {"kind": "loguniform", "lo": 0.0009, "hi": 0.0029}}),
    "k7": dict(kind="param", classes=["comet", "asteroid"], state="stand-in", replacement="GEN.520",
               dists={name: {"kind": "normal", "log10": True, "mu": f"log10({v})", "sigma": 0.6}
                      for name, v in (("CH4", 3.4e-3), ("C2H6", 2.9e-3), ("H2CO", 3.2e-3), ("H2S", 1.1e-2), ("N2", 8.9e-4),
                                      ("O2", 3.1e-2), ("OCS", 4.1e-4), ("SO2", 1.3e-3), ("Ar", 5.8e-6), ("Kr", 4.9e-7),
                                      ("Xe", 2.4e-7))},
               correction="Built from one measured comet (67P, Rubin 2019): a flagged stand-in until GEN.488. Ne is only an "
                          "upper limit (below 5e-8) and is not drawn."),

    # ---- l. white-dwarf and pulsar planets ----
    "l1": dict(state="stand-in", replacement="GEN.513", kind="option", group="l", classes=ROCKY, when=[c("host_type", one_of=["WD", "PSR"])],
               weight="phi((-34.7*exp(-0.00112*ts_k) + 6)/0.5)",
               values={"p_eq_bar": "10**(-34.7*exp(-0.00112*ts_k))"}),
    "l2": dict(state="stand-in", replacement="GEN.513", kind="option", group="l", classes=ROCKY, when=[c("host_type", one_of=["WD", "PSR"])],
               weight="1 - w('l1')"),
    "l3": dict(state="stand-in", replacement="GEN.513", kind="option", group="l", classes=GIANT, when=[c("host_type", equals="WD")], weight="c_l3"),

    # ---- m. trace gases ----
    "m1": dict(kind="binary", classes=AIR, when=[c("o2_pal", minimum=1e-3)], weight=1,
               dist={"kind": "loguniform", "lo": 192, "hi": 279}),
    "m2": dict(state="stand-in", replacement="GEN.524", kind="binary", classes=AIR, when=[c("o2_pal", maximum=1e-3)], constants={"c_m2": [0.1, 0.6, "log"]},
               weight="c_m2 if (co2_rich and h2_poor and host_type in ('M', 'K')) else 0"),
    "m3": dict(state="stand-in", replacement="GEN.525", kind="binary", classes=AIR, constants={"c_m3": [0.05, 0.5, "log"]},
               weight="c_m3 if (life_chemical or sep_reducing) else 0"),
    "m4": dict(state="stand-in", replacement="GEN.525", kind="binary", classes=AIR, constants={"c_m4": [0.0, 0.05, "uniform"]}, weight="c_m4"),
    "m5": dict(state="stand-in", replacement="GEN.524", kind="binary", classes=ROCKY, when=[drawn("c", "c1", "c2", "c3", "c4")],
               constants={"c_m5": [0.1, 0.5, "log"]}, weight="c_m5"),
    "m6": dict(state="stand-in", replacement="GEN.524", kind="binary", classes=ROCKY, when=[drawn("c", "c1", "c2", "c3", "c4")],
               constants={"c_m6": [0.1, 0.5, "log"]}, weight="c_m6"),
    "m7": dict(kind="binary", classes=GIANT + BD + LAVA, weight="phi((t_layer_k - 962)/33)"),
    "m8": dict(kind="binary", classes=GIANT + BD + LAVA, weight="phi((t_layer_k - 775)/26)"),
    "m9": dict(kind="derived", classes=AIR + GIANT + BD, when=[c("gas_has_co_h2s", equals=True)],
               values={"log10_k_eq": "-1.477 + 4/t_layer_k", "k_eq": "10**(-1.477 + 4/t_layer_k)"}),
    "m10": dict(kind="binary", classes=AIR, when=[c("has_solid_ground", equals=True)], weight=1),
    "m11": dict(state="stand-in", replacement="GEN.525", kind="binary", classes=AIR, constants={"c_m11": [0.0, 0.3, "uniform"]},
                weight="c_m11 if (life_chemical or lightning) else 0"),

    # ---- n. clouds, hazes and dust ----
    "n1": dict(kind="windows", classes=GIANT + BD + LAVA,
               windows=[["Al-Ti oxides", 1850, 2000], ["Fe", 1690, 1840], ["forsterite", 1597, 1698],
                        ["enstatite", 1513, 1597], ["Cr", 1416, 1521], ["MnS", 1271, 1343]]),
    "n2": dict(kind="windows", classes=GIANT + BD + LAVA, windows=[["Na2S", 929, 996], ["KCl/ZnS", 749, 801]]),
    "n3": dict(kind="windows", classes=GIANT + BD, windows=[["NH4SH", 189, 223], ["H2S ice", 105, 128]]),
    "n4": dict(kind="windows", classes=ANY, when=[c("h2o_present", equals=True)], windows=[["water cloud or ice", 161, 253]]),
    "n5": dict(kind="windows", classes=GIANT + ["K", "Y", "D", "comet"],
               windows=[["CO2 ice", 114, 195], ["NH3 ice", 135, 179], ["CH4 ice", 53, 90], ["N2 ice", 31, 62],
                        ["CO ice", 34, 66]]),
    "n6": dict(state="partly verified", limit="A_H against Teq is measured (Fu 2017, N=34); the crossover Teq is a model output, not a measurement (GEN.519)", kind="binary", classes=GIANT + BD + LAVA, weight="phi((teq_k - 950)/100)"),
    "n7": dict(state="stand-in", replacement="GEN.519", kind="binary", classes=AIR + GIANT + BD, when=[c("carbon_bearing", equals=True)],
               constants={"y": [0.01, 1.0, "log"]}, weight="(1 - w('n6'))*y"),
    "n8": dict(kind="binary", classes=AIR, weight="1 if u_star > u_star_threshold else 0"),
    "n9": dict(kind="binary", classes=AIR, weight="1 if sulfur_flux_earth >= (1 if h2_air else 100) else 0"),

    # ---- o. ground radioisotopes ----
    "o1": dict(kind="param", classes=ROCKY, when=[c("thin_disc_fgk", equals=True)], bins=[0.7, 1.4],
               dist={"kind": "normal", "log10": True, "mu": 0.08, "sigma": 0.058}),
    "o2": dict(state="partly verified", limit="thin and thick disc [Th/Mg] from Hypatia (N=215, 32); halo, young and M hosts have no data (GEN.514)", kind="param", classes=ROCKY, when=[c("thin_disc_fgk", equals=False)],
               dist={"kind": "loguniform", "lo": "10**-0.23", "hi": "10**0.40"}),
    "o3": dict(state="stand-in", replacement="GEN.514", kind="param", classes=ROCKY, dist={"kind": "loguniform", "lo": 8.9e3, "hi": 1.4e4}),
    "o4": dict(state="stand-in", replacement="GEN.514", kind="param", classes=ROCKY, dist={"kind": "uniform", "lo": 3.8, "hi": 4.5}),
    "o5": dict(kind="derived", classes=ROCKY,
               values={"heating_ratio": "0.384*exp(0.693147/4.468*(4.5 - age_gy)) + 0.016*exp(0.693147/0.7038*(4.5 - age_gy))"
                                        " + 0.402*exp(0.693147/14.05*(4.5 - age_gy)) + 0.193*exp(0.693147/1.248*(4.5 - age_gy))"}),
    "o6": dict(state="partly verified", limit="the enriched share spans 0.04 to 0.66 (clusters, white dwarf pollution); no single distribution (GEN.514)", kind="binary", classes=ROCKY, when=[c("planetesimals_early", equals=True)],
               constants={"c_o6": [0.1, 1.0, "log"]}, weight="c_o6"),

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
    "p6": dict(state="partly verified", limit="percentile tracks of Tu 2015 (READ); the 0.25/0.5/0.25 weights are replaced by a uniform percentile (GEN.523)", kind="choice", classes=ANY, when=[drawn("p", "p4", "p5")],
               options={"slow": {"weight": 0.25, "percentile": 10, "t_sat_myr": 5.7},
                        "median": {"weight": 0.5, "percentile": 50, "t_sat_myr": 23},
                        "fast": {"weight": 0.25, "percentile": 90, "t_sat_myr": 226}}),
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
