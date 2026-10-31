#!/usr/bin/env python3
"""
Builds the list of every stand-in value (GEN.725) from the stand-in registry of the gap-filling plan.

    python scripts/build_stand_in_values.py            # write docs/design/research/stand-in-values.md
    python scripts/build_stand_in_values.py --check    # build in memory and fail if the committed file differs

The registry is the markdown tables of `docs/design/research/gap-filling-plan.md`: the ranked fills of section 5 (5.1 to 5.3), the
generated dial ranges of 6.3a, the fills of the shared download step in 6.6 and the stand-ins of the atmosphere mixture draw in 6.7 and those of the young-group draw in 6.8.
Every one of those tables carries a Tier column (Boss's nine source tiers, T1 to T9) and a Tweak column (L0 to L4). The list is
copied from them and never typed by hand, so it cannot drift; `src/tests/test_stand_in_values.py` fails when it does or when a row
lacks a tier or a tweak.
"""

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN = os.path.join(HERE, "..", "docs", "design", "research", "gap-filling-plan.md")
OUT = os.path.join(HERE, "..", "docs", "design", "research", "stand-in-values.md")
WIRING = os.path.join(HERE, "..", "docs", "design", "research", "wiring-register.md")
BURN_DOWN = os.path.join(HERE, "..", "docs", "design", "research", "burn-down-register.md")

TIERS = (
    ("T1", "Verified", "read in the source and the value is the source's"),
    ("T2", "Confident unverified", "from a source we trust but have not opened or cannot check"),
    ("T3", "Marginal unverified", "from a source we doubt, or a bracket from sources that disagree"),
    ("T4", "Unverified", "no source at all for the value, only a reason for its range"),
    ("T5", "Confident educated guess based on physics", "a law of physics fixes the form and the size to within a factor of a few"),
    ("T6", "Marginal educated guess based on physics", "physics fixes the form, but the size is argued, not derived"),
    ("T7", "Educated guess based on inferred physics", "the physics is itself inferred from other bodies or other regimes"),
    ("T8", "Educated guess based on reputable fringe or new physics", "a published but not yet accepted mechanism"),
    ("T9", "Educated guess based on needed output (full tweak)", "chosen so that the generated sky matches what has to come out"),
)
TWEAKS = (
    ("L0", "As envisioned", "the value is used as its source gives it"),
    ("L1", "Nudged inside the source's uncertainty", "moved within the error bar the source prints"),
    ("L2", "Tuned to a published statistic", "moved until a published count or distribution is reproduced (a tuning target)"),
    ("L3", "Hand-tuned against a sweep", "chosen by running the generator over many draws and looking at the outcome"),
    ("L4", "Fitted to the outcome with no independent source", "the number is whatever makes the output right"),
)

TIER_CODE = re.compile(r"T[1-9]\b")
TWEAK_CODE = re.compile(r"L[0-4]\b")

SECTIONS = (
    ("5.1", "5.1 Galaxy laws", ("Id", "Kind", "Tier", "Tweak", "Fill (ranked)", "Wrong if", "Result lands as")),
    ("5.2", "5.2 Body composition", ("Id", "Kind", "Tier", "Tweak", "Fill (ranked)", "Wrong if", "Result lands as")),
    ("5.3", "5.3 Atmosphere", ("Id", "Kind", "Tier", "Tweak", "Fill (ranked)", "Wrong if", "Result lands as")),
    ("6.3a", "6.3a Ranges of the galaxy shape dials", ("Lever", "Range", "Label", "Tier", "Tweak", "Source")),
    ("6.6", "6.6 Fills of the shared download step", ("fill_id", "kind", "tier", "tweak", "stand-in", "replaced_by", "sources")),
    ("6.7", "6.7 Stand-ins of the atmosphere mixture draw",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.8", "6.8 Stand-ins of the hours outside in a day",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.9", "6.9 Stand-ins of how a rendered body is drawn",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.10", "6.10 Stand-ins of the radius valley",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.11", "6.11 Stand-ins of the young-group draw",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.12", "6.12 Stand-ins of the orbit elements",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "replaced_by", "evidence")),
    ("6.13", "6.13 Stand-ins of the tidal heating",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "replaced_by", "evidence")),
    ("6.14", "6.14 Stand-ins of the spin and its despinning",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "replaced_by", "evidence")),
    ("6.15", "6.15 Stand-ins of the trapped radiation belts",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "replaced_by", "evidence")),
    ("6.16", "6.16 Stand-ins of the interior heat and the r-process uranium",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.18", "6.18 Stand-ins of the life constants",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.19", "6.19 Stand-ins of the life class rules",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.20", "6.20 Stand-ins of the planet mass law",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.21", "6.21 Stand-ins of the life timeline and the biogenic oxygen",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.22", "6.22 Stand-ins of the class cut-over",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.23", "6.23 Stand-ins of the tidal rheology",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.24", "6.24 Stand-ins of the spin torques",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.25", "6.25 Stand-ins of the moon orbits",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.26", "6.26 Stand-ins of the birth spin from the last impacts",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.30", "6.30 Stand-ins of the dynamo field",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.142", "6.142 Stand-ins of the dynamo wiring",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.143", "6.143 Stand-ins of the belt dose a moon takes",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.144", "6.144 Stand-ins of the induced field of cores and planets",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.31", "6.31 Stand-ins of the accreting compact binaries",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "replaced_by", "evidence")),
    ("6.32", "6.32 Stand-ins of the neutron star structure",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.33", "6.33 Stand-ins of the stellar spectra",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.34", "6.34 Stand-ins of the alpha-enhanced tracks and the age-metallicity relation",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.35", "6.35 Stand-ins of the condensates by temperature",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.36", "6.36 Stand-ins of the planet architecture",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.37", "6.37 Stand-ins of the brown dwarf ceiling",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.38", "6.38 Stand-ins of the aluminium-26 and iron-60 at birth",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.39", "6.39 Stand-ins of the X-ray saturation time by mass",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.40", "6.40 Stand-ins of the brackets the atmosphere recipes still draw by design",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.41", "6.41 Stand-ins of the derived taxon kernel",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.42", "6.42 Stand-ins of the brackets the research items GEN.480 and GEN.489 to GEN.493 leave",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.43", "6.43 Stand-ins of the mass left in a belt",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.44", "6.44 Stand-ins of the chondrite-group map",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.45", "6.45 Stand-ins of the rock types from bulk oxides",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.46", "6.46 Stand-ins of the comet ice budget",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.47", "6.47 Stand-ins of the cosmic-ray and particle doses and the ozone spread",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "replaced_by", "evidence")),
    ("6.48", "6.48 Stand-ins of the star counts' composition",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.49", "6.49 Stand-ins of the star class edges and the brown dwarf companion",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.50", "6.50 Stand-ins of the stellar mass function",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.51", "6.51 Stand-ins of the stars after the main sequence",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.52", "6.52 Stand-ins of the stellar [X/Fe] refit to Brewer and Fischer",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.53", "6.53 Stand-ins of the debris disks",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.60", "6.60 Stand-ins of the multiple-star limits and the white-dwarf planets",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.61", "6.61 Stand-ins of the potassium depletion by body class",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.62", "6.62 Stand-ins of the carbon dioxide far wing and the sulfuric acid cloud",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.63", "6.63 Stand-ins of the acid cloud's infrared optics and the minor-gas wings",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.80", "6.80 Stand-ins of the 3+1 quadruple stability network",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.81", "6.81 Stand-ins of the 2+2 quadruple stability network",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.82", "6.82 Stand-ins of the oddity rates confirmed against papers",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.64", "6.64 Stand-ins of the NH3 active-Sun factor",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.56", "6.56 Stand-ins of the registered rows that rest on a Solar System body",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.57", "6.57 Stand-ins of the magma-ocean redox model",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.58", "6.58 Stand-ins of the dry-body carbon",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.65", "6.65 Stand-ins of the arm contrast of each stellar population",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.66", "6.66 Stand-ins of the low-mass main sequence and the star/brown-dwarf boundary",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.67", "6.67 Stand-ins of the stellar wind and flares at a world",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.68", "6.68 Stand-ins of the open-ocean equilibrium tide",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.69", "6.69 Stand-ins of the surface liquids",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.70", "6.70 Stand-ins of the ring origin",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.71", "6.71 Stand-ins of the resonance rings of small bodies",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.72", "6.72 Stand-ins of the oddity table and the weirdness dial",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.73", "6.73 Stand-ins of the disintegrating rocky planet",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.74", "6.74 Stand-ins of atmosphere retention",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.75", "6.75 Stand-ins of the nightside collapse of a locked planet's air",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.86", "6.86 Stand-ins of the rogue planet top-up",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.87", "6.87 Stand-ins of the rogue planet ejection speed",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.88", "6.88 Stand-ins of the rogue planet age",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.89", "6.89 Stand-ins of the rogue planet drift",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.90", "6.90 Stand-ins of the Venus-type spin",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.91", "6.91 Stand-ins of the seasonal swing",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.92", "6.92 Stand-ins of the resonant chains",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.98", "6.98 Stand-ins of the band winds",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.99", "6.99 Stand-ins of the belt density factor",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.76", "6.76 Stand-ins of the super-puff radius",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.77", "6.77 Stand-ins of the helium enrichment of a thinning envelope",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.78", "6.78 Stand-ins of the giant-impact tilt of an ice giant",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.93", "6.93 Stand-ins of the collisional mantle stripping",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.79", "6.79 Stand-ins of the Hycean box",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.83", "6.83 Stand-ins of the scatter step of an outer disc",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.103", "6.103 Stand-ins of the interstellar asteroid fields",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.101", "6.101 Stand-ins of the oddity rates, second batch",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.100", "6.100 Stand-ins of the captured asteroid moons",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.110", "6.110 Stand-ins of the masses of ejected bodies",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.111", "6.111 Stand-ins of the air and water of ejected planets",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.112", "6.112 Stand-ins of the moons carried out with ejected planets",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.113", "6.113 Stand-ins of the planets a birth cluster strips",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.114", "6.114 Stand-ins of the radiogenic heat of a rogue",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.115", "6.115 Stand-ins of the composition of a belt member",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.118", "6.118 Stand-ins of the captured small bodies of a planet",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.119", "6.119 Stand-ins of the composition of captured bodies",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.120", "6.120 Stand-ins of the asteroid belt page",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.121", "6.121 Stand-ins of the Lagrange trojans",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.122", "6.122 Stand-ins of the albedo colour of the star",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.125", "6.125 Stand-ins of the whole bodies at Lagrange points",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.126", "6.126 Stand-ins of the roundness engine and the binding energy",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.134", "6.134 Stand-ins of the structure of a giant",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "replaced_by", "evidence")),
    ("6.135", "6.135 Stand-ins of a giant's internal heat",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "replaced_by", "evidence")),
    ("6.136", "6.136 Stand-ins of the core's thermal history",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.150", "6.150 Stand-ins of the shoulders of the best-fit classes",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.156", "6.156 Stand-ins of the two cuts of the body type",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.162", "6.162 Stand-ins of the temperature a body took its shape at",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.161", "6.161 Stand-ins of the mass of an asteroid field",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.160", "6.160 Stand-ins of the lumpy mesh of a small body",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.159", "6.159 Stand-ins of the star stepper",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.158", "6.158 Stand-ins of the acceleration kernel",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.157", "6.157 Stand-ins of the galactic potential",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.155", "6.155 Stand-ins of the life readings",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.129", "6.129 Stand-ins of the strengths of the ice types",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.140", "6.140 Stand-ins of the gas disc of a satellite galaxy",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.141", "6.141 Stand-ins of the gas of a satellite galaxy by its distance from the host",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.94", "6.94 Stand-ins of the irregular zone",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.102", "6.102 Stand-ins of the smoothed rogue density",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.127", "6.127 Stand-ins of the high-gravity class line",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.128", "6.128 Stand-ins of the class lines of the train 6 epoch",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.105", "6.105 Stand-ins of the field dispersal",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.106", "6.106 Stand-ins of the field rates",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.107", "6.107 Stand-ins of the field composition",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.108", "6.108 Stand-ins of the light top-up field",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.116", "6.116 Stand-ins of the resonance and chaotic-zone map of a belt",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.117", "6.117 Stand-ins of the history of a belt",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.124", "6.124 Stand-ins of the wiring of aluminium-26 at birth",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.130", "6.130 Stand-ins of the decay exponent after X-ray saturation",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.131", "6.131 Stand-ins of the halo mass function and the conditional halo",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.132", "6.132 Stand-ins of the reference black hole draw",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.133", "6.133 Stand-ins of the day-and-night state",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.146", "6.146 Stand-ins of the thick disc's ages",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.147", "6.147 Stand-ins of the chemically peculiar stars",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.148", "6.148 Stand-ins of the brightness of a ring and the dust that darkens it",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.138", "6.138 Stand-ins of the variable-star kinds",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.139", "6.139 Stand-ins of the carbon and S stars",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.145", "6.145 Stand-ins of the cluster maximum mass and cloud-shock disruption of the central ring",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.95", "6.95 Stand-ins of the belt sources",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.96", "6.96 Stand-ins of the irregular orbits",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.97", "6.97 Stand-ins of the belt origins",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.59", "6.59 Stand-ins of the comet species",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
    ("6.127", "6.127 Stand-ins of the architecture's losses",
     ("fill_id", "kind", "tier", "tweak", "constant", "knobs", "validity", "falsifier", "replaced_by", "evidence")),
)
"""tuple: Each registry table: the plan's section number, its title here, the columns that go in the list."""

UNREAD_COLUMN = "Not read yet (item)"
"""str: A column of every table in the list: the item that will wire a row whose constants nothing in the app reads yet
(`wiring-register.md`, GEN.754), or a dash."""

BURN_DOWN_COLUMN = "Burn-down (GEN.946)"
"""str: The last column of every table in the list: what a row at tweak level L3 or L4 needs before it can go (`burn-down-register.md`:
`no data` or `replace` with the item), or a dash."""


def unread_rows(text):
    """
    The registry rows that the wiring register (GEN.754) says nothing in the app reads yet.

    Returns:
        dict: `"{section}:{fill_id}"` to the item that will wire it, or `by design`.
    """
    found = {}
    for line in text.splitlines():
        cells = [cell.strip().strip("`") for cell in line.strip().strip("|").split("|")]
        if line.startswith("|") and len(cells) == 4 and cells[0] == "row":
            found[cells[1]] = cells[2]
    return found


def burn_down_rows(text):
    """
    The burn-down register's lines (GEN.946).

    Returns:
        dict: `"{section}:{fill_id}"` to `(disposition, item, reason)`; the disposition is `no data`, `replace` or `replaced`.
    """
    found = {}
    for line in text.splitlines():
        cells = [cell.strip().strip("`") for cell in line.strip().strip("|").split("|")]
        if line.startswith("|") and len(cells) == 4 and cells[1] in ("no data", "replace", "replaced"):
            found[cells[0]] = (cells[1], cells[2], cells[3])
    return found


def _cells(line):
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def registry_tables(text):
    """
    Every registry table of the plan.

    Returns:
        dict: Section number to `(header, rows)`; `header` the list of column names and `rows` lists of cells. Rows whose cell
            count differs from the header are returned too (the test names them).
    """
    tables, section, lines = {}, None, text.splitlines()
    wanted = {number for number, _, _ in SECTIONS}
    i = 0
    while i < len(lines):
        line = lines[i]
        heading = re.match(r"#{2,4} (\d+\.\d+a?)\b", line)
        if heading:
            section = heading.group(1)
        if section in wanted and section not in tables and line.startswith("|") and i + 1 < len(lines) and lines[i + 1].startswith("|---"):
            header, rows, i = _cells(line), [], i + 2
            while i < len(lines) and lines[i].startswith("|"):
                rows.append(_cells(lines[i]))
                i += 1
            tables[section] = (header, rows)
            continue
        i += 1
    return tables


def problems(tables):
    """The registry rows that lack a tier or a tweak, or that have the wrong number of cells, one sentence each."""
    found = []
    for number, _, _ in SECTIONS:
        if number not in tables:
            found.append(f"section {number} has no registry table")
            continue
        header, rows = tables[number]
        lowered = [name.lower() for name in header]
        if "tier" not in lowered or "tweak" not in lowered:
            found.append(f"section {number} has no Tier and Tweak columns")
            continue
        if not rows:
            found.append(f"section {number} has no rows")
        tier_at, tweak_at = lowered.index("tier"), lowered.index("tweak")
        for row in rows:
            if len(row) != len(header):
                found.append(f"{number} row {row[0]!r}: {len(row)} cells for {len(header)} columns")
                continue
            if not TIER_CODE.match(row[tier_at]) or len(set(re.findall(r"T\d+", row[tier_at])) - {t for t, _, _ in TIERS}):
                found.append(f"{number} row {row[0]!r}: tier {row[tier_at]!r} is not T1 to T9")
            if not TWEAK_CODE.fullmatch(row[tweak_at]):
                found.append(f"{number} row {row[0]!r}: tweak {row[tweak_at]!r} is not L0 to L4")
    return found


def _first_option(text):
    """The first-ranked option of a "Fill (ranked)" cell (`1 ... 2 ...`), or the cell when it is not ranked."""
    match = re.match(r"1 (.*?)(?: 2 [A-Z`(].*)?$", text)
    return match.group(1).strip() if match else text


def _counts(tables):
    """The number of rows at each tier and each tweak, counting the first code of a cell."""
    tier_count = {code: 0 for code, _, _ in TIERS}
    tweak_count = {code: 0 for code, _, _ in TWEAKS}
    cross = {(t, w): 0 for t, _, _ in TIERS for w, _, _ in TWEAKS}
    total = 0
    for number, _, _ in SECTIONS:
        header, rows = tables[number]
        lowered = [name.lower() for name in header]
        for row in rows:
            tier, tweak = TIER_CODE.match(row[lowered.index("tier")]).group(0), row[lowered.index("tweak")]
            tier_count[tier] += 1
            tweak_count[tweak] += 1
            cross[(tier, tweak)] += 1
            total += 1
    return total, tier_count, tweak_count, cross


def render(plan_text, wiring_text="", burn_down_text=""):
    """The text of the list; `wiring_text` is the wiring register, which fills the column of rows nothing reads yet, and
    `burn_down_text` the burn-down register, which fills the last column."""
    tables = registry_tables(plan_text)
    unread = unread_rows(wiring_text)
    burn_down = burn_down_rows(burn_down_text)
    bad = problems(tables)
    if bad:
        raise SystemExit("the registry is not complete:\n  " + "\n  ".join(bad))
    total, tier_count, tweak_count, cross = _counts(tables)
    out = [
        "# Every stand-in value, with its source tier and its tweak level (GEN.725)",
        "",
        "Status: generated from the stand-in registry of [gap-filling-plan.md](gap-filling-plan.md) by `scripts/build_stand_in_values.py`;"
        " do not edit it by hand (`python scripts/build_stand_in_values.py --write` rewrites it, and a test fails when it differs"
        " from the registry). The tier and the tweak level of each row are the registry's; a row's owner confirms them.",
        "",
        "A stand-in is any value that is not a verified measurement of the thing it stands for. Every one carries two labels"
        " (Boss, 3:40 PM EDT on 2026-10-24; the tweak scale decided 4:01 PM EDT): **where its number came from** (the source tier)"
        " and **how far it was bent to give the outcome wanted** (the tweak level). Any change that adds a stand-in adds its row"
        " to the registry with both labels in the same pull request.",
        "",
        f"There are {total} rows below. A row may be a family of values (a recipe table, a set of bands); its tier cell then names the tier"
        " of each part, and its tweak is the highest level any part has.",
        "",
        "## The source tiers",
        "",
        "| Tier | Name | Meaning |",
        "|---|---|---|",
    ]
    out += [f"| {code} | {name} | {meaning} |" for code, name, meaning in TIERS]
    out += ["", "Unverified and guessed tiers may be tweaked to hit target probabilities (Boss, 3:40 PM EDT on 2026-10-24).",
            "", "## The tweak levels", "", "| Level | Name | Meaning |", "|---|---|---|"]
    out += [f"| {code} | {name} | {meaning} |" for code, name, meaning in TWEAKS]
    out += ["", "## How many rows are at each tier and tweak level", "",
            "| Tier | " + " | ".join(code for code, _, _ in TWEAKS) + " | All |",
            "|---|" + "---|" * (len(TWEAKS) + 1)]
    for code, _, _ in TIERS:
        out.append(f"| {code} | " + " | ".join(str(cross[(code, w)]) for w, _, _ in TWEAKS) + f" | {tier_count[code]} |")
    out.append("| All | " + " | ".join(str(tweak_count[w]) for w, _, _ in TWEAKS) + f" | {total} |")
    out += [
        "",
        "## How the labels were given",
        "",
        "- **Tier** follows the kind of the fill and the evidence label of its inputs: a derived (D) fill from a source we read is T1 or"
        " T2; a bracket (B) is T2 or T3 by how many sources agree and whether we opened them; a theory (T) fill is T5 to T7 by how"
        " much of the physics is derived and how much inferred; a placeholder (P) with no source is T4 or, when physics argues its"
        " size, T6 or T7. Rows whose parts differ list each part.",
        "- **Tweak** is L0 unless the value was moved after the fact. The L3 rows (hand-tuned against a sweep) are the mantle redox"
        " laws AT.a1 and a2, the atmosphere draw AT.b5, BC6 and AM19 (the exobase ratios); the L2 row is AM10 (the alkali"
        " activity, tuned to Kite et al. 2016); L1 rows are nudged inside the source's error bar. No row is L4: a row at L4 needs the"
        " item that justified it.",
        "- The recipe rows of the atmosphere draw (`atmosphere_recipes_data.py`) are not rows here one by one. Those that carry a"
        " tuning target in `tuning_targets.json` (a3, a4, e1, f5, f7, k1 to k5, k7, o2) are L2: their weights were tuned to a published"
        " statistic.",
        "",
    ]
    for number, title, columns in SECTIONS:
        header, rows = tables[number]
        out += [f"## {title}", "", "| " + " | ".join(columns + (UNREAD_COLUMN, BURN_DOWN_COLUMN)) + " |", "|" + "---|" * (len(columns) + 2)]
        for row in rows:
            cells = dict(zip(header, row))
            shown = [(_first_option(cells[name]) if name == "Fill (ranked)" else cells[name]) for name in columns]
            shown.append(unread.get(f"{number}:{row[0]}", "-"))
            disposition = burn_down.get(f"{number}:{row[0]}")
            shown.append("-" if disposition is None else disposition[0] + (f" ({disposition[1]})" if disposition[1] != "-" else ""))
            out.append("| " + " | ".join(shown) + " |")
        out.append("")
    return "\n".join(out).rstrip("\n") + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--check", action="store_true", help="fail if the committed list differs from the registry")
    args = parser.parse_args(argv)
    with open(PLAN, encoding="utf-8") as handle, open(WIRING, encoding="utf-8") as wiring, \
            open(BURN_DOWN, encoding="utf-8") as burn_down:
        text = render(handle.read(), wiring.read(), burn_down.read())
    if args.check:
        current = open(OUT, encoding="utf-8").read() if os.path.exists(OUT) else ""
        if current != text:
            print("docs/design/research/stand-in-values.md is out of date; run python scripts/build_stand_in_values.py", file=sys.stderr)
            return 1
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote {os.path.relpath(OUT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
