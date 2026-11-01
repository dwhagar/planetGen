# planetgen/scripts/materials_data.py

"""
The numbers typed in for the materials table (GEN.427), read by `build_materials_table.py`.

Everything here was read from a named source on 2026-10-22 (the Handbook of Mineralogy pages, PubChem annotations, Seager et
al. 2007 Tables 1 to 3) or is taken from a research part that read it (`docs/design/composition/body-composition-parts/`). The
raw pages and papers are not kept in the repository; the values and the place they were read are.
"""

HANDBOOK = "handbook-mineralogy"
WEBMINERAL = "webmineral"

# id: (value g/cm3, measured low, measured high or None, source key, how the value was taken)
# "D(calc.)" is the density computed from the cell and formula; where the Handbook gives none the value is the middle
# of the measured range or the single measured value, and the basis says so.
MINERAL_DENSITY = {
    "troilite": (4.85, 4.67, 4.79, HANDBOOK, "D(calc.)"),
    "schreibersite": (7.28, 7.0, 7.3, HANDBOOK, "middle of D(calc.) 7.12 to 7.44"),
    "cohenite": (7.68, 7.20, 7.65, HANDBOOK, "D(calc.)"),
    "magnetite": (5.20, 5.175, 5.175, HANDBOOK, "D(calc.); page is a scan, read as 5.20"),
    "hematite": (5.255, 5.26, 5.26, HANDBOOK, "D(calc.)"),
    "chromite": (5.12, 4.5, 4.8, HANDBOOK, "D(calc.) of the ideal formula (printed in brackets); page is a scan"),
    "sperrylite": (10.78, 10.58, 10.58, HANDBOOK, "D(calc.)"),
    "cooperite": (10.2, 9.5, 9.5, HANDBOOK, "D(calc.)"),
    "braggite": (9.383, None, None, HANDBOOK, "D(calc.); measured is printed as about 10"),
    "laurite": (6.39, 6.43, 6.43, HANDBOOK, "D(calc.)"),
    "vysotskite": (6.705, None, None, HANDBOOK, "D(calc.); no measurement"),
    "uraninite": (10.88, 10.63, 10.95, HANDBOOK, "D(calc.); falls to 6.5 with oxidation"),
    "thorianite": (9.991, 9.7, 9.8, HANDBOOK, "D(calc.)"),
    "carnotite": (4.91, 4.70, 4.70, HANDBOOK, "D(calc.)"),
    "autunite": (3.14, 3.05, 3.2, HANDBOOK, "D(calc.); varies with hydration"),
    "brannerite": (5.20, 4.2, 5.43, HANDBOOK, "D(calc.) printed in brackets; page is a scan"),
    "coffinite": (5.1, 5.1, 5.1, HANDBOOK, "the one measured value (no calculated density printed)"),
    "torbernite": (3.2, 3.2, 3.2, WEBMINERAL, "the single density listed"),
    "rutile": (4.25, 4.23, 4.23, HANDBOOK, "D(calc.); the measured value rises with Nb and Ta"),
    "ilmenite": (4.789, 4.72, 4.72, HANDBOOK, "D(calc.)"),
    "titanite": (3.53, 3.48, 3.60, HANDBOOK, "D(calc.) printed in brackets"),
    "perovskite": (4.02, 3.98, 4.26, HANDBOOK, "D(calc.)"),
    "brookite": (4.133, 4.08, 4.18, HANDBOOK, "D(calc.)"),
    "anatase": (3.89, 3.79, 3.97, HANDBOOK, "D(calc.) printed in brackets"),
    "quartz": (2.66, 2.65, 2.65, HANDBOOK, "D(calc.)"),
    "corundum": (3.997, 3.98, 4.10, HANDBOOK, "D(calc.)"),
    "beryl": (2.640, 2.63, 2.97, HANDBOOK, "D(calc.); rises with alkali content"),
    "topaz": (3.55, 3.49, 3.57, HANDBOOK, "D(calc.)"),
    "spinel": (3.64, 3.57, 3.72, WEBMINERAL, "the listed average"),
    "zircon": (4.714, 4.6, 4.7, HANDBOOK, "D(calc.)"),
    "schorl": (3.244, 3.18, 3.22, HANDBOOK, "D(calc.)"),
    "orthoclase": (2.563, 2.55, 2.63, HANDBOOK, "D(calc.)"),
    "calcite": (2.711, 2.7102, 2.7102, HANDBOOK, "D(calc.)"),
    "halite": (2.165, 2.168, 2.168, HANDBOOK, "D(calc.)"),
    "fluorite": (3.180, 3.175, 3.184, HANDBOOK, "D(calc.); to 3.56 if rich in rare earths"),
    "fluorapatite": (3.18, 3.1, 3.25, HANDBOOK, "D(calc.)"),
    "turquoise": (2.91, 2.86, 2.86, HANDBOOK, "D(calc.); less if massive"),
    "malachite": (3.983, 4.05, 4.05, HANDBOOK, "D(calc.)"),
    "almandine": (4.313, 4.318, 4.318, HANDBOOK, "D(calc.)"),
    "pyrope": (3.74, 3.65, 3.84, WEBMINERAL, "the listed average"),
    "spessartine": (4.179, 4.190, 4.190, HANDBOOK, "D(calc.)"),
    "grossular": (3.594, 3.594, 3.594, HANDBOOK, "D(calc.)"),
    "muscovite": (2.83, 2.77, 2.88, HANDBOOK, "D(calc.)"),
    "biotite": (3.25, 2.7, 3.3, HANDBOOK, "D(calc.)"),
    "chrysoberyl": (3.69, 3.75, 3.75, HANDBOOK, "D(calc.)"),
    "serpentine": (2.57, 2.55, 2.55, HANDBOOK, "D(calc.) of lizardite"),
    "forsterite": (3.271, 3.275, 3.275, HANDBOOK, "D(calc.)"),
    "fayalite": (4.40, 4.392, 4.392, HANDBOOK, "D(calc.) printed in brackets"),
    "enstatite": (3.189, 3.2, 3.9, HANDBOOK, "D(calc.)"),
    "ferrosilite": (3.87, 3.96, 3.96, HANDBOOK, "D(calc.) printed in brackets"),
    "diopside": (3.278, 3.22, 3.38, HANDBOOK, "D(calc.)"),
    "albite": (2.615, 2.60, 2.65, HANDBOOK, "middle of D(calc.) 2.609 to 2.621"),
    "anorthite": (2.760, 2.74, 2.76, HANDBOOK, "D(calc.)"),
    "gypsum": (2.31, 2.317, 2.317, HANDBOOK, "D(calc.)"),
    "kaolinite": (2.63, 2.61, 2.68, HANDBOOK, "D(calc.)"),
    "goethite": (4.18, 4.28, 4.28, HANDBOOK, "D(calc.)"),
    "kamacite": (7.90, None, None, HANDBOOK, "D(calc.) printed in brackets; no measurement"),
    "taenite": (8.29, 7.8, 8.22, HANDBOOK, "D(calc.) of the ordered phase"),
    "osmiridium": (22.15, None, None, HANDBOOK, "D(calc.) of iridosmine (Os,Ir); no measurement"),
    "moissanite": (3.22, 3.218, 3.22, WEBMINERAL, "the listed range"),
    "periclase": (3.58, 3.56, 3.68, HANDBOOK, "D(calc.)"),
    "wustite": (5.97, 5.88, 5.88, HANDBOOK, "D(calc.)"),
    "bastnasite": (4.97, 4.95, 5.0, WEBMINERAL, "the listed average (bastnasite-(Ce))"),
    "monazite": (5.15, 4.8, 5.5, WEBMINERAL, "the listed average (monazite-(Ce))"),
    "xenotime": (4.75, 4.4, 5.1, WEBMINERAL, "the listed average (xenotime-(Y))"),
    "cerite": (4.86, 4.86, 4.86, WEBMINERAL, "the single density listed (cerite-(Ce))"),
    "gadolinite": (4.2, 4.2, 4.2, WEBMINERAL, "the single density listed (gadolinite-(Ce))"),
    "samarskite": (5.69, 5.6, 5.8, WEBMINERAL, "the listed average (samarskite-(Y))"),
    "fergusonite": (5.05, 4.3, 5.8, WEBMINERAL, "the listed average (fergusonite-(Y))"),
    "euxenite": (4.84, 4.7, 5.0, WEBMINERAL, "the listed average (euxenite-(Y))"),
}

# The pure elements: density in g/cm3 at room temperature of the common form, from the CRC-based "use" row of the Wikipedia
# data page "Densities of the elements (data page)" (partly verified: the compilation was read, not the handbook).
ELEMENT_DENSITY = {
    "carbon_graphite": 2.267, "carbon_diamond": 3.513, "magnesium": 1.738, "aluminum": 2.70, "silicon": 2.33,
    "phosphorus": 2.34, "sulfur": 2.08, "potassium": 0.89, "calcium": 1.55, "titanium": 4.506, "nickel": 8.908,
    "ruthenium": 12.45, "rhodium": 12.41, "palladium": 12.023, "tungsten": 19.25, "osmium": 22.59, "iridium": 22.56,
    "platinum": 21.45, "gold": 19.3,
}

# PubChem "Density" annotations (experimental, from the sources named there). id: (value g/cm3, temperature text)
PUBCHEM_DENSITY = {
    "naphthalene": (1.162, "20 C, solid"),
    "anthracene": (1.25, "27 C, solid"),
    "phenanthrene": (1.179, "25 C, solid"),
    "pyrene": (1.271, "23 C, solid"),
    "fluoranthene": (1.252, "0 C, solid"),
    "co_ice": (0.929, "65 K, hexagonal solid"),
    "n2_ice": (1.0265, "solid, temperature not given"),
    "co2_ice": (1.56, "about 194 K, solid"),
}

# Solid densities computed by part D (section 2.1) from the COD crystal cell of the same formula (label COMPUTED): id:
# (value g/cm3, basis)
COD_DENSITY = {
    "ch3oh_ice": (1.028, "122 K, COD crystal cell (0.991 at 160 K)"),
    "h2s_ice": (1.2, "COD crystal cell, one entry, about"),
}

# [R3992: Table 1]: K0 and K0' of the Vinet (V) or Birch-Murnaghan (BME) fits at
# zero pressure and 300 K, with the rho0 of the table. id: (model, order, rho0 kg/m3, K0 GPa, +-, K0', +-, note)
SEAGER_TABLE_1 = {
    "carbon_graphite": ("bme", 3, 2250.0, 33.8, 3.0, 8.9, 1.0, "Table 1, references 1 and 2"),
    "iron": ("bme", 3, 7860.0, 162.5, 5.0, 5.5, 0.8, "alpha iron, Table 1, references 1 and 3"),
    "iron_epsilon": ("vinet", None, 8300.0, 156.2, 1.8, 6.08, 0.12, "Table 1, reference 4"),
    "troilite": ("bme", 3, 4770.0, 35.0, 4.0, 5.0, 2.0, "FeS, Table 1, references 1 and 5; the low-pressure phase"),
    "periclase": ("bme", 3, 3560.0, 177.0, 4.0, 4.0, 0.1, "MgO, Table 1, reference 8; SLB's 161.14 is the other choice"),
    "enstatite": ("bme", 3, 3220.0, 125.0, None, 5.0, None, "MgSiO3 (en), Table 1; K0' is assumed there"),
    "mg_bridgmanite": ("bme", 4, 4100.0, 247.0, 4.0, 3.97, None, "MgSiO3 (pv), Table 1; fourth-order fit with K0'' -0.016 /GPa"),
    "mgfe_bridgmanite": ("bme", 3, 4260.0, 266.0, 6.0, 3.9, 0.04, "(Mg0.88,Fe0.12)SiO3 (pv), Table 1"),
    "moissanite": ("bme", 3, 3220.0, 227.0, 3.0, 4.1, 0.1, "SiC, Table 1, references 1 and 12"),
    "ice_vii": ("bme", 3, 1460.0, 23.7, 0.9, 4.15, 0.07, "H2O (ice VII), Table 1, reference 6"),
}
SEAGER_K0PP_PEROVSKITE = -0.016
"""float: Seager Table 1 note c: the fourth-order fit of MgSiO3 (pv) used K0'' = -0.016 per GPa."""
SEAGER_WATER_LIQUID_K0_GPA = 2.28
"""float: Table 1: liquid water, K0 2.28 GPa (seawater at 12 K as printed, reference 7)."""
SEAGER_ICE_X_VINET = (1460.0, 14.3771, 6.57972)
"""tuple: The Vinet fit to the DFT ices VIII and X of Seager et al. (rho0 kg/m3, K0 GPa, K0'), valid to 7.7 TPa (section 3.2)."""

# Table 3: the fits rho(P) = rho0 + c P^n merged with the Thomas-Fermi-Dirac EOS, P in Pa, valid below 1e16 Pa.
SEAGER_TABLE_3 = {
    "iron_epsilon": (8300.0, 0.00349, 0.528, "printed as Fe(alpha) in the table; the merge uses the epsilon Vinet fit (K0 156.2)"),
    "mg_bridgmanite": (4100.0, 0.00161, 0.541, "MgSiO3 (perovskite)"),
    "mgfe_bridgmanite": (4260.0, 0.00127, 0.549, "(Mg,Fe)SiO3"),
    "ice_x": (1460.0, 0.00311, 0.513, "H2O"),
    "carbon_graphite": (2250.0, 0.00350, 0.514, "C (graphite)"),
    "moissanite": (3220.0, 0.00172, 0.537, "SiC"),
}

# Table 2: the DFT equation of state of water ice VIII and X: pressure GPa, density kg/m3.
SEAGER_TABLE_2 = (
    (2.320, 1636.617), (4.155, 1725.860), (6.664, 1821.712), (9.823, 1924.800), (18.791, 2155.574), (25.361, 2284.820),
    (33.744, 2424.761), (44.314, 2576.285), (56.970, 2740.606), (74.188, 2919.275), (94.406, 3113.919), (126.815, 3326.163),
    (155.924, 3485.175), (240.696, 3867.031), (351.114, 4290.649), (498.660, 4760.933), (691.938, 5282.621),
    (937.585, 5861.450), (1260.182, 6503.954), (1673.049, 7216.474), (2188.301, 8006.625), (2853.712, 8884.186),
    (3691.387, 9858.497), (4737.211, 10938.603), (6040.611, 12134.882), (7686.171, 13467.635),
)
"""tuple: [R3992: Table 2]. The table prints the density in kg per cubic metre but its numbers are in Mg/m3; they are
multiplied by 1000 here (the table at 2.32 GPa reads 1.6366)."""

# Part B table 2.1b, the other iron and light-element EOS rows. id: (model, rho0 kg/m3, K0 GPa, K0', source key, label, note)
PART_B_IRON = (
    ("iron_epsilon", "vinet", 8264.2, 157.07, 5.6167, "part-b-fit",
     "the single Vinet fit the solver uses, 0.3 percent over 1 to 1000 GPa and 1.9 percent to 20 TPa (COMPUTED in part B)"),
    ("iron_epsilon", "vinet", 8270.0, 163.4, 5.38, "hakim2018", "room-temperature Vinet below 234.4 GPa (Dewaele et al., as in [R1752])"),
    ("iron_epsilon", "vinet", 8300.0, 160.2, 5.82, "hakim2018", "the Valencia parameters, [R1752] Table F.7"),
    ("iron_epsilon", "bme", 8341.0, 173.98, 5.297, "hakim2018", "the Dorn parameters, [R1752] Table F.7"),
    ("iron_epsilon", "bme", 8269.0, 172.7, 4.79, "hakim2018", "the Fei parameters, [R1752] Table F.7"),
    ("troilite", "holzapfel", 6118.0, 144.3, 4.79, "hakim2018", "[R3923] refit, [R1752] Table D.5; the high-pressure phase"),
    ("fesi", "holzapfel", 6522.0, 219.5, 4.29, "hakim2018", "[R3923] refit, [R1752] Table D.5"),
    ("cohenite", "holzapfel", 7980.0, 289.4, 3.81, "hakim2018", "Fe3C, [R3923] refit, [R1752] Table D.5"),
    ("wustite", "holzapfel", 5764.0, 152.7, 4.13, "hakim2018", "Fe0.95O, [R3923] refit, [R1752] Table D.5"),
)

# The ice phases of part C table 2.3 ([R2204] Gibbs splines, evaluated). Each row of the part's table is
# (phase id, P MPa, T K, rho kg/m3, K_T GPa, dK/dP, alpha 1/K) in the order the part prints them.
ICE_POINTS = (
    ("water_ice", 0.1, 273.15, 916.72, 8.49, 4.89, 1.60e-4),
    ("water_ice", 0.1, 200.0, 926.13, 9.46, 6.38, 1.16e-4),
    ("water_ice", 0.1, 100.0, 933.05, 10.45, 8.07, 3.0e-5),
    ("water_ice", 200.0, 251.15, 939.82, 9.86, 5.39, 1.13e-4),
    ("ice_iii", 209.5, 251.15, 1154.5, 9.19, 6.25, 2.2e-4),
    ("ice_iii", 355.0, 256.43, 1170.8, 10.01, 6.00, 2.05e-4),
    ("ice_ii", 300.0, 200.0, 1196.0, 14.89, 5.95, 1.14e-4),
    ("ice_v", 355.0, 256.43, 1244.9, 12.92, 6.15, 1.8e-4),
    ("ice_v", 618.4, 272.73, 1265.7, 14.24, 5.89, 1.68e-4),
    ("ice_vi", 618.4, 272.73, 1327.2, 14.87, 6.73, 2.44e-4),
    ("ice_vi", 1000.0, 300.0, 1351.1, 16.64, 6.34, 2.26e-4),
    ("ice_vi", 2216.0, 355.0, 1422.5, 22.44, 5.53, 1.78e-4),
    ("ice_vii", 2216.0, 355.0, 1476.4, 26.57, 2.11, 2.6e-4),
    ("ice_vii", 3000.0, 300.0, 1537.6, 29.29, 2.37, 2.0e-4),
    ("ice_vii", 10000.0, 300.0, 1845.3, 50.85, 3.53, 1.0e-4),
    ("ice_x", 50000.0, 300.0, 2674.6, 204.8, 5.11, 2.0e-5),
    ("ice_x", 100000.0, 300.0, 3107.1, 503.8, 4.62, 7.4e-6),
    ("ice_x", 500000.0, 300.0, 4669.7, 1708.0, 2.79, 1.6e-6),
    ("ice_x", 1000000.0, 300.0, 5814.4, 2880.0, 1.91, 7.8e-7),
)
ICE_VALID = {
    "water_ice": {"p_mpa": [0.0, 400.0], "t_k": [1.0, 301.0]},
    "ice_ii": {"p_mpa": [0.0, 900.0]},
    "ice_iii": {"p_mpa": [0.0, 500.0]},
    "ice_v": {"p_mpa": [0.0, 1000.0]},
    "ice_vi": {"p_mpa": [0.0, 3000.0], "t_k": [0.0, 400.0]},
    "ice_vii": {"p_mpa": [1700.0, 1000000.0], "t_k": [20.0, 1800.0]},
    "ice_x": {"p_mpa": [1700.0, 1000000.0], "t_k": [20.0, 1800.0]},
}
LIQUID_WATER_POINT = (0.1, 300.0, 996.56)
"""tuple: Liquid water at 1 bar and 300 K, kg/m3 (Bollengier spline, part C): the check value of the table."""

# Class composition texts (tuning.PLANET_CLASSES) and what each names, as ids or legacy names.
CLASS_COMPOSITIONS = {
    "igneous silica and basalt": ["igneous silica", "basalt"],
    "iron, potassium, silicon, and sulfur": ["iron", "potassium", "silicon", "sulfur"],
    "anthracite, basalt, and hydrocarbons": ["anthracite", "basalt", "hydrocarbons"],
    "frozen hydrocarbons and ice": ["frozen hydrocarbons", "ice"],
    "silicon, iron, and magnesium": ["silicon", "iron", "magnesium"],
    "rock, ice, methane, and ammonia": ["rock", "ice", "methane", "ammonia"],
    "hydrogen and helium": ["hydrogen", "helium"],
    "silicon, iron, magnesium, and ice": ["silicon", "iron", "magnesium", "ice"],
    "silicate rock over a large iron core": ["silicate rock", "iron"],
    "hydrogen, helium, and hydrocarbons": ["hydrogen", "helium", "hydrocarbons"],
    "iron, iridium, tungsten, and nickel": ["iron", "iridium", "tungsten", "nickel"],
}
