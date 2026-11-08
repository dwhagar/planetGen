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
    # GEN.1356: the Handbook of Mineralogy pages of the sulfides, nitrides and carbides, read 2026-11-07.
    "pyrite": (5.013, 5.018, 5.018, HANDBOOK, "D(calc.)"),
    "pyrrhotite": (4.69, 4.58, 4.65, HANDBOOK, "D(calc.) of the Fe7S8 cell (a 6.865, b 11.9, c 22.72 A, 8 formula units)"),
    "pentlandite": (4.956, 4.6, 5.0, HANDBOOK, "D(calc.) as printed; the cell (a 9.928 A, Z 4) and a 1:1 iron-nickel formula give about 5.2"),
    "chalcopyrite": (4.283, 4.1, 4.3, HANDBOOK, "D(calc.)"),
    "oldhamite": (2.589, 2.58, 2.58, HANDBOOK, "D(calc.)"),
    "niningerite": (3.40, 3.21, 3.59, HANDBOOK, "middle of D(calc.) 3.21 to 3.59 (no measurement; the range is the composition range)"),
    "alabandite": (4.053, 3.95, 4.04, HANDBOOK, "D(calc.)"),
    "daubreelite": (3.842, 3.81, 3.81, HANDBOOK, "D(calc.)"),
    "osbornite": (5.38, None, None, HANDBOOK, "D(calc.) printed in brackets (synthetic TiN, a 4.24173 A, Z 4); no measurement"),
    "khamrabaevite": (4.97, None, None, HANDBOOK,
                      "COMPUTED from the Handbook's cell (rock salt, a 4.319 A, Z 4) and analysis (metal mass 48.23 g/mol, one carbon per metal); "
                      "the page prints D(calc.) 10.01, the same cell counted with Z 8"),
    "sinoite": (2.84, 2.80, 2.85, HANDBOOK, "D(calc.)"),
}

# GEN.1356 and the mineral rows after it: id -> (crystal system, Mohs middle, Mohs low, Mohs high, color). The Mohs numbers are the
# Handbook's printed hardness range (its middle is the value); None where the page prints none (n.d.). Colors in the page's words.
MINERAL_FIELDS = {
    "pyrite": ("cubic", 6.25, 6.0, 6.5, "pale brass-yellow, tarnishes darker and iridescent"),
    "pyrrhotite": ("monoclinic", 4.0, 3.5, 4.5, "bronze-yellow to pinchbeck-brown"),
    "pentlandite": ("cubic", 3.75, 3.5, 4.0, "pale bronze-yellow"),
    "chalcopyrite": ("tetragonal", 3.75, 3.5, 4.0, "brass-yellow"),
    "oldhamite": ("cubic", 4.0, 4.0, 4.0, "pale chestnut-brown"),
    "niningerite": ("cubic", None, None, None, "gray in reflected light"),
    "alabandite": ("cubic", 3.75, 3.5, 4.0, "iron-black, tarnishing brown"),
    "daubreelite": ("cubic", 5.0, 5.0, 5.0, "black"),
    "osbornite": ("cubic", None, None, None, "golden yellow"),
    "khamrabaevite": ("cubic", 9.5, 9.0, 10.0, "dark gray"),
    "sinoite": ("orthorhombic", None, None, None, "colorless in thin section, light gray in reflected light"),
    "troilite": ("hexagonal", 4.0, 3.5, 4.5, "pale grayish brown, tarnishes rapidly"),
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
    "rock, ice, and frozen hydrocarbons": ["rock", "ice", "frozen hydrocarbons"],
    "silicon, iron, and magnesium": ["silicon", "iron", "magnesium"],
    "rock, ice, methane, and ammonia": ["rock", "ice", "methane", "ammonia"],
    "hydrogen and helium": ["hydrogen", "helium"],
    "silicon, iron, magnesium, and ice": ["silicon", "iron", "magnesium", "ice"],
    "silicon, iron, magnesium, ice, and hydrocarbons": ["silicon", "iron", "magnesium", "ice", "hydrocarbons"],
    "silicon, iron, magnesium, water ice, and ammonia": ["silicon", "iron", "magnesium", "water ice", "ammonia"],
    "silicon, iron, magnesium, and sulfur": ["silicon", "iron", "magnesium", "sulfur"],
    "silicate rock over a large iron core": ["silicate rock", "iron"],
    "hydrogen, helium, and hydrocarbons": ["hydrogen", "helium", "hydrocarbons"],
    "iron, iridium, tungsten, and nickel": ["iron", "iridium", "tungsten", "nickel"],
}

# GEN.1355: the silicates and oxides the ground recipes and the norm name, read from the Webmineral data page of each mineral on 2026-11-07 (the research files
# real-material-data/minerals_missing_from_table.csv and part 3 section 1 carry the same values; the pages were opened again for this table).
# id: (formula the table reads or None, family, aliases, the formula as the page prints it, crystal system, Mohs (middle, low, high) or None, density g/cm3,
#      how the density was taken, measured range or None, color as listed, the page read)
# A formula of None is a solid solution, a clay of variable water and charge, or a mixture the table cannot give one molar mass: the page's formula is the composition.
# Where the table reads a formula it is the end member the page names. Wollastonite and hornblende were recalled in the research file: the page of
# Wollastonite-1A lists an average of 2.84 (not 2.9) and webmineral redirects hornblende to magnesiohornblende, whose page is read for the series.
PART1_MINERALS = {
    'akaganeite': (None, 'oxidized', [], 'Fe+++(O,OH,Cl)', 'monoclinic', None, 3.0, 'meas', None, 'Brown, Rusty brown', 'Akaganeite'),
    'cronstedtite': (None, 'clay', [], 'Fe++2Fe+++(SiFe+++)O5(OH)4', 'trigonal', (3.5, 3.5, 3.5), 3.34, 'meas avg', (3.34, 3.35), 'Brownish black, Greenish black, Dark brown, Black', 'Cronstedtite'),
    'ferrihydrite': (None, 'oxidized', [], 'Fe+++2O3.0.5(H2O)', 'trigonal', None, 3.8, 'meas', None, 'Yellow brown, Dark brown', 'Ferrihydrite'),
    'lepidocrocite': ('FeO(OH)', 'oxidized', [], 'FeO(OH)', 'orthorhombic', (5.0, 5.0, 5.0), 4.0, 'listed value (read; the research file had it recalled)', None, 'Red, Yellowish brown, Blackish brown', 'Lepidocrocite'),
    'maghemite': ('Fe2O3', 'oxidized', [], 'gamma-Fe+++2O3', 'cubic', (6.0, 6.0, 6.0), 4.9, 'meas', None, 'Dark brown', 'Maghemite'),
    'magnesiochromite': ('MgCr2O4', 'oxide', [], 'MgCr2O4', 'cubic', (5.5, 5.5, 5.5), 4.2, 'meas', None, 'Black', 'Magnesiochromite'),
    'schwertmannite': ('Fe16O16(OH)12(SO4)2', 'oxidized', [], 'Fe+++16O16(OH)12(SO4)2', 'tetragonal', (3.0, 2.5, 3.5), 3.88, 'meas avg', (3.77, 3.99), 'Brownish yellow', 'Schwertmannite'),
    'tochilinite': (None, 'sulfide', [], '6Fe0.9S.5(Mg,Fe++)(OH)2', 'monoclinic', (1.25, 1.0, 1.5), 2.97, 'meas avg', (2.96, 2.99), 'Bronze black', 'Tochilinite'),
    'ulvospinel': ('Fe2TiO4', 'oxide', [], 'TiFe++2O4', 'cubic', (5.75, 5.5, 6.0), 4.91, 'calc', None, 'Black', 'Ulvospinel'),
    'antigorite': ('Mg3Si2O5(OH)4', 'clay', [], '(Mg,Fe++)3Si2O5(OH)4', 'monoclinic', (3.75, 3.5, 4.0), 2.54, 'meas avg', (2.5, 2.6), 'Green, Gray, Bluish gray, Brown, Black', 'Antigorite'),
    'beidellite': (None, 'clay', [], 'Na0.5Al2(Si3.5Al0.5)O10(OH)2.n(H2O)', 'monoclinic', (1.5, 1.0, 2.0), 2.15, 'meas avg', (2.0, 2.3), 'White, Gray, Reddish white, Brownish white, Pale yellow', 'Beidellite'),
    'brucite': ('Mg(OH)2', 'oxide', [], 'Mg(OH)2', 'trigonal', (2.75, 2.5, 3.0), 2.39, 'meas avg', (2.39, 2.4), 'Blue, Gray, Gray blue, Yellow, White', 'Brucite'),
    'chamosite': (None, 'clay', [], '(Fe++,Mg,Fe+++)5Al(Si3Al)O10(OH,O)8', 'monoclinic', (3.0, 3.0, 3.0), 3.2, 'meas avg', (3.0, 3.4), 'Gray, Greenish gray, Brown, Blackish green', 'Chamosite'),
    'chrysotile': ('Mg3Si2O5(OH)4', 'clay', [], 'Mg3Si2O5(OH)4', 'monoclinic', (2.5, 2.5, 2.5), 2.53, 'meas', None, 'Green', 'Chrysotile'),
    'clinochlore': ('Mg5Al2Si3O10(OH)8', 'clay', [], '(Mg,Fe++)5Al(Si3Al)O10(OH)8', 'monoclinic', (2.25, 2.0, 2.5), 2.65, 'meas avg', (2.55, 2.75), 'Blackish green, Bluish green, White, Yellowish green, Olive ', 'Clinochlore'),
    'illite': (None, 'clay', [], '(K,H3O)(Al,Mg,Fe)2(Si,Al)4O10[(OH)2,(H2O)]', 'monoclinic', (1.5, 1.0, 2.0), 2.75, 'meas avg', (2.6, 2.9), 'White', 'Illite'),
    'lizardite': ('Mg3Si2O5(OH)4', 'clay', [], 'Mg3Si2O5(OH)4', 'triclinic', (2.5, 2.5, 2.5), 2.57, 'meas avg', (2.55, 2.6), 'Green, Green blue, Yellow, White', 'Lizardite'),
    'montmorillonite': (None, 'clay', [], '(Na,Ca)0,3(Al,Mg)2Si4O10(OH)2.n(H2O)', 'monoclinic', (1.75, 1.5, 2.0), 2.35, 'meas avg', (2.0, 2.7), 'White, Gray white, Yellow, Brownish yellow, Greenish yellow', 'Montmorillonite'),
    'nontronite': (None, 'clay', [], 'Na0.3Fe+++2(Si,Al)4O10(OH)2.n(H2O)', 'monoclinic', (1.75, 1.5, 2.0), 2.3, 'meas', None, 'Greenish yellow, Brownish green, Olive green', 'Nontronite'),
    'palygorskite': (None, 'clay', [], '(Mg,Al)2Si4O10(OH).4(H2O)', 'monoclinic', (2.25, 2.0, 2.5), 2.15, 'meas avg', (2.1, 2.2), 'White, Gray, Brownish white', 'Palygorskite'),
    'saponite': (None, 'clay', [], '(Ca/2,Na)0,3(Mg,Fe++)3(Si,Al)4O10(OH)2.4(H2O)', 'monoclinic', (1.75, 1.5, 2.0), 2.3, 'meas', None, 'White, Yellowish white, Greenish white, Reddish white, Bluis', 'Saponite'),
    'sepiolite': ('Mg4Si6O15(OH)2.6H2O', 'clay', [], 'Mg4Si6O15(OH)2.6(H2O)', 'orthorhombic', (2.0, 2.0, 2.0), 2.0, 'meas', None, 'Grayish white, Yellowish white, Bluish green white, Reddish ', 'Sepiolite'),
    'talc': ('Mg3Si4O10(OH)2', 'clay', [], 'Mg3Si4O10(OH)2', 'monoclinic', (1.0, 1.0, 1.0), 2.75, 'meas avg', (2.7, 2.8), 'Pale green, White, Gray white, Yellowish white, Brownish whi', 'Talc'),
    'cristobalite': ('SiO2', 'oxide', [], 'SiO2', 'tetragonal', (6.5, 6.5, 6.5), 2.27, 'meas', None, 'Blue gray, Brown, Gray, Yellow, White', 'Cristobalite'),
    'lime': ('CaO', 'oxide', ['calcium oxide'], 'CaO', 'cubic', (3.5, 3.5, 3.5), 3.35, 'meas', None, 'White', 'Lime'),
    'lonsdaleite': ('C', 'carbon', [], 'C', 'hexagonal', (7.5, 7.0, 8.0), 3.41, 'meas avg', (3.3, 3.52), 'Brownish black, Light brownish yellow', 'Lonsdaleite'),
    'opal': (None, 'oxide', [], 'SiO2.n(H2O)', 'amorphous', (5.75, 5.5, 6.0), 2.09, 'meas avg', (1.9, 2.3), 'White, Yellow, Red, Brown, Blue', 'Opal'),
    'portlandite': ('Ca(OH)2', 'oxide', [], 'Ca(OH)2', 'trigonal', (2.75, 2.5, 3.0), 2.23, 'meas', None, 'Colorless', 'Portlandite'),
    'tridymite': ('SiO2', 'oxide', [], 'SiO2', 'triclinic', (6.75, 6.5, 7.0), 2.3, 'meas avg', (2.28, 2.33), 'Colorless, White, Yellowish white, Gray', 'Tridymite'),
    'akermanite': ('Ca2MgSi2O7', 'silicate', [], 'Ca2MgSi2O7', 'tetragonal', (5.5, 5.0, 6.0), 2.94, 'meas', None, 'Colorless, Grayish, Green, Brown', 'Akermanite'),
    'gehlenite': ('Ca2Al2SiO7', 'silicate', [], 'Ca2Al(AlSi)O7', 'tetragonal', (5.5, 5.0, 6.0), 2.98, 'meas avg', (2.9, 3.07), 'Brown, Colorless, Greenish, Gray, Yellow', 'Gehlenite'),
    'grossmanite': ('CaTiAlSiO6', 'silicate', [], 'CaTi+++AlSiO6', 'monoclinic', None, 3.41, 'calc', None, '', 'Grossmanite'),
    'hibonite': (None, 'oxide', [], '(Ca,Ce)(Al,Ti,Mg)12O19', 'hexagonal', (7.75, 7.5, 8.0), 3.84, 'meas', None, 'Black, Black, Reddish brown', 'Hibonite'),
    'melilite': (None, 'silicate', [], '(Ca,Na)2(Al,Mg,Fe++)(Si,Al)2O7', 'tetragonal', (5.25, 5.0, 5.5), 2.95, 'meas avg', (2.9, 3.0), 'White, Yellowish white, Brown, Gray, Greenish gray', 'Melilite'),
    'analcime': ('NaAlSi2O6.H2O', 'silicate', [], 'NaAlSi2O6.(H2O)', 'triclinic', (5.0, 5.0, 5.0), 2.3, 'meas', None, 'White, Grayish white, Greenish white, Yellowish white, Reddi', 'Analcime'),
    'augite': (None, 'silicate', [], '(Ca,Na)(Mg,Fe,Al,Ti)(Si,Al)2O6', 'monoclinic', (5.75, 5.0, 6.5), 3.4, 'meas avg', (3.2, 3.6), 'Brown green, Green, Light brown, Dark brown, Black', 'Augite'),
    'epidote': (None, 'silicate', [], 'Ca2(Fe+++,Al)3(SiO4)3(OH) = Ca2(Fe,Al)Al2(SiO4)(Si2O7)O(OH)', 'monoclinic', (7.0, 7.0, 7.0), 3.45, 'meas avg', (3.3, 3.6), 'Yellowish green, Brownish green, Black, Yellow, Gray', 'Epidote'),
    'hornblende': (None, 'silicate', [], 'Ca2(Mg,Fe,Al)5(Si,Al)8O22(OH)2', 'monoclinic', (5.5, 5.0, 6.0), 3.23, 'listed average of 3 to 3.47 (the page of magnesiohornblende, where webmineral redirects hornblende)', (3.0, 3.47), 'Brown, Green, Greenish brown, Greenish black, Dark green', 'Magnesiohornblende'),
    'hypersthene': (None, 'silicate', [], '(Mg,Fe++)2Si2O6', 'orthorhombic', (5.75, 5.5, 6.0), 3.55, 'meas avg', (3.2, 3.9), 'Grayish white, Greenish white, Yellowish white, Bronzy brown', 'Hypersthene'),
    'jadeite': ('NaAlSi2O6', 'silicate', [], 'Na(Al,Fe+++)Si2O6', 'monoclinic', (6.5, 6.5, 6.5), 3.3, 'meas avg', (3.25, 3.35), 'Green, White, Pale bluish gray, Grayish green, Pale purple', 'Jadeite'),
    'leucite': ('KAlSi2O6', 'silicate', [], 'KAlSi2O6', 'tetragonal', (6.0, 6.0, 6.0), 2.47, 'meas', None, 'Colorless, Gray, Yellow gray, White', 'Leucite'),
    'nepheline': ('NaAlSiO4', 'silicate', [], '(Na,K)AlSiO4', 'hexagonal', (6.0, 6.0, 6.0), 2.59, 'meas avg', (2.55, 2.65), 'White, Gray, Brown, Brownish gray, Reddish white', 'Nepheline'),
    'pigeonite': (None, 'silicate', [], '(Mg,Fe++,Ca)(Mg,Fe++)Si2O6', 'monoclinic', (6.0, 6.0, 6.0), 3.38, 'meas avg', (3.3, 3.46), 'Brown, Greenish brown, Light purple brown, Black', 'Pigeonite'),
    'tremolite': ('Ca2Mg5Si8O22(OH)2', 'silicate', [], '[ ]Ca2Mg5Si8O22(OH)2', 'monoclinic', (5.5, 5.0, 6.0), 3.05, 'meas avg', (2.9, 3.2), 'Brown, Colorless, Gray, White, Light green', 'Tremolite'),
    'wollastonite': ('CaSiO3', 'silicate', [], 'CaSiO3', 'triclinic', (5.0, 5.0, 5.0), 2.84, 'listed average of 2.8 to 2.9 (Wollastonite-1A page; the research file had a recalled 2.9)', (2.8, 2.9), 'White, Yellow, Gray, Red, Brown', 'Wollastonite-1A'),
}

# Groups of the rows above (and the amphibole, chlorite and smectite names the recipes use): group id -> (name, member ids). A group carries no number.
PART1_GROUPS = {
    "smectite": ("smectite", ("montmorillonite", "saponite", "nontronite", "beidellite")),
    "chlorite": ("chlorite", ("clinochlore", "chamosite")),
    "amphibole": ("amphibole", ("tremolite", "hornblende")),
}

# GEN.1357: the carbonates, salts and phosphates the evaporite and ore recipes will draw, read from the Webmineral data page of each mineral on 2026-11-07 (the two
# perchlorates, which webmineral does not list, from PubChem). Same fields as PART1_MINERALS; a page starting "PubChem" is read from the PubChem record.
# Values the research file had recalled and the pages changed: rhodochrosite Mohs 3 (not 4); chlorapatite 3.15 (the average of 3.1 to 3.2, not 3.17) and
# hydroxylapatite 3.08 (not 3.16), on the pages of Apatite-(CaCl) and Apatite-(CaOH) where webmineral redirects the older names; magnesium perchlorate 2.21 (PubChem,
# USCG), not 2.6. Sylvite 1.99, langbeinite 2.83, rhodochrosite 3.69 and calcium perchlorate 2.651 were opened and stand.
PART3_MINERALS = {
    'ankerite': (None, 'carbonate', [], 'Ca(Fe++,Mg,Mn)(CO3)2', 'trigonal', (3.75, 3.5, 4.0), 3.05, 'the listed average of 3 to 3.1', (3.0, 3.1), 'Brown, Brownish yellow, Gray, Gray brown, White', 'Ankerite'),
    'aragonite': ('CaCO3', 'carbonate', [], '', 'orthorhombic', (3.75, 3.5, 4.0), 2.93, 'the single value listed', None, 'Colorless, White, Gray, Yellowish white, Reddish white', 'Aragonite'),
    'dolomite': ('CaMg(CO3)2', 'carbonate', [], '', 'trigonal', (3.75, 3.5, 4.0), 2.84, 'the listed average of 2.8 to 2.9', (2.8, 2.9), 'White, Gray, Reddish white, Brownish white, Gray', 'Dolomite'),
    'hydromagnesite': ('Mg5(CO3)4(OH)2.4H2O', 'carbonate', [], 'Mg5(CO3)4(OH)2.4(H2O)', 'monoclinic', (3.5, 3.5, 3.5), 2.18, 'the listed average of 2.16 to 2.2', (2.16, 2.2), 'Colorless, White', 'Hydromagnesite'),
    'magnesite': ('MgCO3', 'carbonate', [], '', 'trigonal', (4.0, 4.0, 4.0), 3.0, 'the single value listed', None, 'Colorless, White, Grayish white, Yellowish white, Brownish white', 'Magnesite'),
    'nesquehonite': ('Mg(HCO3)(OH).2H2O', 'carbonate', [], 'Mg(HCO3)(OH).2(H2O)', 'monoclinic', (2.5, 2.5, 2.5), 1.84, 'the listed average of 1.83 to 1.85', (1.83, 1.85), 'Colorless, White', 'Nesquehonite'),
    'rhodochrosite': ('MnCO3', 'carbonate', [], '', 'trigonal', (3.0, 3.0, 3.0), 3.69, 'the single value listed', None, 'Pinkish red, Red, Rose red, Yellowish gray, Brown', 'Rhodochrosite'),
    'siderite': ('FeCO3', 'carbonate', [], 'Fe++CO3', 'trigonal', (3.5, 3.5, 3.5), 3.96, 'the single value listed', None, 'Yellowish brown, Brown, Gray, Yellowish gray, Greenish gray', 'Siderite'),
    'alunite': ('KAl3(SO4)2(OH)6', 'evaporite', [], '', 'trigonal', (3.75, 3.5, 4.0), 2.74, 'the listed average of 2.59 to 2.9', (2.59, 2.9), 'White, Gray, Yellowish gray, Reddish gray, Yellowish white', 'Alunite'),
    'anhydrite': ('CaSO4', 'evaporite', [], '', 'orthorhombic', (3.5, 3.5, 3.5), 2.97, 'the listed average of 2.96 to 2.98', (2.96, 2.98), 'Colorless, White, Bluish white, Violet white, Dark gray', 'Anhydrite'),
    'bassanite': ('Ca2(SO4)2.H2O', 'evaporite', [], '2CaSO4.(H2O)', 'orthorhombic', None, 2.7, 'the single value listed', None, 'Colorless, Snow white', 'Bassanite'),
    'bischofite': ('MgCl2.6H2O', 'evaporite', [], 'MgCl2.6(H2O)', 'monoclinic', (1.75, 1.5, 2.0), 1.56, 'the single value listed', None, 'Colorless, White', 'Bischofite'),
    'burkeite': ('Na6(CO3)(SO4)2', 'carbonate', [], '', 'orthorhombic', (3.5, 3.5, 3.5), 2.57, 'the single value listed', None, 'Gray, Gray white, White, Light brown, Pink', 'Burkeite'),
    'carnallite': ('KMgCl3.6H2O', 'evaporite', [], 'KMgCl3.6(H2O)', 'orthorhombic', (2.5, 2.5, 2.5), 1.6, 'the single value listed', None, 'Blue, Colorless, Yellow, White, Red', 'Carnallite'),
    'epsomite': ('MgSO4.7H2O', 'evaporite', [], 'MgSO4.7(H2O)', 'orthorhombic', (2.25, 2.0, 2.5), 1.67, 'the listed average of 1.67 to 1.68', (1.67, 1.68), 'Colorless, White, Yellowish white, Greenish white, Pinkish white', 'Epsomite'),
    'glauberite': ('Na2Ca(SO4)2', 'evaporite', [], '', 'monoclinic', (2.75, 2.5, 3.0), 2.77, 'the listed average of 2.7 to 2.85', (2.7, 2.85), 'Colorless, Gray, Gray white, Yellow, Yellowish white', 'Glauberite'),
    'hexahydrite': ('MgSO4.6H2O', 'evaporite', [], 'MgSO4.6(H2O)', 'monoclinic', (2.25, 2.0, 2.5), 1.76, 'the single value listed', None, 'Colorless, White, Greenish white', 'Hexahydrite'),
    'jarosite': ('KFe3(SO4)2(OH)6', 'evaporite', [], 'KFe+++3(SO4)2(OH)6', 'trigonal', (3.0, 2.5, 3.5), 3.09, 'the listed average of 2.9 to 3.3', (2.9, 3.3), 'Brown, Yellow, Yellow brown, Light yellow', 'Jarosite'),
    'kainite': ('KMgSO4Cl.3H2O', 'evaporite', [], 'MgSO4.KCl.3(H2O)', 'monoclinic', (3.0, 3.0, 3.0), 2.1, 'the single value listed', None, 'Blue, Colorless, Violet, Gray, Gray blue', 'Kainite'),
    'kieserite': ('MgSO4.H2O', 'evaporite', [], 'MgSO4.(H2O)', 'monoclinic', (3.5, 3.5, 3.5), 2.57, 'the single value listed', None, 'Colorless, Green yellow, Gray white, Yellowish, White', 'Kieserite'),
    'langbeinite': ('K2Mg2(SO4)3', 'evaporite', [], '', 'cubic', (3.75, 3.5, 4.0), 2.83, 'the single value listed', None, 'Colorless, Violet, Greenish, Gray, Yellowish', 'Langbeinite'),
    'mirabilite': ('Na2SO4.10H2O', 'evaporite', [], 'Na2SO4.10(H2O)', 'monoclinic', (1.75, 1.5, 2.0), 1.49, 'the single value listed', None, 'Colorless, White, Yellowish white, Greenish white', 'Mirabilite'),
    'nahcolite': ('NaHCO3', 'carbonate', [], '', 'monoclinic', (2.5, 2.5, 2.5), 2.21, 'the single value listed', None, 'Colorless, White, Gray, Brown, Black', 'Nahcolite'),
    'natrite': ('Na2CO3', 'carbonate', [], '', 'monoclinic', (3.5, 3.5, 3.5), 2.54, 'the single value listed', None, 'Colorless, Yellow orange, White, Pink', 'Natrite'),
    'polyhalite': ('K2Ca2Mg(SO4)4.2H2O', 'evaporite', [], 'K2Ca2Mg(SO4)4.2(H2O)', 'triclinic', (3.0, 2.5, 3.5), 2.77, 'the listed average of 2.77 to 2.78', (2.77, 2.78), 'Colorless, Gray, Gray white, White, Pink', 'Polyhalite'),
    'sylvite': ('KCl', 'evaporite', [], '', 'cubic', (2.5, 2.5, 2.5), 1.99, 'the single value listed', None, 'White, Yellowish white, Reddish white, Bluish white, Brownish white', 'Sylvite'),
    'thenardite': ('Na2SO4', 'evaporite', [], '', 'orthorhombic', (2.5, 2.5, 2.5), 2.68, 'the listed average of 2.67 to 2.7', (2.67, 2.7), 'White, Grayish white, Yellowish white, Reddish white, Brownish white', 'Thenardite'),
    'thermonatrite': ('Na2CO3.H2O', 'carbonate', [], 'Na2CO3.(H2O)', 'orthorhombic', (1.0, 1.0, 1.0), 2.255, 'the single value listed', None, 'Colorless, Grayish yellow, Grayish, Yellowish, White', 'Thermonatrite'),
    'trona': ('Na3H(CO3)2.2H2O', 'carbonate', [], 'Na3(CO3)(HCO3).2(H2O)', 'monoclinic', (2.5, 2.5, 2.5), 2.13, 'the listed average of 2.11 to 2.17', (2.11, 2.17), 'Colorless, Gray white, Grayish, Yellowish, Yellowish white', 'Trona'),
    'chlorapatite': ('Ca5(PO4)3Cl', 'phosphate', [], '', 'hexagonal', (5.0, 5.0, 5.0), 3.15, 'the listed average of 3.1 to 3.2', (3.1, 3.2), 'White, Light yellow, Pink white', 'Apatite-(CaCl)'),
    'hydroxylapatite': ('Ca5(PO4)3(OH)', 'phosphate', [], '', 'hexagonal', (5.0, 5.0, 5.0), 3.08, 'the single value listed', None, 'Colorless, White, Gray, Yellow, Yellowish green', 'Apatite-(CaOH)'),
    'merrillite': ('Ca18Na2Mg2(PO4)14', 'phosphate', [], '', 'trigonal', None, 3.1, 'the specific gravity listed (the page also gives a bulk electron density of 3.07)', None, '', 'Merrillite'),
    'whitlockite': ('Ca9Mg(PO4)6(PO3OH)', 'phosphate', [], 'Ca9(Mg,Fe++)(PO4)6(PO3OH)', 'trigonal', (5.0, 5.0, 5.0), 3.13, 'the single value listed', None, 'Colorless, Gray, Gray white, Yellowish, White', 'Whitlockite'),
    'mg_perchlorate': ('Mg(ClO4)2', 'evaporite', ['magnesium perchlorate'], '', '', None, 2.21, 'PubChem density annotation 2.21 at 68 F (a US Coast Guard record); the research file had a recalled 2.6', None, '', 'PubChem CID 24840'),
    'ca_perchlorate': ('Ca(ClO4)2', 'evaporite', ['calcium perchlorate'], '', '', None, 2.651, 'PubChem record text: density 2.651 g/cm3', None, '', 'PubChem CID 61629'),
}

# The apatite family (fluorapatite is a GEN.427 row): group id -> (name, member ids).
PART3_GROUPS = {
    "apatite": ("apatite", ("fluorapatite", "chlorapatite", "hydroxylapatite")),
}

# Rows no draw reads (TEST.685, tests/test_materials_reachability.py). A row is reached when a module that uses the table names it (its id, name or an alias),
# or a reached group lists it. A row that is neither is listed here with why, so it is vocabulary a later item may wire, not a row someone forgot; a row listed
# here that a draw reads, or a row that is not in the table, fails the test. Reasons are plain words.
_HP_SILICATE = "high-pressure silicate phase: no draw reads it, the mantle density comes from a fixed table in physics/bulk_structure.py"
_ICE_PHASE = "ice phase kept for the phase diagram: the water layer reads only ice VI, VII and X (physics/bulk_structure.py)"
_SPECIES_ROW = "species row added from the comet and giant-planet inventories (GEN.1359): no draw names it yet"
VOCABULARY_ONLY = {
    "mg_wadsleyite": _HP_SILICATE, "fe_wadsleyite": _HP_SILICATE, "mg_ringwoodite": _HP_SILICATE, "fe_ringwoodite": _HP_SILICATE,
    "mg_bridgmanite": _HP_SILICATE, "fe_bridgmanite": _HP_SILICATE, "mgfe_bridgmanite": _HP_SILICATE, "mg_post_perovskite": _HP_SILICATE,
    "fe_post_perovskite": _HP_SILICATE, "ca_perovskite": _HP_SILICATE, "al_perovskite": _HP_SILICATE, "al_post_perovskite": _HP_SILICATE,
    "mg_akimotoite": _HP_SILICATE, "mg_majorite": _HP_SILICATE, "hp_clinoenstatite": _HP_SILICATE, "clinoenstatite": _HP_SILICATE,
    "coesite": _HP_SILICATE, "stishovite": _HP_SILICATE, "seifertite": _HP_SILICATE,
    "ice_ii": "ice phase kept for the phase diagram: the water layer reads only ice VI, VII and X (physics/bulk_structure.py)",
    "ice_iii": "ice phase kept for the phase diagram: the water layer reads only ice VI, VII and X (physics/bulk_structure.py)",
    "ice_v": "ice phase kept for the phase diagram: the water layer reads only ice VI, VII and X (physics/bulk_structure.py)",
    "hydrogen": "envelope row: the giant-planet draws read the gas table (physics/gases), not this row",
    "helium": "envelope row: the giant-planet draws read the gas table (physics/gases), not this row",
    "potassium": "element row: the element draws go through physics/elements, not this row",
    "tungsten": "element row: the element draws go through physics/elements, not this row",
    "fesi": "light-element compound of a core: no draw names it",
    "carbon_diamond": "carbon polymorph: the carbon draws name graphite, black and amorphous carbon only",
    "amorphous_carbon": "carbon form: no draw names it",
    "wustite": "iron oxide: no draw names it",
    "goethite": "iron hydroxide: only the band albedo of the row exists, no recipe names it",
    "komatiite": "rock type: no recipe or taxon names it",
    "hematite_lt10um": "a grain size of hematite (a double of the hematite row): no draw names it",
    "hydrocarbons": "class text of the planet classes (legacy_names.class_compositions): no draw reads the mapping",
    "frozen_hydrocarbons": "class text of the planet classes (legacy_names.class_compositions): no draw reads the mapping",
    # GEN.1359: the species rows the inventories name; a later item wires each to its draw
    "ice_ic": _ICE_PHASE, "ice_iv": _ICE_PHASE, "ice_ix": _ICE_PHASE, "ice_xi": _ICE_PHASE, "ice_xv": _ICE_PHASE,
    "o2_ice": _SPECIES_ROW, "xe_ice": _SPECIES_ROW, "ar_ice": _SPECIES_ROW, "kr_ice": _SPECIES_ROW, "c2h2_ice": _SPECIES_ROW,
    "nh4cl": _SPECIES_ROW, "nh4hco3": _SPECIES_ROW, "nh4sh": _SPECIES_ROW, "ammonium_formate": _SPECIES_ROW, "ammonium_sulfate": _SPECIES_ROW,
    "ammonia_monohydrate": _SPECIES_ROW, "sulfuric_acid_monohydrate": _SPECIES_ROW, "sulfuric_acid_tetrahydrate": _SPECIES_ROW,
    "tholin": "organic solid of the haze: physics/aerosols.py carries its own density constant and no draw reads this row yet",
    "insoluble_organic_matter": "organic solid of chondrites and comets: the comet draw reads the lumped complex_organics row, not this one",
    **{row: "added by GEN.1356; the recipe swap of GEN.1384 reads it (the reduced-body taxa)" for row in (
        "pyrrhotite", "pentlandite", "chalcopyrite", "oldhamite", "niningerite", "alabandite", "daubreelite", "osbornite",
        "sinoite", "khamrabaevite")},
}

# The rows GEN.1355 adds as data. The norm (physics/rock_types.py) names nepheline, leucite, hypersthene and wollastonite, so those four are reached; the rest no
# recipe or taxon names yet. They are tagged here, as the item says the tags come with the rows, and lose the tag when a draw reads them.
_NORM_NAMED = ("nepheline", "leucite", "hypersthene", "wollastonite")
_ADDED_AS_DATA = "mineral added as table data (GEN.1355): no recipe or taxon names it yet"
VOCABULARY_ONLY.update({row_id: _ADDED_AS_DATA for row_id in (*PART1_MINERALS, *PART1_GROUPS) if row_id not in _NORM_NAMED})
VOCABULARY_ONLY["smectite"] = ("group added as table data (GEN.1355): surface_composition.py names Hu's Mg-smectite only in a comment, kaolinite stands in, and swapping "
                               "it is a change to drawn output (a later epoch item)")

# GEN.1357 rows: no recipe or taxon names them yet (the evaporite placeholder is halite, gypsum and calcite; GEN.1366 draws the salts).
_ADDED_BY_1357 = "mineral added as table data (GEN.1357): no recipe or taxon names it yet; the evaporite draw of GEN.1366 reads the salts"
_NORM_NAMED_1357 = ("apatite", "chlorapatite", "hydroxylapatite")   # the norm (physics/rock_types.py) names apatite, and the group reaches its members
VOCABULARY_ONLY.update({row_id: _ADDED_BY_1357 for row_id in (*PART3_MINERALS, *PART3_GROUPS) if row_id not in _NORM_NAMED_1357})
