> **Imported research, kept as received (banner added 2026-10-15; the text below is unchanged).** It is Boss's document, not a statement of what planetGen does. The design notes that use it: [activity-magnetism-radiation-hydrosphere](../../habitability/activity-magnetism-radiation-hydrosphere.md), [habitability-index](../../habitability/habitability-index.md), [units-and-number-formatting](../../architecture/units-and-number-formatting.md).

Chemical Habitability, Geochemical Toxicity, and Environmental Margins of Terrestrial Exoplanets
================================================================================================

Geochemical Determinants of Terrestrial and Aquatic Toxicity
------------------------------------------------------------

The habitability of a terrestrial exoplanet is conventionally evaluated through macroscopic astrophysical parameters: stellar irradiance, planetary mass, orbital eccentricity, and the presence of liquid water within the circumstellar habitable zone [R0843, R3986]. However, the geochemical inventory and surface thermodynamic state determine whether a planetary surface is biologically viable, biochemically inert, or actively cytotoxic [R0843]. Toxicity in planetary environments is defined by the disruption of macromolecular integrity, the inhibition of energetic metabolism, or the destruction of cellular membranes in carbon-based life [R1763]. These toxic states emerge through specific chemical species, extreme ionic conditions, and thermodynamic constraints on water availability [R0325].

Heavy metals and metalloids—most prominently arsenic, cadmium, lead, mercury, and hexavalent chromium—act as potent cellular poisons across phylogenetic domains [R0843]. Their toxicity stems from their behavior as soft or intermediate Lewis acids, showing high affinities for sulfhydryl ($-SH$) groups in cysteine residues [R0843]. This affinity drives the non-functional displacement of essential transition metal cofactors, such as ferrous iron ($\text{Fe}^{2+}$), zinc ($\text{Zn}^{2+}$), and manganese ($\text{Mn}^{2+}$), in catalytic metalloenzymes, which permanently arrests core metabolic pathways.

Oxychlorine species, including perchlorates ($\text{ClO}_4^-$), chlorates ($\text{ClO}_3^-$), and hypochlorites ($\text{ClO}^-$), introduce powerful oxidative and metabolic hazards ([R0690]; NASA, 2024). While kinetically metastable in dilute aqueous solutions at ambient temperatures, perchlorates become reactive oxidants when exposed to ionizing radiation or thermal activation [R0690]. Physiologically, oxychlorines compete with iodide ($\text{I}^-$) at the sodium-iodide symporter in vertebrate systems, disrupting thyroid hormonogenesis, basal metabolic regulation, and systemic organ performance [R0690, R1019]. In bacterial and archaeal systems lacking specialized dissimilatory perchlorate reduction pathways, oxychlorine compounds promote severe intracellular oxidative stress, mutating nucleic acids and degrading cell membranes (Curr Trends Biol Sci, 2025; [R4580]).

Reactive oxygen species, such as hydrogen peroxide ($\text{H}_2\text{O}_2$), the superoxide radical ($\text{O}_2^{\bullet-}$), and hydroxyl radicals ($\text{OH}^\bullet$), represent unavoidable intermediates of photolytic and radiolytic processes in volatilized water and oxidized regoliths [R0690]. Hydroxyl radicals possess oxidation potentials capable of non-selectively abstracting hydrogen from organic backbones, inducing fatal cross-linking in proteins and double-strand breaks in nucleic polymers [R0690]. Volatile chalcogens, including hydrogen sulfide ($\text{H}_2\text{S}$) and sulfur dioxide ($\text{SO}_2$), act as double-edged environmental drivers; while capable of supporting specialized chemolithotrophic respiration, elevated partial pressures in anoxic atmospheres irreversibly arrest aerobic respiration by outcompeting oxygen at cytochrome $c$ oxidase [R0631, R1488].

Beyond isolated chemical poisons, physical chemistry defines global habitability through water activity ($a_w$) and solution chaotropicity [R0242, R0915]. Water activity, defined as the ratio of the vapor pressure of water in a solution to the vapor pressure of pure water ($a_w = \frac{p}{p_0}$) at constant temperature, governs thermodynamic water availability [R4222]. Terrestrial microbial life is bounded by an absolute lower threshold for cell division at $a_w \approx 0.605 - 0.611$, observed in extremely xerophilic fungi and halophilic archaea [R4222, R4396]. Below this limit, intracellular hydration shells collapse, halting enzyme activity and osmotic turgor maintenance [R1763, R4222].

Chaotropic solutes, characterized by large, low-charge-density ions such as magnesium ($\text{Mg}^{2+}$), calcium ($\text{Ca}^{2+}$), perchlorate ($\text{ClO}_4^-$), and iodide ($\text{I}^-$), weaken electrostatic and hydrogen-bonding networks in aqueous matrices [R0242]. By disrupting water ordering, chaotropes lower enthalpy barriers for denaturation, unfolding functional protein structures and fluidizing lipid bilayers to the point of lytic collapse [R0242, R1763]. In contrast, kosmotropic solutes, such as sulfate ($\text{SO}_4^{2-}$) and phosphate ($\text{HPO}_4^{2-}$), stabilize water structures and preserve macromolecular configurations, though excessive concentrations restrict the conformational flexibility required for enzymatic turnover [R1763]. The balance between water activity, ionic strength, and chaotropicity frequently creates uninhabitable environments within bodies of liquid water, demonstrating that liquid stability alone is insufficient for biological viability [R0325, R0843].

| **Chemical or Thermodynamic Stressor**        | **Primary Environmental Phase**                                                               | **Principal Abiotic Source / Driving Mechanism**                                                  | **Specific Biochemical / Physiological Toxicity Mechanism**                                                                                        | **Universal Biological Limit**                                                                                         | **Human Physiological Threshold**                                                          |
| --------------------------------------------- | --------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| Perchlorate (ClO4-) / Chlorate (ClO3-)        | Surface regolith, eutectic brines [R0690]                                      | Photochemical oxidation of chlorine, surface cosmic ray radiolysis [R0690]         | Oxidative denaturation under UV; competitive inhibition of sodium-iodide symporter [R1019, R4580]                | Bactericidal under solar UV; tolerated by select perchlorate reducers [R1019, R4580] | 24.5 ug/L chronic drinking water limit [R1019]                               |
| Reactive Oxygen Species (H2O2, OH*)           | Regolith crust, atmospheric boundary layer [R0690]                             | Water photolysis, mineral triboelectric fracturing, UV ionization [R0690]          | Non-selective free radical oxidation of lipids, protein backbone scission, DNA degradation [R0690]                                  | Bounded by enzymatic antioxidant capacity (catalase, SOD) [R0690]                                       | Severe mucosal necrosis and respiratory failure at >1 ppm vapor                            |
| Chaotropic Cations (Mg2+, Ca2+)               | Hypersaline exo-oceans, eutectic ground brines [R0242, R0915]   | Evaporative concentration, basalt/ultramafic leaching [R4396]     | Weakening of hydrophobic interactions; structural collapse of lipid bilayers and protein unfolding [R0325, R1763] | Absolute boundary observed at chaotropic activities >73.8 kJ/kg [R0242, R0915]           | Systemic osmotic dehydration, electrolyte dysregulation                                    |
| Depressed Water Activity (aw)                 | Dense saline water bodies, cryogenic brines [R1763, R4222] | High solute saturation (NaCl, CaCl2, sulfates) [R1763, R4222]  | Loss of macromolecular hydration; inability to maintain cellular turgor; metabolic arrest [R1763, R4222]        | Cell division ceases globally below aw = 0.605 - 0.611 [R4222]                                        | Lethal through rapid tissue desiccation upon ingestion                                     |
| Carbon Monoxide (CO)                          | Lower atmosphere, shallow marine photic zone [R3986]                      | Photolysis of CO2 with deficient OH sink around M-dwarfs [R3986]              | High-affinity coordination to iron porphyrin centers; irreversible carboxyhemoglobin complexation [R3986]                      | Metabolically utilized by carboxydotrophs; highly toxic to aerobic respiration [R3986]             | Lethal at sustained atmospheric partial pressures >100-200 ppm [R3986, R3987] |
| Sulfur Dioxide (SO2) / Hydrogen Sulfide (H2S) | Volcanic atmospheres, euxinic oceans [R1488]                                  | Mantle outgassing dictated by oxygen fugacity and low degassing pressures [R1488] | Hydrolytic formation of sulfurous acids; irreversible inhibition of cytochrome c oxidase [R1488]                                   | Chemolithoautotrophic substrate at low levels; cytotoxic at elevated concentrations [R0631]        | Acute systemic pulmonary toxicity and cellular asphyxiation at >50 ppm                     |

Endogenic Geochemical and Abiotic Drivers of Toxic Solutes
----------------------------------------------------------

The generation, transport, and concentration of toxic chemical suites on a terrestrial exoplanet are governed by fundamental endogenic processes: mantle thermodynamics, magmatic outgassing, tectonic regimes, hydrothermal circulation, and chemical weathering [R0843, R1488].

The composition of a planet's atmospheric and oceanic volatile inventory is coupled to the mantle oxidation state, quantified as oxygen fugacity ($f_{\text{O}_2}$) relative to standard mineral buffers such as Fayalite-Magnetite-Quartz (FMQ) or Iron-Wüstite (IW) [R1488]. On planets with oxidized upper mantles ($f_{\text{O}_2} \ge \text{FMQ}$), volcanic degassing releases sulfur primarily as sulfur dioxide ($\text{SO}_2$) and carbon as carbon dioxide ($\text{CO}_2$) [R1488]. Dissolution of subaerial and submarine $\text{SO}_2$ plumes into surface waters triggers disproportionation and hydrolysis:

$4\text{SO}_2 + 4\text{H}_2\text{O} \to 3\text{H}_2\text{SO}_4 + \text{H}_2\text{S}$

This reaction delivers massive sulfuric acid fluxes to surface reservoirs, overwhelming local buffering systems and generating hyperacidic marine environments [R1488]. Conversely, reduced mantles ($f_{\text{O}_2} \le \text{IW}$) outgas sulfur predominantly as hydrogen sulfide ($\text{H}_2\text{S}$) and carbon as carbon monoxide ($\text{CO}$) or methane ($\text{CH}_4$) [R0843, R1488]. Under high submarine lithostatic pressures, volatile solubility models show that $\text{H}_2\text{S}$ is preferentially partitioned into hydrothermal fluids [R1488]. This sustains widespread euxinia in the absence of biological oxidation sinks, enriching water columns with toxic dissolved sulfides [R1488].

Hydrothermal circulation across ultramafic lithospheres drives serpentinization, which dictates both solution pH and redox conditions [R0843]. When oceanic fluids interact with olivine-rich peridotites at temperatures between $150^\circ\text{C}$ and $350^\circ\text{C}$, the oxidation of ferrous iron in the fayalite component by water produces magnetite, molecular hydrogen, and serpentine-group minerals:

$3\text{Fe}_2\text{SiO}_4\text{ (fayalite)} + 2\text{H}_2\text{O} \to 2\text{Fe}_3\text{O}_4\text{ (magnetite)} + 3\text{SiO}_2 + 2\text{H}_2$

$6\text{Mg}_2\text{SiO}_4\text{ (forsterite)} + \text{SiO}_2 + 8\text{H}_2\text{O} \to 2\text{Mg}_3\text{Si}_2\text{O}_5(\text{OH})_4\text{ (serpentine)}$

This process consumes hydronium ions, releasing hydroxide and generating hyperalkaline effluent fluids ($\text{pH } 9-12$) enriched in $\text{H}_2$, dissolved methane via abiotic Sabatier and Fischer-Tropsch-type syntheses, and calcium cations ($\text{Ca}^{2+}$) [R0843]. If serpentinization occurs on a global scale—such as in an ocean basin lacking subaerial granitic continents—marine chemistry transitions toward an alkaline state [R2313]. The extreme pH of these waters denatures standard nucleic acid structures and saponifies biological lipids, presenting severe challenges to cellular integrity despite providing abundant chemical energy [R0843].

Chemical weathering cycles act as secondary geochemical filters on planetary toxicity [R2488]. On worlds where atmospheric $p\text{CO}_2$ is maintained at high levels to counteract faint stellar conditions, carbonic acid leaching accelerates the denudation of basaltic crusts [R0692]:

$\text{CaSiO}_3 + 2\text{CO}_2 + \text{H}_2\text{O} \to \text{Ca}^{2+} + 2\text{HCO}_3^- + \text{SiO}_2$

While this process provides climate-stabilizing negative feedbacks via carbonate precipitation, an imbalance in cation extraction can saturate surface reservoirs with heavy divalent metals [R2488].

In arid conditions with thin atmospheres, aqueous weathering is replaced by photochemical oxidation and ionizing radiolysis [R0690]. Solar ultraviolet flux ($\lambda < 200\text{ nm}$) and galactic cosmic rays penetrate directly to the regolith, breaking surface water vapor down into $\text{OH}^\bullet$ and $\text{HO}_2^\bullet$ radicals [R0690]. These oxidants react with atmospheric and detrital chloride salts to form chlorate ($\text{ClO}_3^-$) and perchlorate ($\text{ClO}_4^-$) through pathways mediated by ozone [R0690]:

$\text{Cl}^- + \text{O}_3 \to \text{ClO}^- + \text{O}_2 \xrightarrow{h\nu\text{, ROS}} \text{ClO}_4^-$

Over geological timescales, this abiotic photochemistry transforms the upper meters of surface regolith into a caustic, oxidizing layer that destroys organic compounds and poses a chronic toxicity hazard across planetary terrains [R0690].
Plausible Chemistries and Architectures of Terrestrial Exo-Oceans

-----------------------------------------------------------------

Exo-ocean chemistry reflects the interaction between initial volatile inventories, mantle-crustal petrology, atmospheric composition, and ocean depth [R2381, R2488]. Deviations from Earth-like baseline conditions give rise to several distinct oceanic archetypes within the ecosphere.

### Hyperacidic Sulfate Oceans

On volatile-rich worlds where mantle outgassing is dominated by sulfur species ($\text{SO}_2$, $\text{H}_2\text{S}$) rather than carbon, or where mafic continental weathering is outpaced by volcanic acid flux, oceans become massive sulfuric acid reservoirs [R1488, R4580]. Characterized by pH values between $1.0$ and $4.0$, these waters are sustained by high concentrations of dissolved bisulfate ($\text{HSO}_4^-$) and sulfate ($\text{SO}_4^{2-}$). The low pH accelerates the dissolution of crustal minerals, mobilizing toxic transition metals such as $\text{Fe}^{3+}$, $\text{Al}^{3+}$, $\text{Cu}^{2+}$, and $\text{Zn}^{2+}$ into solution at concentrations orders of magnitude higher than Earth baselines [R0843, R2626]. In such waters, proton motive forces are disrupted; terrestrial acidophilic archaea rely on inverted membrane potentials to withstand comparable environments, but high heavy-metal concentrations establish a steep evolutionary and physiological hurdle [R0817].

### Alkaline Soda Oceans

When a carbon-rich planetary atmosphere contacts an ultramafic or basaltic lithosphere undergoing widespread serpentinization without subaerial granite exposure, oceans evolve into alkaline soda reservoirs [R0843, R2313]. These waters are characterized by high pH ($9.0 - 11.5$) and elevated concentrations of sodium bicarbonate ($\text{NaHCO}_3$) and sodium carbonate ($\text{Na}_2\text{CO}_3$) [R2313]. This chemistry fundamentally alters nutrient availability, particularly for phosphorus. In neutral to mildly acidic waters containing free $\text{Ca}^{2+}$, orthophosphate precipitates as insoluble apatite ($\text{Ca}_5(\text{PO}_4)_3(\text{OH,F,Cl})$), leading to the prebiotic "phosphate problem" where phosphorus becomes scarce for biological synthesis [R4399].

In a soda ocean, carbonate ions outcompete phosphate for calcium, precipitating calcite or aragonite ($\text{CaCO}_3$) and stripping the water column of free $\text{Ca}^{2+}$ [R4399]:

$\text{Ca}^{2+} + \text{CO}_3^{2-} \to \text{CaCO}_3\downarrow$

This leaves dissolved orthophosphate unconstrained in solution, allowing phosphorus concentrations to rise into the millimolar and molar ranges [R4399]. While this chemical state resolves a major nutrient constraint, its high alkalinity destabilizes nucleotide ester bonds and RNA backbones, requiring specialized cellular envelopes and stable internal pH regulators to avoid alkaline cell lysis.

### Hypersaline Chaotropic Brines

On worlds with low hydrospheric water mass fractions or endorheic continental hydrological basins subject to prolonged arid desiccation, oceans and epicontinental seas concentrate into hypersaline brines [R0325, R4396]. Unlike marine systems dominated by kosmotropic sodium chloride ($\text{NaCl}$), advanced geochemical evaporation pathways precipitate halite early, leaving residual brines dominated by calcium chloride ($\text{CaCl}_2$) and magnesium chloride ($\text{MgCl}_2$)—analogous to terrestrial features like Don Juan Pond and the Dallol geothermal complex [R0325, R4396]. These systems combine low water activities ($a_w \ll 0.60$) with high chaotropicity ($>70\text{ kJ/kg}$) [R0242, R0915]. The divalent $\text{Mg}^{2+}$ and $\text{Ca}^{2+}$ cations disrupt hydrogen bonding in cellular hydration shells, destabilizing proteins and membrane structures and establishing a physicochemical barrier that arrests native biological function [R0242, R0915].

### High-Pressure Ice-Sealed Oceanworlds

Planets possessing deep hydrospheres—where water comprises $>1-10\text{ wt}\%$ of planetary mass—develop oceans hundreds of kilometers deep [R2381]. At abyssal depths, immense hydrostatic pressures ($>1\text{ GPa}$) force liquid water to freeze into high-pressure ice polymorphs (Ice VI, VII, and X) regardless of temperature [R0777, R2381]. This high-pressure ice mantle forms a physical barrier that isolates the liquid ocean from the underlying silicate mantle, shutting down direct water-rock interactions and hydrothermal circulation [R2381, R2626].

Without seafloor silicate weathering or hydrothermal inputs, the ocean is cut off from primary geochemical buffering and mantle phosphorus fluxes [R2381, R2626]. Atmospheric carbon dioxide equilibrates with this isolated water column without being sequestered into mineral carbonates, causing sustained acidification down to $\text{pH } 3-5$ [R2381, R2626]. These oceans become dilute, nutrient-poor reservoirs characterized by thermodynamic disequilibrium, where native abiogenesis is limited by severe macronutrient starvation [R2381].

| **Exo-Ocean Classification** | **Dominant Dissolved Solutes**                                           | **Representative pH Range**           | **Typical Water Activity (aw)**  | **Key Biogeochemical Dynamics**                                                                                                  | **Primary Bottlenecks for Habitability**                                                                                      |
| ---------------------------- | ------------------------------------------------------------------------ | ------------------------------------- | -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Acidic Sulfate-Rich          | HSO4-, SO4(2-), Fe3+, Mg2+ [R4580]                   | 1.0 - 4.5 [R4580] | 0.90 - 0.98                      | High crustal metal leaching; outgassed sulfur dominates over carbon buffering [R1488, R4580] | Metal toxicity; proton motive force disruption; extreme acid stress [R0843]             |
| Alkaline Soda Ocean          | Na+, HCO3-, CO3(2-), Cl- [R2313]    | 9.0 - 11.5 [R0690]     | 0.92 - 0.99                      | Apatite inhibition via CaCO3 precipitation; high phosphate enrichment [R4399]                                    | RNA phosphodiester bond hydrolysis; cytoplasmic alkalinization [R0690, R2381]                      |
| Chaotropic Chloride Brine    | Ca2+, Mg2+, Cl-, ClO4- [R1763, R4222] | 3.5 - 6.5 [R0325]      | 0.30 - 0.65 [R4396] | Halite precipitation concentrates divalent cations; low freezing depression [R1763, R4222]    | Protein unfolding; membrane fluidization; severe thermodynamic water restriction [R0242, R0915] |
| Ice-Sealed Waterworld        | H2CO3, HCO3-, trace salts [R2381]                            | 3.0 - 5.5 [R2626]        | >0.99                            | High-pressure ice mantle isolates water from the core; closed-system carbon storage [R2381]                          | Severe phosphorus and transition-metal starvation; unbuffered carbonic acidification [R2381]                      |

Atmospheric Pressure and Photochemical Regimes in Marginal Habitability
-----------------------------------------------------------------------

Atmospheric column depth and mass dictate the surface radiation environment, gas retention, and chemical equilibrium of surface waters, creating distinctive toxicity regimes in marginal atmospheres ([R0690]; [R3986, R3987]).

On planets with thin atmospheres ($P_{\text{surf}} \le 0.1\text{ bar}$), low column depths allow stellar ultraviolet (FUV and NUV) radiation and galactic cosmic rays to penetrate unhindered to the surface [R0690, R4580]. This high radiation flux drives the photolysis of volatile water and carbon dioxide, producing reactive oxidants without supporting an ozone protective layer [R0690]. The resulting surface environment features dry deposition of superoxides, hydrogen peroxide, and oxychlorines, producing a reactive chemical blanket across topsoils [R0690]. Any standing water remains vulnerable to rapid atmospheric escape via thermal jeans and non-thermal hydrodynamic escape, raising residual salinities and concentrating toxic mineral components [R2678, R4774]. While unshielded microbial life is rapidly sterilized by combined UV and chemical oxidation, technological survival remains viable through pressurized shielding, utilizing ambient oxychlorines as consumable oxidizers and oxygen reservoirs [R1019, R4580].

Conversely, thick atmospheres ($P_{\text{surf}} \ge 2 - 10\text{ bar}$) generate distinct toxicity profiles driven by gas partial pressures and photochemical stagnation [R3986]. To maintain surface temperatures near the outer boundary of the habitable zone, climate models require high carbon dioxide partial pressures ($p\text{CO}_2 \ge 1 - 5\text{ bar}$) [R3986]. While simple unicellular organisms can adapt to high dissolved inorganic carbon across diverse pH ranges, these atmospheres severely limit complex aerobic life [R3986, R3987]. For mammalian physiology, carbon dioxide levels exceeding $p\text{CO}_2 \approx 0.005 - 0.01\text{ bar}$ trigger hypercapnia, systemic acidosis, and neurological dysfunction, with lethal thresholds around $0.05\text{ bar}$ [R3986]. This physiological constraint restricts the Habitable Zone for Complex Life (HZCL) to roughly $20-30\%$ of the traditional microbial habitable zone [R3986].

Stellar spectral irradiance strongly influences atmospheric chemical balance [R3986]. Terrestrial planets orbiting low-mass M-dwarf stars receive a stellar flux with an elevated ratio of far-ultraviolet to near-ultraviolet radiation (FUV/NUV) relative to Sun-like G-dwarfs [R3986, R3987]. Near-ultraviolet photons ($\lambda \approx 290 - 330\text{ nm}$) are required to photolyze tropospheric ozone ($\text{O}_3$), producing excited oxygen atoms ($\text{O}(^1\text{D})$) that react with water vapor to generate hydroxyl radicals ($\text{OH}^\bullet$):

$\text{O}_3 + h\nu \to \text{O}_2 + \text{O}(^1\text{D})$

$\text{O}(^1\text{D}) + \text{H}_2\text{O} \to 2\text{OH}^\bullet$

Because M-dwarfs emit low NUV fluxes relative to their FUV output, tropospheric $\text{OH}^\bullet$ concentrations are reduced [R3986, R3987]. The reaction with hydroxyl radicals represents the primary atmospheric sink for carbon monoxide:

$\text{CO} + \text{OH}^\bullet \to \text{CO}_2 + \text{H}^\bullet$

Without adequate hydroxyl production, this sink fails [R3986, R3987]. As FUV photolysis continues to split atmospheric $\text{CO}_2$ into $\text{CO}$ and atomic oxygen, carbon monoxide accumulates abiotically, reaching concentrations between several hundred parts per million and a few percent of total atmospheric volume [R3986]. While primitive chemotrophs can exploit $\text{CO}$ as an energetic substrate, these atmospheres are toxic to complex aerobic organisms [R3986]. In mammalian physiology, carbon monoxide binds to hemoglobin with an affinity approximately 200 times higher than that of oxygen, blocking cellular respiration and leading to asphyxiation even in oxygenated environments [R3986, R3987].

## Thermodynamic Extremes: Eutectic Cryo-Brines and Moist-Greenhouse Leaching

Planetary habitability is bounded by thermal extremes that transform surface water from a balanced solvent into a hostile, toxic medium [R1763, R4396].

At the sub-freezing cryo-margin, mean planetary surface temperatures fall well below $273\text{ K}$ [R0843, R4396]. Under these conditions, pure liquid water cannot persist, and aqueous environments are restricted to concentrated eutectic brines [R4396]. Different salt mixtures depress the freezing point of water to varying extents [R4396]. For example, the eutectic temperature ($T_{\text{eut}}$) drops to $-21.1^\circ\text{C}$ for $\text{NaCl}-\text{H}_2\text{O}$, $-36^\circ\text{C}$ for $\text{FeCl}_2-\text{H}_2\text{O}$, $-50.4^\circ\text{C}$ for $\text{Mg(ClO}_4)_2-\text{H}_2\text{O}$, and $-70^\circ\text{C}$ for $\text{Ca(ClO}_4)_2-\text{H}_2\text{O}$ [R2880, R4396].

While eutectic brines allow liquid solutions to remain stable in cold environments, their physical chemistry introduces severe biological stresses [R4396]. Freezing depression requires high solute concentrations, which sharply lowers water activity [R1763, R4222]. A perchlorate or chloride brine near its eutectic point exhibits water activity levels well below the terrestrial biological minimum ($a_w \approx 0.20 - 0.55$) [R4222]. Concurrently, the chaotropic effects of concentrated divalent cations ($\text{Mg}^{2+}, \text{Ca}^{2+}$) destabilize intracellular structures [R0242, R0915]. Unshielded cells experience rapid dehydration and structural degradation, rendering cold eutectic brines chemically uninhabitable despite their liquid state [R0325, R4222].

At the high-temperature margin, increasing stellar irradiance drives a planet toward the moist-greenhouse and runaway-greenhouse thresholds [R2085, R3236]. The runaway greenhouse state is governed by the Komabayashi-Ingersoll and Simpson-Nakajima limits [R1608, R2085]. These establish an upper limit on the outgoing longwave radiation (OLR) that a water-saturated troposphere can emit to space:

$\text{OLR}_{\text{limit}} \approx 280 - 310\text{ W m}^{-2}$

When absorbed stellar flux exceeds this ceiling, the planet cannot radiate away incoming energy, and surface temperatures rise continuously [R1608, R2085, R3236]. Before entering a full runaway state, a planet passes through a moist greenhouse phase [R0843, R2275]. In this regime, surface temperatures range between $50^\circ\text{C}$ and $90^\circ\text{C}$, causing water vapor to saturate the troposphere and bypass the cold-trap into the stratosphere [R0843, R2275]. There, stellar UV radiation photolyzes water molecules, driving rapid hydrogen loss to space and gradual planetary desiccation [R0843, R2275].

Geochemically, this elevated thermal state intensifies chemical weathering and mineral leaching [R2626]. Kinetic weathering rates follow the Arrhenius relationship:

$k = A \exp\left(-\frac{E_a}{R T}\right)$

Elevated temperatures accelerate the extraction of silicate and trace elements from continental and seafloor rocks [R0692, R2381]. Intense evaporation strips freshwater systems, driving residual oceans toward saturation in heavy metals, silica ($\text{SiO}_2SiO2\text{SiO}_2SiO2​SiO2SiO2SiO2SiO2​\text{SiO}_2SiO2\text{SiO}_2SiO2​SiO2​SiO2SiO2​SiO2​$), borates, and sulfates [R2626]. The combination of thermal stress, high mineral loads, and elevated salinity imposes severe evolutionary constraints, precluding complex multi-tissue organisms and restricting any indigenous ecology to hyperthermophilic single-celled life [R0843, R1763].

Astrobiological Divergence: Native Abiogenesis Versus Technological Habitation
------------------------------------------------------------------------------

The environmental parameters required for the origin of life (_abiogenesis_) differ substantially from those necessary for the survival of complex organisms, and both differ from the requirements for technological human survival [R0843, R3986]. A planet chemically hostile to the emergence of indigenous biology can nonetheless support long-term technological habitation [R3986, R3987].

Abiogenesis requires specific environmental conditions to support the assembly of primitive metabolic cycles and replicating informational polymers [R0843]. These include concentrated organic precursors—such as hydrogen cyanide ($\text{HCN}$), formaldehyde ($\text{HCHO}$), and formamide—maintained in localized aqueous systems [R4399]. Prebiotic chemistry further demands accessible macronutrients, particularly soluble phosphate not locked in insoluble apatite minerals [R4399]. Thermodynamic driving forces, such as wet-dry cycles on exposed landmasses, are essential to facilitate condensation reactions and polymerize peptides and oligonucleotides [R0843, R0989]. Finally, moderate pH and redox gradients are necessary to establish primitive proton motive forces without degrading exposed catalytic structures [R0690].

Many planetary environments fail these criteria completely [R0843]. In high-pressure ice-sealed waterworlds, the absence of subaerial land eliminates the wet-dry cycles required for nucleotide polymerization, while the ice barrier shuts off access to crustal catalytic surfaces and mineral phosphorus [R2381]. On thin-atmosphere worlds, surface oxidizing agents destroy organic precursors before self-assembling prebiotic structures can mature [R0690]. In hyperalkaline soda systems, the same conditions that concentrate phosphorus also accelerate the alkaline hydrolysis of RNA phosphodiester backbones, preventing genetic replication without protective mechanisms [R4399]. Consequently, an exoplanet may remain completely sterile throughout its evolutionary history due to these prebiotic bottlenecks, while remaining physically located within the circumstellar habitable zone ([R0843]; [R3986, R3987]).

For technological human habitation, these prebiotic geochemical constraints are largely irrelevant. A technological civilization relies on enclosed systems and engineered life-support architectures, utilizing local raw materials through in-situ resource utilization (ISRU) rather than depending on native biospheric cycles (Curr Trends Biol Sci, 2025; [R1019]). Environmental toxins that prevent native biogenesis can frequently be processed into useful resources [R1019].

Perchlorates ($\text{ClO}_4^-$), which are toxic to both plants and human thyroid function at low concentrations, can be treated using perchlorate-reducing bacteria or synthetic molybdenum catalysts (Curr Trends Biol Sci, 2025; [R1019]):

$\text{ClO}_4^- \xrightarrow{\text{perchlorate reductase}} \text{ClO}_2^- \xrightarrow{\text{chlorite dismutase}} \text{Cl}^- + \text{O}_2$

This bioprocess detoxifies local soils and water supplies while generating elemental oxygen for life support and propellant systems [R1019]. Chaotropic and hypersaline waters can be processed using multistage nanofiltration and high-pressure reverse osmosis, stripping out divalent cations and restoring water activity to potable levels while generating mineral salts for construction and metallurgical extraction [R1019, R1763]. Atmospheres containing elevated levels of carbon dioxide, carbon monoxide, and sulfur dioxide can be managed within enclosed habitats using regenerative chemical scrubbing and catalytic Sabatier reactors, which reduce carbon oxides into methane and potable water ([R1019]; [R3986, R3987]). As a result, a planet classified as barren or hostile to primitive life can remain fully viable for sustained technological settlement, demonstrating that biological habitability and technological habitability represent distinct regimes in planetary science ([R0843]; [R3986, R3987]).

| **Planetary Regime Archetype** | **Atmospheric State & Pressure**                                                                              | **Marine / Aqueous Chemistry**                                                                                 | **Dominant Geochemical Toxicity Drivers**                                                                         | **Native Abiogenesis Potential**                                                                                          | **Technological Habitation Outlook**                                                                                      |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Cryogenic Perchlorate World    | Thin (P <= 0.1 bar); CO2-dominated; unshielded UV/GCR [R0690, R4580]       | Eutectic brines; subzero liquid state; depressed aw << 0.60 [R4222, R4396]       | Oxychlorines (ClO4-), free radicals, high chaotropicity [R0690]                 | Prevented: prebiotic organics oxidized; sub-threshold water activity [R4222]                 | High: surface water extraction viable; perchlorates reduced to metabolic O2 [R1019]                         |
| Super-Earth Soda World         | Moderate (P ~ 1 - 3 bar); neutral to reducing outgassing [R1488]               | Hyperalkaline (pH 9.5 - 11.5); Na2CO3-rich; high dissolved P [R0690, R4399]     | Severe alkalinity; rapid RNA phosphodiester bond hydrolysis [R0690, R2381]             | Moderate: high phosphorus availability, but requires alkaline stabilization [R0690, R4399] | High: abundant accessible water; simple acid-base neutralization required                                                 |
| Outer-Edge Hypercapnic World   | Thick (P >= 2 - 5 bar); elevated pCO2; dense cloud cover ([R0692]; [R3986, R3987]) | Mildly acidic carbonic waters; high mineral bicarbonate [R0692]                              | Severe systemic hypercapnia; blood acidosis in complex life [R3986, R3987]                           | High: stable clement temperatures, active silicate-carbonate cycles [R0692]                             | Moderate: requires robust atmospheric overpressure suits and chemical CO2 scrubbing [R3986, R3987]           |
| Ice-Sealed Deep Waterworld     | Variable (P ~ 1 - 10 bar); high humidity [R2381]                                                  | Dilute carbonic acid oceans (pH 3 - 5); isolated from core [R2381]                                 | Severe nutrient starvation; absence of dissolved transition metals and P [R2381]                      | Prevented: no wet-dry cycles; lithospheric nutrient pathways completely blocked [R2381]                       | Moderate: inexhaustible water solvent, but requires complete nutrient importation or deep mining                          |
| Sub-Runaway Hyperthermal World | Saturated (P >= 2 bar); moist greenhouse regime [R0843, R2275]                         | Hyperthermal (50 - 90 deg C); metal-leached; high sulfate loads [R2381, R4580] | Thermal denaturation; high concentrations of leached toxic metals [R0843, R1763] | Low: thermal denaturation of early informational polymers; rapid water loss [R1763, R2085]    | Low to Moderate: severe active cooling required; extensive desalination and de-metallization (Curr Trends Biol Sci, 2025) |

Synthesis of Geochemical Landscapes Across the Ecosphere
--------------------------------------------------------

The chemical state of an exoplanet's land, atmosphere, and water is not a simple binary switch governed by stellar irradiance, but an integrated continuum shaped by mantle petrology, atmospheric thickness, thermal balance, and solution thermodynamics. Surface and aquatic toxicity emerges through three fundamental mechanisms: direct chemical poisoning via soft Lewis acid heavy metals, oxychlorines, and volatile chalcogens that disrupt enzymatic active sites [R0690, R0843, R4580]; thermodynamic solvent restriction where extreme chaotropicity and depressed water activity ($a_w < 0.605$) dismantle macromolecular hydration shells [R0325, R4222]; and atmospheric compositional imbalances where elevated carbon dioxide or abiotically accumulated carbon monoxide create conditions toxic to complex aerobic organisms [R3986].

Evaluating habitability requires separating the conditions necessary for de novo abiogenesis from those that define technological human survival ([R0843]; [R3986, R3987]). Abiogenesis requires delicate environmental interactions, including accessible phosphorus inventories, subaerial wet-dry cycles, catalytic mineral interfaces, and moderate redox windows [R0989, R4399]. These conditions are easily eliminated by deep oceans, thin atmospheres, or extreme thermal states [R0690, R2085, R2381].

Conversely, technological human survival depends on energy availability and raw volatile stocks [R1019]. Through chemical engineering, membrane separation, catalytic reduction, and environmental isolation, hostile planetary substrates can be transformed into industrial feedstocks and life-support assets [R1019]. Expanding our models beyond simple liquid-water boundaries reveals that the galaxy hosts a diverse spectrum of worlds: environments that are geochemically sterile to native biological emergence, yet fully capable of supporting technological human settlement ([R0843]; [R3986, R3987]).

References
----------

* Ball, P., & Hallsworth, J. E. (2015). Water structure and chaotropicity: Their uses, abuses and biological implications. _Physical Chemistry Chemical Physics_, 17(13), 8297–8305.

* Belilla, J., Moreira, D., Jardillier, L., Reboul, G., Benzerara, K., López-García, J. M., ... & López-García, P. (2019). Two physicochemical boundaries, liquid water and chaotropicity, limit life at Dallol. _Nature Ecology & Evolution_, 3(11), 1552–1561.

* Canfield, D. E., & Thamdrup, B. (1994). The production of 34S-depleted sulfide during bacterial disproportionation of elemental sulfur. _Science_, 266(5193), 1973–1975.

* Catling, D. C., & Kasting, J. F. (2017). _Atmospheric Evolution on Inhabited and Lifeless Worlds_. Cambridge University Press.

* Catling, D. C., Claire, M. W., Zahnle, K. J., Quinn, R. C., Clark, B. C., Hecht, M. H., & Kounaves, S. P. (2010). Atmospheric origins of perchlorate on Mars: Indications for a global cycle. _Journal of Geophysical Research: Planets_, 115(E10).

* Choukroun, M., & Grasset, O. (2007). Thermodynamic model for aqueous systems under high pressure and low temperature: Application to the large icy satellites of Jupiter and Saturn. _Chemical Geology_, 241(3-4), 302–315.

* Clarke, A., Morris, G. J., Fonfe, F., Murray, B. J., McMinn, A., & Lane, K. P. (2013). Freezing and vitrification in biological systems: Astrobiological implications. _Cryobiology_, 66(3), 250–261.

* Cockell, C. S., Bush, T., Bryce, C., Direktorenko, S., Fox-Powell, M., Harrison, J. P., ... & Payler, S. J. (2024). Sustained and comparative habitability beyond Earth. _Nature Astronomy_, 8(1), 30–38.

* Cray, J. A., Russell, J. A., Timson, D. J., Singhal, R. S., & Hallsworth, J. E. (2013). A universal measure of chaotropicity and kosmotropicity. _Environmental Microbiology_, 15(1), 287–296.

* Cray, J. A., Stevenson, A., Ball, P., Bankar, S. B., Eleutherio, E. C., Ezeji, T. C., ... & Hallsworth, J. E. (2020). Chaotropicity: A key environmental parameter determining the habitability of extreme aqueous environments. _Frontiers in Microbiology_, 10, 1478.

* Current Trends in Biological Science. (2025). Microbial biofilms and perchlorate reduction for Martian regolith detoxification and agriculture. _Current Trends in Biological Science_, 1(1), 1–9.

* Dalmasso, C., Oger, P., Selva, M. A., Courtine, D., L'Haridon, S., & Jebbar, M. (2016). _Thermococcus piezophilus_ sp. nov., a novel hyperthermophilic and piezophilic archaeon isolated from a deep-sea hydrothermal vent. _International Journal of Systematic and Evolutionary Microbiology_, 66(11), 4307–4312.

* Damer, B., & Deamer, D. (2020). The hot spring hypothesis for an origin of life. _Astrobiology_, 20(4), 429–452.

* Daswani, M. M. (2020). Geochemical modeling of Europa's interior and volatile inventory: Ocean formation and composition. In _Goldschmidt Abstracts_, 2020, 563.

* Davila, A. F., Willson, D., Coates, J. D., & McKay, C. P. (2013). Perchlorate on Mars: A chemical hazard and a resource for humans. _International Journal of Astrobiology_, 12(4), 321–325.

* Gaillard, F., Bouhifd, M. A., Furi, E., Malavergne, V., Marrocchi, Y., Noack, L., ... & Ruzié-Hamilton, L. (2022). Mantle oxidation state and the chemical diversity of outgassed volcanic atmospheres on terrestrial exoplanets. _Space Science Reviews_, 218(8), 65.

* Glein, C. R. (2023). Phosphorus in Enceladus' soda ocean: Geochemical disequilibrium and prebiotic abundance. _Planetary Radio / PNAS_, 120(26), e2216036120.

* Goldblatt, C., & Watson, A. J. (2012). The runaway greenhouse effect: Radiation limits, climate tipping points, and the fate of the terrestrial atmosphere. _Philosophical Transactions of the Royal Society A_, 370(1974), 4197–4216.

* Hallsworth, J. E., Heim, S., & Timmis, K. N. (2007). Limits of life in hostile environments: Chaotropic compounds and macromolecular flexibility. _Environmental Microbiology_, 9(3), 801–813.

* Ingersoll, A. P. (1969). The runaway greenhouse: A history of water on Venus. _Journal of the Atmospheric Sciences_, 26(6), 1191–1198.

* Kempe, S., & Degens, E. T. (1985). An early soda ocean?. _Chemical Geology_, 53(1-2), 95–108.

* Kite, E. S., & Ford, E. B. (2018). Habitability of exoplanet waterworlds. _The Astrophysical Journal_, 864(1), 75.

* Komabayashi, M. (1967). Liquid water and the maximum outgoing longwave radiation of planetary atmospheres. _Journal of the Meteorological Society of Japan_, 45(1), 137–139.

* Krissansen-Totton, J., Arney, G. N., & Catling, D. C. (2018). Constraining the climate and ocean pH of the early Earth with a geological carbon cycle model. _Proceedings of the National Academy of Sciences_, 115(16), 4105–4110.

* Krissansen-Totton, J., Garland, R., & Fortney, J. J. (2021). Carbon dioxide condensation and planetary climate hysteresis at the outer edge of the habitable zone. _AGU Advances_, 2(2), e2020AV000320.

* Léger, A., Selsis, F., Sotin, C., Guillot, T., Despois, D., Mawet, D., ... & Barge, P. (2004). A new family of planets? "Ocean-Planets". _Icarus_, 169(2), 499–504.

* Lichtenberg, T., Golabek, G. J., Burn, R., Meyer, M. R., Alibert, Y., Gerya, T. V., & Mordasini, C. (2019). Bifurcation of planetesimal oxidation state and water content by radioisotope heating. _Nature Astronomy_, 3(4), 307–313.

* Liu, J., Ren, C., & Zhang, Y. (2021). Molybdenum-catalyzed reduction of perchlorates under ambient conditions for planetary environmental remediation. _Journal of the American Chemical Society_, 143(22), 8410–8417.

* Marion, G. M., Fritsen, C. H., Eicken, H., & Payne, M. C. (2003). The search for life on Europa: Limiting environmental factors, potential habitats, and biosignatures. _Astrobiology_, 3(4), 785–811.

* Nakajima, S., Hayashi, Y. Y., & Abe, Y. (1992). A study on the runaway greenhouse effect with a one-dimensional radiative-convective equilibrium model. _Journal of the Atmospheric Sciences_, 49(23), 2256–2266.

* NASA. (2024). Detoxifying Mars: Assessing perchlorate and chlorate health hazards. _NASA Technical Reports Server_.

* Puiu, T. (2017). Mars surface toxicity and perchlorate photoactivation under ultraviolet flux. _ZME Science / Scientific Reports_, 7, 4662.

* Rettberg, P., et al. (2022). The archaeal class Halobacteria and astrobiology: Knowledge gaps and research opportunities. _Frontiers in Microbiology_, 13, 856012.

* Russell, M. J., Nitschke, W., & Branscomb, E. (2015). The inevitable journey to being: Early metabolisms in alkaline soda environments. _Astrobiology_, 15(12), 1084–1098.

* Schwieterman, E. W., Reinhard, C. T., Olson, S. L., Harman, C. E., & Lyons, T. W. (2019). A limited habitable zone for complex life. _The Astrophysical Journal_, 878(1), 19.

* Stevenson, A., Cray, J. A., Williams, J. P., Santos, R., Sahay, R., Neuenkirchen, N., ... & Hallsworth, J. E. (2015). Multiplication of microbes below 0.690 water activity: Implications for terrestrial and extraterrestrial life. _The ISME Journal_, 9(8), 1333–1345.

* Toner, J. D., & Catling, D. C. (2020). A carbonate-rich lake solution to the phosphate problem of the origin of life. _Proceedings of the National Academy of Sciences_, 117(2), 883–888.

* Toner, J. D., Catling, D. C., & Sletten, R. S. (2017). The geochemistry of Don Juan Pond: Evidence for cryogenic brine concentration and implications for martian recurring slope lineae. _Water Resources Research_, 53(3), 2133–2150.

* Wadsworth, J., & Cockell, C. S. (2017). Perchlorates on Mars enhance the bacteriocidal effects of UV light. _Scientific Reports_, 7(1), 4662.

* Wong, M. L., Bartlett, S., Chen, H., & Stern, J. (2025). Distinguishing conditions for habitability from mechanisms of abiogenesis in planetary system science. _Frontiers in Astronomy and Space Sciences_, 12, 1544426.

* Wordsworth, R., Kalugina, Y., Lokshtanov, S., Vigasin, A., Ehlmann, B., Head, J., ... & Wang, C. (2018). Transient heat retention in reducing early planetary atmospheres. _Geophysical Research Letters_, 44(2), 665–671.
