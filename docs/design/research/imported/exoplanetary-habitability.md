> **Imported research, kept as received (banner added 2026-10-15; the text below is unchanged).** It is Boss's document, not a statement of what planetGen does. The design notes that use it: [habitability-index](../../habitability/habitability-index.md).

# Exoplanetary Habitability and Human Sustainability: An Ecological and Physiological Classification Architecture

The assessment of extrasolar planetary habitability requires a rigorous analytical division between microbial viability and human sustainability. While astrobiological metrics such as the Earth Similarity Index (ESI) and the Planetary Habitability Index (PHI) evaluate whether a planetary body can sustain liquid surface water or host carbon-based extremophiles, human physiological requirements impose far narrower physical, atmospheric, and biochemical boundaries [R0841, R3979].

Human metabolic homeostasis requires a precise partial pressure envelope of oxygen and inert buffer gases, tight thermodynamic bounds, continuous hydration, protection from acute and chronic ionizing radiation, and biochemically compatible caloric input.

To systematically evaluate exoplanetary environments, planetary bodies can be classified along a four-tier scale: Red (completely hostile/non-sustainable), Yellow (marginal/high-dependency survival), Green (semi-autonomous/partially sustainable), and Blue (optimal/superhabitable). Evaluating these environments along axes of water availability, natural shelter geomorphology, and ecological complexity reveals the theoretical survival limits of human populations across individual, expeditionary group, and permanent colonial scales.

---

## 1. Physiological Baselines and Astrobiological Boundaries

Human physiological vulnerability defines the absolute boundary conditions for habitability classification. Atmospheric composition and barometric pressure govern alveolar gas exchange through Dalton’s and Henry’s laws. At total ambient pressures below $6.3\text{ kPa}$ (the Armstrong limit), water vapor reaches its boiling point at normal human core temperature ($37^\circ\text{C}$), precipitating immediate bodily fluid ebullism, subcutaneous tissue expansion, pulmonary atelectasis, alveolar rupture, and acute circulatory collapse.

Even when ambient pressure exceeds this mechanical threshold, arterial oxygen saturation depends entirely on the alveolar partial pressure of oxygen ($pO_2$). An ambient $pO_2$ below $6.0\text{ kPa}$ prevents hemoglobin saturation and reverses diffusion across the alveolar-capillary membrane, venting oxygen from venous blood directly into the atmosphere and inducing hypoxic blackout within seconds. Conversely, sustained exposure to ambient $pO_2$ levels exceeding $45\text{--}50\text{ kPa}$ induces central nervous system toxicity, membrane lipid peroxidation, and fatal pulmonary edema.

Inert buffer gases—predominantly diatomic nitrogen ($N_2$)—must supply between $10\text{ kPa}$ and $30\text{ kPa}$ of partial pressure to maintain alveolar structural inflation and suppress hyperoxic flammability hazards. Simultaneously, the ambient carbon dioxide partial pressure ($pCO_2$) must remain below $0.5\text{--}1.0\text{ kPa}$; levels above $1.0\text{ kPa}$ cause progressive respiratory acidosis, cognitive decline, hypercapnic narcosis, and eventual unconsciousness.

Thermal tolerance is dictated by the thermodynamic limits of metabolic heat rejection. The human body generates continuous metabolic heat that must be dissipated into the surrounding atmosphere via convection, radiation, and evaporative cooling through perspiration. This heat rejection mechanism ceases when environmental conditions exceed a critical wet-bulb temperature ($T_w$) of $31\text{--}35^\circ\text{C}$. Beyond this threshold, core body temperature increases uncontrollably, precipitating fatal heat stroke and cardiovascular collapse within hours, regardless of hydration or convective airflow.

At the cold extreme, dry-bulb temperatures below $0^\circ\text{C}$ trigger peripheral vasoconstriction, rapid caloric depletion, tissue ischemia, and lethal systemic hypothermia unless continuous external thermal management and insulation are maintained.

Ionizing radiation presents another fundamental barrier to long-term survival. Baseline terrestrial exposure averages under $1\text{ mSv/year}$, with standard industrial occupational limits set between $20\text{ mSv}$ and $50\text{ mSv/year}$. Acute surface exposures exceeding $0.75\text{--}1.0\text{ Sv}$ trigger acute radiation syndrome (ARS), resulting in bone marrow destruction, mucosal sloughing, and systemic immune collapse [R0952]. Doses exceeding $4.5\text{--}5.0\text{ Sv}$ represent an unassisted lethal dose for fifty percent of an exposed cohort ($LD_{50/30}$).

Multi-generational colonial settlement requires annual background exposure to remain below $50\text{--}100\text{ mSv/year}$ to prevent reproductive sterility, elevated teratogenesis, and embryonic lethality, necessitating substantial mass shielding against both galactic cosmic rays (GCR) and stellar energetic particles (SEP).

---

## 2. Four-Tier Classification Architecture

The transition between complete environmental non-viability and ecological optimization can be categorized into four distinct planetary tiers. The physical parameters defining these regimes reflect the complex interplay among stellar spectral irradiance, planetary mass, atmospheric loss mechanisms, and internal geodynamic longevity.

### Table 1: Physical and Astrophysical Parameter Matrix Across Tiers

| Parameter | Red Tier (Hostile) | Yellow Tier (Marginal) | Green Tier (Viable) | Blue Tier (Optimal) |
| :--- | :--- | :--- | :--- | :--- |
| **Earth Similarity Index (ESI)**<br>[R3979] | $< 0.40$ | $0.40\text{--}0.74$ | $0.75\text{--}0.89$ | $\ge 0.90$ |
| **Total Surface Pressure**<br> | $< 6.3\text{ kPa}$ or $> 500\text{ kPa}$ | $20\text{--}70\text{ kPa}$ or $150\text{--}300\text{ kPa}$ | $70\text{--}150\text{ kPa}$ | $90\text{--}120\text{ kPa}$ |
| **$pO_2$ Availability**<br>[R0841] | $0.0\text{ kPa}$ (Trace) | $< 6.0\text{ kPa}$ | $12.0\text{--}18.0\text{ kPa}$ | $19.0\text{--}24.0\text{ kPa}$ |
| **$pCO_2$ Concentration**<br>[R0841] | Unbounded / $> 10\text{ kPa}$ | $0.5\text{--}5.0\text{ kPa}$ | $0.04\text{--}0.5\text{ kPa}$ | $< 0.05\text{ kPa}$ |
| **Mean Surface Temperature**<br>[R0841] | $<-80^\circ\text{C}$ or $>+70^\circ\text{C}$ | $-50^\circ\text{C}\text{ to }+45^\circ\text{C}$ | $-10^\circ\text{C}\text{ to }+35^\circ\text{C}$ | $+10^\circ\text{C}\text{ to }+25^\circ\text{C}$ |
| **Surface Radiation Flux**<br>[R0952, R3933] | $> 1.0\text{ Sv/yr}$ (Unshielded) | $0.1\text{--}1.0\text{ Sv/yr}$ | $5\text{--}50\text{ mSv/yr}$ | $< 5\text{ mSv/yr}$ |
| **Host Stellar Spectral Class**<br>[R2769, R4064, R4308] | Active Flare M-dwarf, O, B, A | Active M-dwarf, Early F | Late K-dwarf, Inactive M | Stable K-dwarf, G-dwarf |
| **Hydrosphere State**<br>[R0841, R2281] | Cryogenic ice, perchlorates, vapor | Permafrost, sub-ice ocean, brine | Open rivers, lakes, active rain | Pervasive oceans, shallow seas |

### Red Tier: Non-Sustainable and Completely Hostile Environments
Red-tier worlds represent environments completely devoid of biological survivability. These planets are characterized by atmospheric pressures falling below the Armstrong limit or exceeding hundreds of bars of superheated, corrosive greenhouse gases. Solar system analogues include Mercury, Venus, Mars, and Jupiter’s moon Io.

Surface radiation fluxes on Red-tier bodies regularly exceed occupational limits by three to four orders of magnitude due to the lack of an intrinsic magnetic dipole or atmospheric mass column, or because the body orbits within the severe charged-particle trapping belts of a gas giant primary. Native surface biological activity is entirely absent, with the theoretical exception of endolithic, chemoautotrophic micro-organisms inhabiting deep hydrothermal fracture zones kilometers beneath the crust.

### Yellow Tier: Marginal and High-Dependency Regimes
Yellow-tier planets maintain sufficient mass and barometric density to prevent bodily fluid ebullism and support stable liquid phases of volatiles, but they lack a breathable atmospheric gas mixture. These environments parallel Archean Earth, contemporary Titan, or a geologically warm, wet Mars.

The atmospheric column is dominated by molecular nitrogen, carbon dioxide, or methane, and frequently contains toxic trace gases such as hydrogen sulfide, carbon monoxide, or sulfur dioxide. While the dense atmosphere shields the planetary surface from solar ultraviolet-C radiation and micrometeoroid ablation, a weak magnetic dipole may permit substantial secondary neutron cascades at the surface. Native life on Yellow-tier worlds is restricted to single-celled anaerobic archaea, methanogens, and anoxygenic phototrophs that generate biosignatures without producing significant free atmospheric oxygen.

### Green Tier: Viable and Semi-Autonomous Ecosystems
Green-tier planets exhibit an active hydrologic cycle, an endogenous core geodynamo, and a biogenically oxygenated atmosphere. These worlds reflect Earth during the late Proterozoic to early Paleozoic transitions.

Photosynthetic primary producers actively maintain a global carbon-oxygen-nitrogen cycle, producing surface oxygen levels sufficient to sustain human metabolism ($pO_2$ between $12\text{ kPa}$ and $18\text{ kPa}$). However, these atmospheres may also exhibit localized hypoxia or mild baseline hypercapnia.

Indigenous autotrophic life includes complex multicellular algae, mats of bryophyte analogues, and primitive vascular flora. Heterotrophic fauna is predominantly invertebrate-grade, marine, or primitive detritivorous organisms. Human beings can traverse the surface of Green-tier worlds without heavy pressure suits, requiring only lightweight respiratory assist systems, particulate filters, or localized oxygen enrichment masks.

### Blue Tier: Optimal and Superhabitable Habitats
Blue-tier worlds match or exceed Earth’s baseline ecological carrying capacity, often meeting the criteria for superhabitability established by [R1897]. These planets are typically super-Earths ($1.2\text{--}2.0\text{ }M_\oplus$, $1.1\text{--}1.3\text{ }R_\oplus$) orbiting long-lived, quiescent K-type dwarf stars.

Their greater planetary mass sustains internal radiogenic heat engines longer than Earth, driving stable plate tectonics and preserving magnetic field protection over $5\text{--}8\text{ billion years}$. Their surface topography is dominated by distributed continental archipelagos and warm, shallow epicontinental seas rather than deep ocean basins, promoting exceptional biological diversity and stable climate feedback loops.

Surface pressures, atmospheric gas mixes, and mean global temperatures ($15\text{--}20^\circ\text{C}$) align precisely with human homeostatic needs. Stratified ecosystems of complex autotrophs and high-order animal analogues populate every biome, allowing unassisted human survival, foraging, and multi-generational colonial growth.

---

## 3. Environmental Resource Dimensions

The sustainability of an exoplanet is governed by the operational availability of three interdependent physical resources: liquid water, natural shielding geomorphology, and biochemically compatible biomass. The interplay among these variables determines whether an environment requires energy-intensive mechanical synthesis or provides open, naturally buffered ecosystem services.

### Hydrological Distribution and Extraction Energetics
Water availability transitions from severe energetic deficit to an ambient, self-replenishing state across the four planetary tiers. On Red-tier worlds, liquid water cannot exist stably on the surface due to thermodynamic constraints; it is confined to cryogenic permafrost deposits, deep sub-surface aquifers, hydrated phyllosilicate clays, or high-salinity eutectic brines rich in toxic magnesium and sodium perchlorates. Extracting potable water requires mechanical mining, high-temperature thermal extraction, catalytic perchlorate destruction, and multistage membrane filtration.

On Yellow-tier planets, water is abundant in the form of extensive glacial caps, surface ice sheets, and subglacial oceans. However, low ambient temperatures and the absence of an open hydrologic cycle require dedicated thermal melting systems and desalination processing to remove heavy metal complexes and volatile sulfides.

Green-tier planets maintain active hydrological networks featuring liquid rivers, freshwater lakes, and regular precipitation. Water purification on these worlds requires only biological filtration to denature alien organic macromolecules or cellular fragments, eliminating the need for industrial distillation infrastructure.

Blue-tier worlds possess an ideal hydrologic system, offering clean, naturally filtered freshwater that requires zero technological intervention for human consumption.

### Geomorphological Shelter and Radiative Attenuation
Planetary shelter availability is determined by the presence of natural landforms that mitigate extreme thermal swings and ambient radiation. In Red and early Yellow environments, the absence of thick atmospheres and strong magnetospheres leaves the surface vulnerable to primary galactic cosmic rays and dangerous solar particle events.

Survival on these worlds depends on using geological formations as passive radiation shields. Subsurface basaltic lava tubes—formed during ancient effusive volcanic eruptions—provide cavernous void spaces beneath roof layers several meters to dozens of meters thick [R1840, R3933, R4332]. These subterranean caverns provide natural protection, reducing radiation levels to safe, earth-normal limits while providing stable ambient temperatures that isolate habitats from surface thermal extremes. In areas lacking volcanic tubes, surface habitats require mechanical burial under $300\text{--}500\text{ g/cm}^2$ of excavated regolith.

On Green-tier worlds, an atmospheric mass column exceeding $1000\text{ g/cm}^2$ shields the surface from cosmic particles, while an active stratospheric ozone layer filters harmful ultraviolet radiation. Human shelters can be constructed directly on the surface using locally sourced stone, earth, and structural timber analogues.

Blue-tier planets provide complete natural environmental protection, permitting unshielded above-ground habitation across most latitudes.

### Autotrophic Dynamics, Xenobiology, and Caloric Viability
The presence of native flora and fauna introduces the critical challenge of biological compatibility. Terrestrial biology is defined by strict homochirality: functional enzymatic proteins rely on left-handed (L-enantiomer) amino acids, while metabolic sugars and nucleic acids use right-handed (D-enantiomer) conformations [R0220, R2207, R3024].

On independently evolved Yellow and Green exoplanets, native life is equally likely to exhibit reversed chirality (D-amino acids and L-sugars). Ingesting mirrored organic matter provides zero nutritional value; human digestive enzymes cannot hydrolyze D-peptide bonds, and metabolic processing of foreign enantiomers can produce toxic byproducts. Even when chiral symmetry matches Earth biology, foreign flora may incorporate non-canonical amino acids outside the 22 utilized by terrestrial life, or synthesize defensive xenobiotic alkaloids and accumulate toxic metalloproteins.

Photosynthetic processes on planets orbiting M-dwarf stars adapt to spectral energy distributions dominated by long-wavelength photons [R4064]. Autotrophs on these worlds use specialized pigments such as chlorophyll d and chlorophyll f, shifting their absorption spectra into the far-red and near-infrared bands ($700\text{--}750\text{ nm}$) ([R0760, R2344], 2007b). This spectral adaptation alters the characteristic "vegetation red edge" seen in reflectance spectroscopy. While these autotrophs can produce substantial atmospheric oxygen, their thermodynamic efficiency and protein profiles make them poor nutritional substitutes for terrestrial crops without secondary refining.

Only on Blue-tier worlds—where evolutionary history has converged on compatible organic structures or where planetary engineering has established terrestrial biomes—can human beings forage native flora and fauna safely without risking severe immunological rejection or fatal metabolic toxicity.

### Table 2: Ecological Resource Availability and Xenobiotic Profiles Across Tiers

| Resource Dimension | Red Tier (Hostile) | Yellow Tier (Marginal) | Green Tier (Viable) | Blue Tier (Optimal) |
| :--- | :--- | :--- | :--- | :--- |
| **Water Accessibility**<br>[R0841, R0952] | Deep hydration bounds; perchlorate brines; sublimation | Glacial ice sheets; subglacial liquid reservoirs | Flowing rivers; freshwater lakes; rain cycles | Global freshwater networks; optimal precipitation |
| **Shelter Geomorphology**<br>[R1840, R3933, R4332] | Deep basalt lava tubes; sintered regolith mounds | Volcanic tubes; deep tectonic rift valleys | Natural caves; rock overhangs; surface construction | Unrestricted surface terrain; structural timber |
| **Flora & Autotrophs**<br>([R0760, R2344], 2007b) | Completely abiotic surface | Chemosynthetic mats; far-red/NIR cyanobacteria | Multicellular bryophytes; red-edge shifted canopies | Fully developed forest canopies; high biomass yield |
| **Biochemical Compatibility**<br>[R0220, R2207, R3024] | Sterile mineral substrate | Incompatible enantiomers; abiotic xenotoxins | Enantiomeric divergence; novel secondary metabolites | Homochiral convergence (L-amino acids, D-sugars) |
| **Fauna & Heterotrophs**<br>[R0841, R4064, R4621] | None | Sub-surface anaerobic extremophiles | Invertebrate-grade detritivores; marine macrofauna | Complex ecological trophic webs; vertebrate analogues |

---

## 4. Human Survivability and Population Scaling Dynamics

The duration of human survival across different planetary environments depends on population size and life-support technology. Analyzing survival limits across the individual, the expeditionary team, and the permanent settlement reveals the primary failure modes of each operational scale.

### Individual Exposure: Unassisted Decompression Versus Technical Augmentation
An unprotected human exposed to the surface of a Red-tier world experiences immediate physiological failure. In a near-vacuum environment below the Armstrong limit, the Time of Useful Consciousness (TUC) is limited to $9\text{--}15\text{ seconds}$. Because the partial pressure of ambient oxygen is zero, the concentration gradient across the pulmonary alveolar-capillary membrane reverses, stripping oxygen from venous blood and venting it directly into the environment. Deoxygenated blood circulates to the brain within seconds, inducing loss of consciousness.

Simultaneously, ebullism causes rapid vapor bubble formation within the bloodstream and soft tissues, leading to circulatory collapse and irreversible cerebral anoxia within $90\text{ to }180\text{ seconds}$ [R3997]. When equipped with a high-pressure extravehicular activity (EVA) suit, survival is extended to $6\text{--}8\text{ hours}$, governed by the chemical saturation of carbon dioxide scrubbers, battery reserves, and water sublimation cooling loops [R4685].

On Yellow-tier planets where total atmospheric pressure exceeds the Armstrong limit but lacks oxygen (such as a cold $CO_2\text{--}N_2$ atmosphere), ebullism does not occur. However, inhaling gas with a high $pCO_2$ causes rapid hypercapnic narcosis and asphyxia, producing unconsciousness within $30\text{--}60\text{ seconds}$ and death from respiratory acidosis within minutes. Equipped with a basic closed-circuit rebreather or an oxygen concentrator and thermal protective garments, an individual can operate without heavy, pressurized suits. In this mode, individual survival extends to days or weeks, limited only by the availability of rations and clean water.

On Green-tier worlds, an unassisted individual can breathe ambient air directly, but long-term survival is bounded by nutritional and chemical factors. The individual may suffer from chronic hypoxia or mild acidosis if the gas mixture deviates from Earth-normal levels. The primary limit to unassisted survival is the 2-to-4-week starvation window; native plants and animals cannot provide usable calories if they possess incompatible chirality or contain novel biochemical toxins.

On Blue-tier worlds, unassisted survival extends to the full human lifespan ($70\text{--}90\text{ years}$), with mortality driven entirely by standard terrestrial health factors and environmental hazards.

### Expeditionary Group Logistics: Months to Multi-Year Sourcing
A small expeditionary group of four to twelve individuals requires complex environmental management and consumable replenishment. On Red-tier worlds, survival depends on deploying hermetically sealed, pressurized habitat modules equipped with physicochemical Environmental Control and Life Support Systems (ECLSS). The habitat must incorporate catalytic Sabatier reactors to recycle hydrogen and carbon dioxide into potable water, alongside solid polymer water electrolysis to recover oxygen.

Because each crew member requires approximately $2.5\text{--}3.0\text{ kg}$ of water daily for hydration and metabolic support, missions without closed-loop recycling face a hard logistics limit of $30\text{ to }90\text{ days}$. Extending mission duration beyond this threshold requires active In-Situ Resource Utilization (ISRU) to mine surface permafrost or hydrated minerals for habitat water and oxygen replenishment.

On Yellow-tier worlds, expeditionary teams face lower structural and logistical mass requirements. Because the ambient atmospheric pressure prevents rapid explosive decompression, habitats can be constructed using deployable, lightweight membrane structures rather than high-mass pressure shells. Oxygen generation does not require recycling internal metabolic waste; compressors can harvest ambient carbon dioxide directly from the atmosphere and crack it into breathable oxygen using solid-oxide electrolysis systems. The primary barrier to multi-month missions shifts from atmospheric support to nutritional autonomy, as foreign soil chemistry and a lack of bio-available nitrogen necessitate importing fertilizers or maintaining closed hydroponic systems [R4685].

On Green-tier planets, an expeditionary team can achieve continuous survival across months and years. Natural caves or basic surface structures built from native stone provide sufficient shelter without radiation-shielding berms. Water is directly sourced from surface rivers and lakes, requiring only microbial sterilization. Small-scale agriculture can be initiated by treating local soils to remove novel plant toxins and inoculating the substrate with terrestrial mycorrhizal fungi and nitrogen-fixing bacteria, providing sustainable caloric replenishment.

On Blue-tier worlds, an expeditionary group operates with complete environmental autonomy from arrival, utilizing standard field equipment and foraging native biomes directly.

### Multi-Generational Colonial Settlement: Closed-Loop BLSS and Biome Integration
Long-term planetary colonization involving hundreds to thousands of individuals shifts the survival model from open supply chains to fully closed Bioregenerative Life Support Systems (BLSS) [R1166]. On Red and Yellow worlds, life support cannot rely indefinitely on mechanical consumables and replacement parts; it requires biological cycles to process waste and produce food.

Higher plants consume carbon dioxide and human metabolic wastes, producing oxygen, purifying water through transpiration, and generating edible biomass [R4685]. Life-support models show that nutritional and atmospheric autonomy requires approximately $40\text{--}50\text{ m}^2$ of illuminated crop growing area per person, alongside $15\text{--}20\text{ m}^2$ for respiratory gas balancing [R1166]. A self-sustaining colony of 1,000 people on a Yellow-tier world requires at least $50,000\text{ m}^2$ (5 hectares) of high-yield hydroponic or aeroponic growing surfaces protected within subterranean lava tubes or heavily shielded surface domes. Aquaculture systems (such as tilapia, carp, or microalgae) are integrated into these loops to achieve higher feed conversion efficiency than traditional terrestrial livestock. The colony's long-term survival is limited by mechanical failure rates, degradation of synthetic polymers, trace heavy metal accumulation in recycled water loops, and the long-term genetic viability of closed agricultural crops.

On Green-tier planets, a permanent colony avoids the severe mass penalties of hermetic pressure hulls and subterranean living. Colonists can build standard surface housing using local timber, stone, and bio-composites. Agriculture can be conducted outdoors, provided the native soil microbiome is pre-treated or outcompeted with terrestrial microbes. The primary threat to generational survival is biological: novel planetary allergens, cross-reactive toxins, and xenobiotic parasites that bypass terrestrial immune defenses [R2006].

On Blue-tier worlds, colonies can expand openly across the landscape. Agriculture, ecological management, and urban development mirror terrestrial patterns without the structural mass, power consumption, or life-support overhead required on harsher planets.

### Table 3: Population Survivability Timelines and Primary Failure Modes

| Population Scale | Red Tier (Hostile) | Yellow Tier (Marginal) | Green Tier (Viable) | Blue Tier (Optimal) |
| :--- | :--- | :--- | :--- | :--- |
| **Individual (Unassisted)**<br>[R3997] | 9–15 seconds (TUC); fatal ebullism/asphyxia in 3 min | 1–10 minutes; hypercapnic narcosis; acute anoxia; thermal shock | 2–4 weeks; chronic hypoxia/hypercapnia; chiral food starvation | Full Natural Lifespan ($70\text{--}90\text{ years}$); ambient equilibrium |
| **Individual (Assisted)**<br>[R4685] | 6–8 hours; EVA life-support limits (scrubber, $O_2$, power) | Days to Weeks; limited by $O_2$ concentrates, battery storage, and rations | Months to Years; requires simple respiratory filters and micro-rations | Indefinite; technology serves as comfort enhancement rather than life support |
| **Expeditionary Group**<br>[R4685] | 1–3 months; requires rigid ECLSS modules; rapid consumable depletion | 6–12 months; deployable habitats; ISRU gas cracking; water ice melting | Indefinite (Years); field cabins; surface water sourcing; basic farming | Indefinite; immediate autonomous survival and field foraging |
| **Permanent Colony**<br>[R1166, R4685] | Fragile / Finite; dependent on $100\%$ closed-loop BLSS and structural hulls | Sustainable with Heavy ISRU; $50\text{ m}^2$ crop area/person; lava tube city | Highly Autonomous; open-air agriculture; simple soil bio-remediation | Fully Autonomous; open biome expansion; multi-generational growth |
| **Primary Mode of Mortality**<br>[R0952, R2006, R3997] | Explosive decompression; ebullism; acute lethal radiation | Suffocation ($CO_2/N_2$); hypothermia; mechanical ECLSS failure | Xenobiotic toxicity; allergic shock; novel pathobiology | Terrestrial senescence; ecological predation; standard pathology |

---

## 5. Comparative Ecological Synthesis

The boundary between planetary tiers reflects a transition from complete technological life support to self-sustaining biological systems.

On Red-tier planets, survival is entirely technological: every liter of water, gram of oxygen, and calorie of nutrition must be mechanically processed and recycled using dedicated energy sources. The environment offers no supportive ecosystem services, meaning any system failure directly threatens human life.

Yellow-tier worlds introduce a hybrid operational profile. Planetary mass and atmospheric density provide passive radiation shielding and thermal buffering, substantially reducing the energy footprint required for environmental maintenance. This stability allows colonies to implement large-scale bioregenerative life support systems. However, true ecological independence remains impossible on Yellow worlds. The absence of an indigenous, oxygenated nitrogen-carbon cycle forces agriculture to remain isolated in sealed, artificial environments.

The transition to Green and Blue tiers marks an ecological threshold: the exoplanet's native biosphere actively maintains atmospheric, hydrological, and thermal equilibrium. Yet, this biological presence introduces a unique ecological challenge: the xenobiological barrier. While sterile worlds (Red and early Yellow) present only physical engineering and chemical hazards, biotically active worlds (Green and non-terrestrial Blue) introduce complex organic risks, including biochemical toxicity, chiral incompatibilities, and novel pathobiological reactions.

Long-term survival on these planets depends on overcoming the immunological and metabolic mismatches between distinct evolutionary lineages. When this biological hurdle is overcome—or on Blue-tier worlds where the underlying biochemistry is natively compatible—human survival transitions from life-support maintenance to open ecological integration, enabling permanent, multi-generational habitation.

---

### References

[R4408]

[R0760]

[R1159]

[R4092]

[R1840]

[R1897]

[R2006]

[R2207]

[R2281]

[R2344]

[R2345]

[R2769]

[R3024]

[R3933]

[R3979]

[R3997]

[R4064]

[R4308]

[R4332]

[R4621]

[R5003]

[R4685]