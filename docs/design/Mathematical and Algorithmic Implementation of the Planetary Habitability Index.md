# Mathematical and Algorithmic Implementation of the Planetary Habitability Index (PHI-4)

## 1. Architectural Partitioning and Pipeline Architecture

To calculate or generate scientifically consistent values across the four core domains of the Planetary Habitability Index (PHI-4)—Pressure, Temperature, Chemistry, and Radiation—a mathematical and physical framework must bridge fundamental planetary astrophysics with geophysical and biophysical constraints. A computational pipeline must distinguish between stochastic boundary conditions (randomly generated primary cosmological and bulk properties) and deterministic planetary physics (outright calculated properties). 

### 1.1 Category A: Primary Stochastic Parameters (Randomly Sampled Boundary Conditions)

These parameters cannot be derived from first principles in a 1D pipeline but have well-studied prior distributions:

* **Host Star Characteristics:** Stellar mass ($M_* \in [0.08, 1.4]\,M_\odot$) is sampled via a Kroupa initial mass function (IMF), along with effective temperature ($T_*$), metallicity ($[\text{Fe}/\text{H}]$), and system age ($t_{\text{sys}} \in [0.5, 10]\text{ Gyr}$).
* **Orbital Architecture:** The semi-major axis ($a$) is drawn from a log-uniform distribution restricted to the broadened habitable zone, and orbital eccentricity ($e \in [0, 0.4]$) is sampled.
* **Planetary Bulk Endowments:** Planetary mass ($M_p \in [0.1, 5.0]\,M_\oplus$), core mass fraction ($\text{CMF} \in [0.2, 0.7]$), and total volatile water mass fraction ($\text{WMF} \in [10^{-5}, 0.1]$) are assigned using truncated distributions.
* **Planetary Rotation Period ($\Omega$) and Tidal Locking:** The rotation rate $\Omega = 2\pi / P_{\text{rot}}$ dictates Coriolis forces, atmospheric circulation regimes, and the geodynamo dipole moment $\mathcal{M}_p$. If system age $t_{\text{sys}} \ge t_{\text{sync}}$, the planet is tidally locked; otherwise, $P_{\text{rot}}$ is sampled from a Log-Normal distribution.
* **Stellar Activity and Spectra:** Flare kinetic energies follow a cumulative power-law frequency distribution, while Far-UV to Near-UV stellar flux ratios ($\mathcal{R}_{\text{UV}}$) are sampled from conditional Normal distributions fitted to empirical observations.
* **Surface Hypsometry and Cloud Coverage:** Subaerial land exposure fraction ($f_{\text{land}}$) requires modeling crustal elevation probability density $P(z)$ as a bimodal Gaussian mixture. Global cloud fraction $f_{\text{cloud}}$ is drawn from a Beta distribution to modify Bond albedo dynamically.
* **Aqueous Solution Ionic Molalities ($m_i$):** Ionic breakdowns are drawn using a Dirichlet distribution anchored by the mantle oxidation state ($\Delta\text{FMQ}$) to establish brine compositions.

### 1.2 Category B: External Inputs (Contextual Parameters)

Certain parameters require explicit upstream configuration to prevent arbitrary or mathematically unphysical results:

* **Target Organism Biological Profile ($\vec{\theta}_{\text{bio}}$):** Boundary values ($T_{\text{wb}} \le 35^\circ\text{C}$, $p_A\text{O}_2 \ge 10\text{ kPa}$, $\dot{D}_{\text{surf}} \le 10\text{ Gy yr}^{-1}$) require explicit organism definition (e.g., complex aerobic endotherms vs. extremophiles) to define the transfer functions.
* **Domain Importance Weights ($\vec{w}$):** The synthesis vector $\vec{w} = (w_{\text{Pressure}}, w_{\text{Temperature}}, w_{\text{Chemistry}}, w_{\text{Radiation}})$ is required to configure assessment priorities across domains.
* **Planetary Tectonic Regime ($\Theta_{\text{tect}}$):** A categorical state flag determines mantle outgassing fluxes and subduction recycling rates over gigayear timescales.
* **Primordial Volatile Inventory ($\vec{N}_{\text{vol}}$):** Initial molar endowments of $\text{N}_2$, $\text{CO}_2$, and $\text{H}_2\text{O}$ are required from upstream accretion models to prevent unrealistic runaway or airless conditions.
* **3D Atmospheric Advection Efficiency ($\epsilon_{\text{adv}}$):** Essential for heat redistribution on slow-rotating or tidally locked worlds, mapped via interpolated lookup matrices derived from Global Climate Models (GCMs).

## 2. Mathematical Framework by PHI-4 Domain (Deterministically Calculated Outputs)

### 2.1 Domain 1: Pressure & Planetary Structure

* **Hydrostatic Equilibrium and Planetary Mass-Radius Equations of State (EOS):** Requires non-linear differential equations and Birch-Murnaghan or Vinet EOS modeling to determine planetary internal density profiles, radius, and surface gravity ($g = GM/R^2$), which sets the surface weight of the volatile column. Interior structure models solve radial mass conservation:
  $$\frac{dM(r)}{dr} = 4\pi r^2 \rho(r), \quad \frac{dP(r)}{dr} = -\frac{G M(r) \rho(r)}{r^2}$$
* **Atmospheric Column Mass Integration and Barometric Formulations:** Uses the barometric formula and vertical column mass integration ($X = \int_{0}^{\infty} \rho(z)\,dz = P_0/g$) to establish surface barometric pressure ($P_0$), scale height ($H = k_B T / \mu g$), and total column depth for radiation/respiratory modeling.
* **Hydrodynamic and Thermal Atmospheric Escape:** Employs Maxwell-Boltzmann velocity distributions, thermal Jeans escape parameter ($\lambda_{\text{esc}, i}$), and energy-limited hydrodynamic blow-off equations to model whether an atmosphere thins toward the Armstrong limit ($6.3\text{ kPa}$) or retains multi-bar envelopes. 
  $$\dot{M}_{\text{loss}} = \frac{\epsilon \pi R_p R_{\text{XUV}}^2 F_{\text{XUV}}(t)}{G M_p K_{\text{tide}}}$$
* **Deep Ocean Hydrostatic Pressure Profiles and Phase Transitions:** Uses depth-integrated hydrostatic pressure ($P(z) = P_0 + \int \rho(z) g\,dz$) and thermodynamic Clapeyron phase diagrams to calculate whether ocean-floor pressures exceed $1\text{--}1.2\text{ GPa}$, precipitating high-pressure ice polymorphs (Ice VI, VII) that mechanically seal the mantle.
* **Fluid Dynamics, Compressibility, and Gas Viscosity:** Utilizes Navier-Stokes equations, Reynolds numbers, and fluid density equations ($\rho = \mu P_0 / R T$) to calculate airway resistance, establishing the limit where diaphragmatic muscular ventilation collapses around $10\text{--}12\text{ MPa}$.

### 2.2 Domain 2: Thermal Dynamics & Climate Modeling

* **Orbital Mechanics and Insolation Geometry:** Applies Kepler’s laws, orbital eccentricity equations, and geometric inverse-square flux decay ($S_0 = L_* / 4\pi a_{\text{m}}^2$) to calculate top-of-atmosphere stellar flux and orbit-averaged insolation.
  $$\langle S \rangle = \frac{1}{2\pi} \int_0^{2\pi} \frac{L_*}{4\pi r(\theta)^2}\,d\theta = \frac{S_0}{\sqrt{1 - e^2}}$$
* **Planetary Energy Balance and 1D Radiative-Convective Equilibrium:** Solves the Stefan-Boltzmann radiation law ($\sigma T_{\text{eq}}^4 = \frac{\langle S \rangle}{4}(1 - A_B)$) paired with multi-stream radiative transfer equations and dry/moist adiabatic lapse rates to generate surface equilibrium temperatures ($T_{\text{surf}}$).
* **Greenhouse Radiative Forcing and Collision-Induced Absorption (CIA):** Requires quantum mechanical vibrational-rotational band models and statistical collision mechanics to quantify infrared absorption, scaling quadratically with pressure ($\tau_{\text{CIA}} \propto P_0^2/g$), modeling greenhouse warming in dense atmospheres.
* **Moist Greenhouse and Runaway Greenhouse Limits:** Solves the Simpson-Nakajima limits on outgoing longwave radiation ($\text{OLR}_{\text{limit}} \approx 280\text{--}310\text{ W m}^{-2}$) to establish the exact stellar flux threshold where oceans boil and thermal runaway occurs.
* **Psychrometric Equations and Wet-Bulb Thermodynamics:** Employs Antoine and Clausius-Clapeyron equations to model saturation vapor pressure and thermodynamic wet-bulb temperature ($T_{\text{wb}}$) equations to compute the critical boundary ($T_{\text{wb}} \le 31^\circ\text{C}$) for evaporative cooling and homeostatic heat rejection.
* **Cryogenic Eutectic Phase Diagrams and Freezing-Point Depression:** Uses colligative thermodynamics, the Pitzer ion-interaction model, and cryohydric phase equilibria to calculate liquid-brine stability curves ($T_{\text{eut}}$) down to $-50^\circ\text{C}$ to $-70^\circ\text{C}$.

### 2.3 Domain 3: Geochemical & Atmospheric Chemistry

* **Chemical Thermodynamics and Mantle Oxygen Fugacity:** Applies Gibbs free energy minimization ($\Delta G = \Delta G^\circ + RT \ln Q$) and chemical equilibrium constants to mineral redox buffers to derive volcanic gas speciation ($\text{CO}/\text{CO}_2$, $\text{H}_2/\text{H}_2\text{O}$) from silicate melts.
* **Atmospheric Photochemical Kinetics and Radical Networks:** Solves systems of non-linear stiff differential equations for atmospheric reaction networks ($d[X]/dt = P_X - L_X$) under stellar ultraviolet spectra, predicting photochemical runaway thresholds for toxic species like carbon monoxide.
* **Dalton's Law of Partial Pressures and Alveolar Gas Diffusion:** Integrates Dalton’s Law ($P_{\text{tot}} = \sum p_i$) and the physiological alveolar gas equation to calculate inspired partial pressures and evaluate hypoxia, hyperoxia, and hypercapnic blood acidosis ($p\text{CO}_2 > 0.93\text{--}5.0\text{ kPa}$).
  $$p_A\text{O}_2 = f_{\text{O}_2}(P_0 - p\text{H}_2\text{O}^*) - \frac{p_A\text{CO}_2}{\text{RER}}$$
* **Thermodynamics of Solutions:** Calculates water activity ($a_w = p / p_0$) using solute activity coefficients and models solution chaotropicity via enthalpy of dissolution to evaluate the mitotic limit ($a_w \ge 0.605$) and cellular lysis limits.
* **Aqueous Equilibrium Chemistry:** Employs carbonic/sulfuric acid-base dissociation equilibria, carbonate saturation states, and carbonate-silicate weathering kinetics to determine ocean pH, phosphorus solubility, and metal mobilization. *(Note: Closed-form negative feedback equations balancing $\text{CO}_2$ drawdown with surface temperature remain an unmodeled data gap requiring external input).*

### 2.4 Domain 4: Ionizing Radiation & Particle Physics

* **Relativistic Particle Motion and Magnetospheric Deflection Dynamics:** Uses the Lorentz force equation ($\mathbf{F} = q(\mathbf{E} + \mathbf{v} \times \mathbf{B})$) and Størmer theory to calculate magnetic rigidity cutoff thresholds ($R_c$) and Chapman-Ferraro magnetopause standoff distances, determining which galactic cosmic rays penetrate planetary magnetic fields.
* **Hadronic Air Showers and Nuclear Spallation Cascades:** Applies relativistic collision kinematics, inelastic hadronic cross-sections, and particle decay lifetimes to calculate Pion ($\pi^0, \pi^\pm$) production, electromagnetic cascades, and relativistic muon propagation down to the surface.
* **Charged Particle Stopping Power and Matter Attenuation:** Implements the Bethe-Bloch formula for continuous ionization energy loss ($dE/dx$) and the Bragg peak, modeling the exponential and ballistic stopping power of atmospheric gas columns, liquid water, and silicate regolith.
* **Radiolytic Reaction Kinetics:** Uses radiolytic chemical yield equations ($G$-values) and linear energy transfer (LET) physics to calculate the generation rate of radiolytic oxidants and reductants that power the Radiolytic Habitable Zone.
* **Statistical Biomacromolecular Degradation:** Applies Manfred Eigen’s quasispecies error threshold formula ($L_{\max} < \ln \sigma / \mu$) and Poisson single-hit target theory to quantify the lethal radiation dose rate ($\dot{D}_{\text{surf}} \ge 10\text{ Gy yr}^{-1}$) that destabilizes nucleic acid polymerization and collapses informational fidelity.

## 3. Mathematical Synthesis: Calculating the PHI-4 Index

Each domain $k \in \{\text{Pressure}, \text{Temperature}, \text{Chemistry}, \text{Radiation}\}$ maps calculated variables through an asymmetrical, continuous trapezoidal-logistic transfer function $S_k(x)$ bounded on $[0, 1]$:
$$S_k(x) = \left[1 + \exp\left(-\kappa_{L} (x - x_{\min})\right)\right]^{-1} \times \left[1 + \exp\left(\kappa_{R} (x - x_{\max})\right)\right]^{-1}$$

Because biological viability across life support systems is strictly non-compensatory, the final habitability index is synthesized using a weighted geometric mean combined with binary geophysical survival gates ($\Theta_j$):
$$\text{PHI-4} = \left( \prod_{k=1}^4 S_k^{w_k} \right)^{\frac{1}{\sum w_k}} \times \prod_{j=1}^m \Theta_j$$
This mathematical construction ensures the habitability score falls continuously to zero if an absolute biophysical limit (e.g., vacuum ebullism, thermal runaway, lethal dose) is breached.

## 4. Algorithmic Pipeline Integration (Pseudocode)

The following basic pseudocode demonstrates the execution architecture of the mathematical framework:

```python
def calculate_phi4_index(system_inputs, bio_profile, weights):
    # --- 1. STOCHASTIC GENERATION & UPSTREAM CONTEXT ---
    star = generate_star(system_inputs.mass, system_inputs.age)
    planet = generate_planet_bulk(system_inputs.radius, system_inputs.cmf, system_inputs.wmf)
    tectonic_state = get_tectonic_regime(planet)

    # --- 2. DETERMINISTIC PHYSICS CALCULATIONS ---
    # Domain 1: Pressure
    gravity = G * planet.mass / (planet.radius ** 2)
    atm_loss = calculate_hydrodynamic_escape(star.xuv_flux, planet)
    P_0 = max((planet.initial_volatiles - atm_loss) * gravity, 0)
    seafloor_pressure = calculate_hydrostatic_ocean_depth(P_0, planet.wmf)

    # Domain 2: Temperature
    insolation = (star.luminosity) / (4 * PI * planet.semi_major_axis ** 2)
    T_surf = compute_radiative_convective_eq(insolation, P_0, planet.albedo)
    runaway_gate = check_simpson_nakajima_limit(insolation, T_surf)
    T_wetbulb = compute_psychrometrics(T_surf, P_0)

    # Domain 3: Chemistry
    mantle_fugacity = compute_gibbs_free_energy(tectonic_state)
    pA_O2, pA_CO2 = alveolar_gas_equation(P_0, planet.gas_fractions)
    water_activity = compute_pitzer_ion_activity(planet.brine_molalities)

    # Domain 4: Radiation
    magnetopause = compute_lorentz_deflection(planet.rotation, star.wind)
    surface_dose = compute_hadronic_cascades(magnetopause, P_0)

    # --- 3. TRANSFER FUNCTIONS & SCORING ---
    S_press = logistic_transfer(P_0, bio_profile.P_min, bio_profile.P_max)
    S_temp  = logistic_transfer(T_wetbulb, bio_profile.Twb_min, bio_profile.Twb_max)
    S_chem  = logistic_transfer(pA_O2, bio_profile.O2_min, bio_profile.O2_max) * \
              logistic_transfer(water_activity, 0.605, 1.0)
    S_rad   = logistic_transfer(surface_dose, bio_profile.Rad_min, bio_profile.Rad_max)

    # --- 4. GEOMETRIC MEAN SYNTHESIS ---
    # Apply Absolute Geophysical Survival Gates (Theta_j)
    if not runaway_gate or seafloor_pressure > 1.2e9:
        return 0.0

    phi4_score = (S_press**weights.w1 * S_temp**weights.w2 * 
                  S_chem**weights.w3 * S_rad**weights.w4) ** (1 / sum(weights.values()))

    return phi4_score
```
