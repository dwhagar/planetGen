> **Imported research, kept as received (banner added 2026-10-15; the text below is unchanged).** It is Boss's document, not a statement of what planetGen does. The design notes that use it: [orbital-solvers-and-integrators](../../orbits/orbital-solvers-and-integrators.md), [orbital-updates](../../orbits/orbital-updates.md), [sky-view](../../orbits/sky-view.md).

Observational Kinematics, Mathematical Physics, and Procedural Generation of Astronomical Rotational Vectors
============================================================================================================

## Observational Kinematics Across Astronomical Populations

The rotational dynamics of celestial objects encompass a broad spectrum of physical regimes, ranging from non-rotating stellar interiors to extreme near-breakup rotators, tidally synchronized exoplanets, and dynamically spun-up compact remnants. Quantifying actual rotational velocities requires distinguishing between the true equatorial velocity ($v_{\text{eq}}$), the projected stellar rotational velocity ($v \sin i$), the rotational period ($P_{\text{rot}}$), and the dimensionless spin parameter ($a^*$). High-resolution spectroscopic surveys and wide-field photometric missions have compiled extensive empirical catalogs that constrain these quantities across stellar, planetary, and compact-object populations (Li et al., 2024; [R4623]).

### Stellar Rotational Velocities and Evolutionary Regimes

Spectroscopic catalogues derived from large-scale astronomical surveys—such as the Gaia Data Release 3 (Gaia DR3), the Apache Point Observatory Galactic Evolution Experiment (APOGEE), and the Large Sky Area Multi-Object Fiber Spectroscopic Telescope (LAMOST)—measure stellar rotation via the line-broadening parameter $v \sin i$, where $i$ represents the inclination of the stellar rotational axis relative to the line of sight (Li et al., 2024). Observed projected rotational velocities span several orders of magnitude, ranging from $v \sin i \approx 0\text{ km/s}$ (observed in stars viewed directly pole-on or in fully despun systems) to values exceeding 300 to 400 km/s in early-type stars.

Massive, hot main-sequence stars (spectral types O, B, and early A) retain a significant fraction of their primordial angular momentum due to the absence of deep outer convective envelopes [R1429]. Equatorial velocities in these early-type stars typically range from 100 to 350 km/s. Extreme sub-classes, such as Be stars, exhibit projected velocities exceeding 300 km/s, rotating at 30% to over 80% of their critical breakup limit ($v_{\text{crit}}$) (Roche, 1873). At these near-critical speeds, centrifugal forces cause severe equatorial oblateness and gravity darkening, triggering circumstellar decretion disk formation (Roche, 1873).

In contrast, low-mass, cool main-sequence stars (late F, G, K, and M dwarfs) possess deep convective envelopes that power magnetic dynamos [R1429, R4544]. Magnetic fields anchor into escaping stellar winds, establishing an efficient magnetic torque that sheds rotational angular momentum over time [R4544]. Consequently, adult solar-type stars are slow rotators, exhibiting equatorial speeds of $v_{\text{eq}} \approx 1.0\text{--}5.0\text{ km/s}$ (the Sun rotates at $v_{\text{eq}} \approx 2.0\text{ km/s}$, corresponding to $P_{\text{rot}} \approx 25\text{--}35\text{ days}$) [R2352, R4544]. The kinematics of these stars within the Galactocentric reference frame are evaluated assuming a local circular velocity $V_{\text{LSR}} \approx 229.0\text{--}232.0\text{ km/s}$ at a galactocentric distance $R_0 \approx 8.2\text{ kpc}$, with a peculiar solar motion $(U_0, V_0, W_0) = (5.0, 17.5, 6.46)\text{ km/s}$ (Li et al., 2024; [R2352]).

As low-mass stars evolve off the main sequence onto the Red Giant Branch (RGB), conservation of angular momentum during massive core collapse and envelope expansion causes surface rotational speeds to decline to $v \sin i < 2.0\text{ km/s}$ (Li et al., 2024). However, high-resolution spectroscopic surveys have identified anomalous populations of rapidly rotating giant stars (Li et al., 2024). Medium-resolution spectra from LAMOST DR9 reveal that among metal-poor giant stars ($[\text{Fe/H}] \in [-4.0, -1.0]$), a critical rotational velocity plateau emerges between 40 and 90 km/s (Li et al., 2024). Stars in this plateau frequently exhibit strong lithium enrichment (abundances ranging from 1.02 to 1.82 dex), indicating that a critical threshold near 40 km/s triggers internal mass transport and activates chromospheric activity as interior material migrates outward (Li et al., 2024).

| **Stellar Class / Population** | **Mass Range (M_sun)** | **Typical v sin i or v_eq (km/s)** | **Typical Rotation Period** | **Angular Momentum Transport Mechanism**                    |
| ------------------------------ | ---------------------- | ---------------------------------- | --------------------------- | ----------------------------------------------------------- |
| O/B Main-Sequence              | 3.0 - 50.0             | 100 - 350                          | 0.5 - 3.0 days              | Weak magnetic braking, high primordial spin retention       |
| Be Stars                       | 3.0 - 20.0             | 250 - 450                          | 0.2 - 1.5 days              | Near-critical equatorial spin-up, decretion disk drives     |
| G-Dwarf (Solar-Type)           | 0.8 - 1.1              | 1.0 - 5.0                          | 20 - 40 days                | Magnetized wind braking (Skumanich spin-down)               |
| M-Dwarf (Fully Convective)     | 0.08 - 0.35            | 0.1 - 15.0                         | 0.5 - 100 days              | Dynamo-saturated magnetic wind braking                      |
| Standard RGB Giants            | 0.8 - 2.0              | 1.0 - 10.0                         | 100 - 1000 days             | Angular momentum conservation during expansion              |
| Metal-Poor RGB Giants          | 0.6 - 1.0              | 40 - 90                            | 10 - 50 days                | Internal mixing, binary engulfment, chromospheric transport |

### Small Solar System Bodies and Planetary Systems

Systematic photometric lightcurve inversions archived in the Asteroid Lightcurve Database (LCDB) show that the rotational periods of minor planets are dictated by body size and structural cohesion [R4623]. Quality codes ($U$) categorize period reliability, where $U = 3$ denotes unambiguous period solutions and $U = 1$ indicates tentative estimates [R4623].

For asteroids larger than $D \approx 200\text{ meters}$, the rotation frequency distribution displays a sharp structural boundary known as the cohesionless spin barrier at $P_{\text{rot}} \approx 2.2\text{ hours}$ (corresponding to a maximum rotation frequency $f \approx 11\text{ rev/day}$) [R4623]. Because bodies larger than 200 m are generally gravity-dominated "rubble piles" held together by self-gravity rather than tensile strength, rotational periods shorter than 2.2 hours generate centrifugal accelerations that exceed self-gravity, causing equatorial mass shedding or structural disruption [R4623]. Super-fast rotators ($P_{\text{rot}} < 2.2\text{ hours}$, such as Kamo'oalewa with $P = 0.465\text{ h}$) exist almost exclusively among small monolithic Near-Earth Asteroids ($D < 200\text{ m}$), where internal material strength withstands centrifugal forces [R4623].

Long-term vector evolution in small bodies is governed by the Yarkovsky-O'Keefe-Radzievskii-Paddack (YORP) effect, a thermal torque resulting from the anisotropic absorption and re-radiation of solar photon momentum [R4568, R4623]. The YORP torque alters both the spin period and the pole orientation $(\lambda, \beta)$, systematically driving orbital obliquities toward asymptotic values near $\beta \approx \pm 90^\circ$ relative to the orbital plane, creating a bi-modal distribution of rotational pole axes (Li et al., 2024; [R4568]).

For close-in planets and natural satellites, rotational velocities are governed by tidal dissipation [R1584]. Differential gravitational forces exert a torque on tidal bulges, transferring spin angular momentum into orbital angular momentum until the rotational period synchronizes with the orbital period ($P_{\text{rot}} = P_{\text{orb}}$) [R1584]. In compact multi-planet systems such as TOI-700, planets TOI-700 b, c, d, and e operate in tidally locked states with zero obliquity relative to their orbital planes.

### Compact Remnant Dynamics and Black Hole Spin Distributions

The rotational state of a black hole is characterized by its dimensionless spin magnitude parameter:

$$a^* = \chi = \frac{c J}{G M^2} \quad \text{where} \quad 0 \le a^* \le 1$$

Gravitational-wave observations by the LIGO-Virgo-KAGRA (LVK) collaboration demonstrate that stellar-mass black hole spins fall into distinct statistical distributions based on their evolutionary assembly channels [R0251]:

First-generation (1G) black holes formed via the direct core collapse of isolated massive stars typically possess low birth spins centered around $a^* \approx 0.1\text{--}0.2$, well-modeled by a Beta distribution $\text{Beta}(\alpha = 1.4, \beta = 3.6)$. High initial spins ($a^* > 0.8$) in isolated 1G black holes occur rarely and require tight binary systems undergoing chemically homogeneous evolution.

Hierarchical merger remnants (2G+) assembled through successive binary coalescences in dense environments (such as globular clusters or AGN accretion disks) display a universal spin distribution [R0251, R2753]. Because the orbital angular momentum of a merging binary converts into the spin of the final remnant, equal-mass non-precessing mergers produce remnants with a characteristic dimensionless spin $a^* \approx 0.69\text{--}0.70$ [R1584]. Subsequent 2G+2G mergers maintain high spin magnitudes spanning $a^* \in [0.6, 0.9]$, subject to retention conditions dictated by gravitational recoil kicks that can reach velocities up to 5000 km/s [R0251].

| **Black Hole Generation** | **Primary Assembly Channel**                  | **Typical Spin Parameter (a*) Distribution** | **Dominant Alignment Characteristics**                         |
| ------------------------- | --------------------------------------------- | -------------------------------------------- | -------------------------------------------------------------- |
| 1G Stellar Collapse       | Isolated binary evolution / Field             | Low spin (mean ~ 0.1 - 0.2), Beta(1.4, 3.6)  | Preferentially aligned with binary orbital axis                |
| 1G Chemically Homogeneous | Tidal synchronization in ultra-tight binaries | High spin (a* ~ 0.8 - 0.95)                  | Strongly aligned with binary orbital plane                     |
| 2G Dynamically Assembled  | Hierarchical mergers in globular clusters     | Moderate-to-high (a* ~ 0.65 - 0.75)          | Isotropically distributed tilt angles due to dynamical capture |
| 2G AGN Disk Remnants      | Gas-driven capture and accretion in AGN disks | Very high spin (a* ~ 0.8 - 0.95)             | Preferentially aligned with galactic disk rotation axis        |

## Mathematical Physics of Rotational Dynamics

### IAU Cartographic Coordinate Systems and Vector Transformation

To define the three-dimensional rotation axis and temporal orientation of an astronomical body, the International Astronomical Union (IAU) Working Group on Cartographic Coordinates parameterizes the direction of the north pole of rotation in Equatorial coordinates relative to the International Celestial Reference Frame (ICRF) via Right Ascension ($\alpha_0$) and Declination ($\delta_0$) [R0135]. The orientation of the prime meridian is tracked as a time-dependent angle $W(t)$ [R0135]:

$$\alpha_0 = \alpha_{\text{ref}} + \dot{\alpha} T$$

$$\delta_0 = \delta_{\text{ref}} + \dot{\delta} T$$

$$W(t) = W_0 + \dot{W} d + \frac{1}{2} \ddot{W} d^2$$

where $T$ represents time in Julian centuries from the J2000.0 epoch ($T = (d - 2451545.0) / 36525$), $d$ represents time in Julian days, $W_0$ is the prime meridian position at epoch, and $\dot{W} = \frac{360^\circ}{P_{\text{rot}}}$ represents the rotational frequency in degrees per day [R0135].

The unit vector specifying the orientation of the rotation axis $\hat{\mathbf{\omega}}$ in Cartesian coordinates is expressed as:

$$\hat{\mathbf{\omega}} = \begin{bmatrix} \cos\delta_0 \cos\alpha_0 \\ \cos\delta_0 \sin\alpha_0 \\ \sin\delta_0 \end{bmatrix}$$

The transformation matrix $\mathbf{R}_{\text{body}\to\text{ICRF}}$ mapping coordinates from the body-fixed coordinate frame $[x_b, y_b, z_b]^T$ to the inertial coordinate system is constructed through sequential Euler angle rotations [R0135]:

$$\mathbf{R}(\alpha_0, \delta_0, W) = \mathbf{R}_z(90^\circ + \alpha_0) \cdot \mathbf{R}_x(90^\circ - \delta_0) \cdot \mathbf{R}_z(W)$$

Expanding this matrix yields the complete transformation tensor:

$$\mathbf{R} = \begin{bmatrix} -\sin\alpha_0 \sin W + \cos\alpha_0 \sin\delta_0 \cos W & -\sin\alpha_0 \cos W - \cos\alpha_0 \sin\delta_0 \sin W & \cos\alpha_0 \cos\delta_0 \\ \cos\alpha_0 \sin W + \sin\alpha_0 \sin\delta_0 \cos W & \cos\alpha_0 \cos W - \sin\alpha_0 \sin\delta_0 \sin W & \sin\alpha_0 \cos\delta_0 \\ -\cos\delta_0 \cos W & \cos\delta_0 \sin W & \sin\delta_0 \end{bmatrix}$$

### Roche Equilibrium Models and Rotational Breakup Limits

For fluid self-gravitating bodies, rapid rotation induces equatorial expansion (Roche, 1873). Under the Roche equilibrium model—which assumes that the internal mass distribution is strongly centrally condensed—the effective gravitational-centrifugal potential $\Phi_{\text{eff}}$ at an equatorial surface distance $r_{\text{eq}}$ is formulated as (Roche, 1873):

$$\Phi_{\text{eff}}(r_{\text{eq}}, \theta = \pi/2) = -\frac{GM}{r_{\text{eq}}} - \frac{1}{2} \omega^2 r_{\text{eq}}^2$$

Critical rotation occurs when the outward centrifugal acceleration matches the inward gravitational acceleration at the equator (Roche, 1873):

$$\omega_{\text{crit}}^2 r_{\text{crit}} = \frac{GM}{r_{\text{crit}}^2} \implies \omega_{\text{crit}} = \sqrt{\frac{GM}{r_{\text{crit}}^3}}$$

Due to rotational distortion, the critical equatorial radius $r_{\text{crit}}$ expands to 1.5 times the polar radius $R_{\text{p}}$ ($r_{\text{crit}} = \frac{3}{2} R_{\text{p}}$) (Roche, 1873). Substituting this geometric limit into the velocity equation determines the critical equatorial velocity $v_{\text{crit}}$ (Roche, 1873):

$$v_{\text{crit}} = \omega_{\text{crit}} r_{\text{crit}} = \sqrt{\frac{GM}{\frac{3}{2} R_{\text{p}}}} = \sqrt{\frac{2}{3} \frac{GM}{R_{\text{p}}}}$$

For solid or cohesionless self-gravitating bodies (such as asteroids), the critical rotational period $P_{\text{crit}}$ below which equatorial shedding occurs depends directly on the mean mass density $\rho$:

$$P_{\text{crit}} = \sqrt{\frac{3\pi}{G\rho}}$$

For a typical asteroid density $\rho = 2000\text{ kg/m}^3$, $P_{\text{crit}} \approx 2.33\text{ hours}$, providing the theoretical origin of the observed 2.2-hour spin barrier [R4623].

### Magnetized Wind Braking and Gyrochronology

Low-mass main-sequence stars shed angular momentum via magnetized stellar winds [R1429, R4114]. The torque $\tau_{\text{wind}}$ exerted on a star by mass loss anchored in its magnetic field lines is expressed as [R1429]:

$$\tau_{\text{wind}} = \frac{d J}{d t} = \frac{2}{3} \dot{M}_* R_{\text{A}}^2 \omega$$

where $\dot{M}_*$ represents the mass-loss rate, $R_{\text{A}}$ represents the Alfvén radius, and $\omega$ represents the stellar angular velocity [R1429]. Parameterizing $R_{\text{A}}$ using magnetohydrodynamic (MHD) formulations yields the spin-down equation [R1429]:

$$\tau_{\text{wind}} = K_1^2 \left(2GM_*\right)^m B_*^{4m} \dot{M}_*^{1-2m} R_*^{5m+2} M_*^{-m} \omega \left(K_2^2 + \frac{\omega^2 R_*^3}{2GM_*}\right)^{-1}$$

where $K_1 = 1.30$, $K_2 = 0.056$, and $m = 0.22$ are dimensionless wind efficiency constants derived from numerical MHD simulations [R1429]. In the unsaturated dynamo regime where magnetic field strength scales with rotation ($B_* \propto \omega$), this torque simplifies to $d\omega / dt \propto -\omega^3$. Integrating this differential equation over time yields the **Skumanich law** [R4114, R4544], where rotational speed decays inversely with the square root of stellar age $t$:

$$\omega(t) \propto t^{-1/2} \implies P_{\text{rot}}(t) \propto t^{1/2}$$

Modern empirical gyrochronology models refine this relation by incorporating color index (such as $B-V$ or $G_{\text{BP}} - G_{\text{RP}}$) to account for mass-dependent convective envelope depths [R1429, R4114]:

$$P_{\text{rot}}(t, B-V) = A \cdot t^n \cdot (B - V - c)^b$$

Fitting this model to open cluster populations yields parameter values $n = 0.55_{-0.09}^{+0.02}$, $b = 0.31_{-0.02}^{+0.05}$, $c = 0.495 \pm 0.010$, and $a = 0.40_{-0.05}^{+0.30}$ [R1429]. Stellar spin behavior is further parameterized by the dimensionless **Rossby number** ($Ro$) [R1429]:

$$Ro = \frac{P_{\text{rot}}}{\tau_{\text{cz}}}$$

where $\tau_{\text{cz}}$ represents the convective turnover timescale [R1429]. For $Ro < Ro_{\text{crit}} \approx 0.1$, stellar activity saturates and spin-down proceeds slowly; for $Ro > 0.1$, stars follow unsaturated Skumanich magnetic braking [R1429].

### Tidal Dissipation and Synchronization Timescales

Gravitational interactions in close binary systems induce tidal bulges that exchange angular momentum between spin and orbital states [R1584]. The rotational synchronization timescale ($\tau_{\text{lock}}$) for a planet or satellite despinning under tidal forces is derived from equilibrium tide theory:

$$\tau_{\text{lock}} = \frac{4}{9} \frac{\omega_0 C Q a^6}{k_2 G M_*^2 R_p^5}$$

where $\omega_0$ is initial spin frequency, $C$ is moment of inertia, $Q$ is the tidal quality factor, $k_2$ is the Love number of degree 2, $a$ is the semi-major axis, $M_*$ is stellar mass, and $R_p$ is planetary radius.

Because $\tau_{\text{lock}}$ scales with semi-major axis to the sixth power ($a^6$), tidal despinning drops precipitously with distance. Bodies within a critical distance despin into synchronized states ($P_{\text{rot}} = P_{\text{orb}}$) on megayear timescales, whereas outer objects retain their primordial rotations over multi-gigayear lifetimes [R1584].

| **Physical Mechanism** | **Governing Master Equation**                                          | **Primary Sensitivity**                       | **Astrophysical Application**                             |
| ---------------------- | ---------------------------------------------------------------------- | --------------------------------------------- | --------------------------------------------------------- |
| Roche Breakup Limit    | v_crit = sqrt((2/3) * G * M / R_p)                                     | Mass density or M^(1/2) R^(-1/2)              | Upper rotational speed cutoff for stars and asteroids     |
| Magnetized Wind Torque | d(omega)/dt = -K * omega^3                                             | Convective turnover timescale                 | Gyrochronology age dating for cool dwarfs                 |
| Tidal Despinning       | tau_lock proportional to (a^6 * Q * M_p) / (k_2 * M_*^2 * R_p^3)       | Orbital separation (a^6)                      | Orbital synchronization of exoplanets and moons           |
| YORP Acceleration      | d(omega)/dt proportional to (Phi_solar / (a^2 * rho * D^2)) * f(shape) | Inversely proportional to surface area (D^-2) | Spin rate evolution and pole alignment of small asteroids |

## Procedural Engine Logic for Simulated Galactic Maps

### Algorithmic Determination of Rotational Speeds

The rotational speed generator evaluates object properties using a conditional physical cascade:

First, the algorithm evaluates tidal synchronization for any orbiting body by calculating $\tau_{\text{lock}}$. If system age $t_{\text{sys}} > \tau_{\text{lock}}$, the object is tidally locked:

$$P_{\text{rot}} = P_{\text{orb}}$$

Its axial obliquity is set to $\varepsilon = 0.0^\circ$ (or locked in a spin-orbit resonance such as 3:2 if orbital eccentricity $e > 0.1$).

Second, for low-mass main-sequence stars ($M_* \le 1.3 M_\odot$) that are not tidally locked, rotational periods are determined via empirical gyrochronology equations [R1429]:

$$P_{\text{rot}}(t, B-V) = \left[ A \cdot t^n \cdot (B - V - c)^b \right] \cdot \exp\left(\mathcal{N}(0, \sigma_{\text{gyro}}^2)\right)$$

where $\mathcal{N}(0, \sigma_{\text{gyro}}^2)$ represents a Gaussian scatter term ($\sigma_{\text{gyro}} \approx 0.10$) accounting for initial angular momentum dispersion [R4114].

Third, for massive hot stars ($M_* > 1.3 M_\odot$), equatorial velocities are sampled from a log-normal distribution:

$$\ln(v_{\text{eq}}) \sim \mathcal{N}(\mu_{v}, \sigma_{v}^2) \quad \text{with } \mu_v = \ln(180\text{ km/s}), \; \sigma_v = 0.45$$

The sampled velocity is checked against the Roche critical breakup limit, applying an upper bound (Roche, 1873):

$$v_{\text{eq, max}} = 0.85 \times \sqrt{\frac{2}{3} \frac{GM_*}{R_{\text{p}}}}$$

Fourth, for small solar system bodies, if diameter $D > 0.2\text{ km}$, $P_{\text{rot}}$ is sampled from a Maxwellian/log-normal distribution truncated at the cohesionless spin barrier $P_{\text{min}} = 2.2\text{ hours}$ [R4623]:

$$P_{\text{rot}} = \max\left(2.2\text{ hours}, \; \text{LogNormal}(\mu = 2.1, \sigma = 0.65)\text{ hours}\right)$$

Fifth, for stellar-mass black holes, 1G black holes draw spins from $a^* \sim \text{Beta}(\alpha = 1.4, \beta = 3.6)$, while hierarchical merger remnants draw spins from $a^* \sim \mathcal{N}(\mu = 0.69, \sigma = 0.04)$.

### Vector Orientations and Axial Tilt Modeling

A body's rotational vector orientation $\hat{\mathbf{\omega}}$ is defined by its axial obliquity $\varepsilon$ relative to its parent frame's orbital normal vector $\hat{\mathbf{L}}_{\text{orbit}}$.

For stars and unperturbed planets, rotational axes align closely with the primordial angular momentum of the local galactic neighborhood:

$$\varepsilon \sim \text{Rayleigh}(\sigma_{\varepsilon} = 15^\circ)$$

For solid planets subject to late-stage collisions, planetary obliquities are sampled from an impact-modified tilt distribution:

$$p(\varepsilon) \propto \frac{1}{2} \sin\varepsilon \left(1 + \gamma \cos^2\varepsilon\right)$$

For asteroids affected by long-term thermal radiation, YORP torques force obliquities toward the ecliptic poles [R4568]:

$$\varepsilon \sim \begin{cases} \mathcal{N}(\mu = 10^\circ, \sigma = 8^\circ) & \text{with probability } 0.5 \\ \mathcal{N}(\mu = 170^\circ, \sigma = 8^\circ) & \text{with probability } 0.5 \end{cases}$$

Given an orbital normal vector $\hat{\mathbf{L}} = [L_x, L_y, L_z]^T$, the final rotational vector $\hat{\mathbf{\omega}}$ is computed combining obliquity $\varepsilon$ and precession angle $\phi \sim U(0, 2\pi)$:

$$\hat{\mathbf{\omega}} = \cos(\varepsilon) \hat{\mathbf{L}} + \sin(\varepsilon) \left[ \cos(\phi) \hat{\mathbf{u}} + \sin(\phi) \hat{\mathbf{v}} \right]$$

The unit vector components $\hat{\mathbf{\omega}} = [\hat{\omega}_x, \hat{\omega}_y, \hat{\omega}_z]^T$ are converted into standard IAU equatorial coordinates $(\alpha_0, \delta_0)$ [R0135]:

$$\delta_0 = \arcsin(\hat{\omega}_z)$$

$$\alpha_0 = \text{atan2}(\hat{\omega}_y, \hat{\omega}_x)$$

| **Object Domain**         | **Procedural Sampling Function**                        | **Target Astrophysics / Observational Anchor**      |
| ------------------------- | ------------------------------------------------------- | --------------------------------------------------- |
| Cool Main-Sequence Spin   | P_rot = A * t^0.55 * (B - V - c)^b * exp(N(0, sigma^2)) | Magnetic wind braking / Gyrochronology              |
| Massive Star Spin         | v_eq = min(v_crit, LogNormal(180 km/s, 0.45))           | High initial angular momentum & Roche limit         |
| Rubble-Pile Asteroid Spin | P_rot = max(2.2 h, LogNormal(13 h, 0.65))               | Cohesionless spin barrier (D > 200 m)               |
| 1G Black Hole Spin        | a* ~ Beta(alpha = 1.4, beta = 3.6)                      | Stellar collapse birth spin parameters              |
| 2G Black Hole Spin        | a* ~ N(mu = 0.69, sigma = 0.04)                         | Universal hierarchical merger dynamics              |
| Stellar Axial Obliquity   | epsilon ~ Rayleigh(sigma = 15 deg)                      | Galactic disk primordial angular momentum alignment |
| Asteroid Pole Obliquity   | epsilon ~ Bi-modal Gaussian(10 deg or 170 deg)          | YORP thermal torque asymptotic pole states          |

### Integrated Procedural Workflow Walk-Through

1. **Galactic Frame Setup**: Define local galactic orbital plane normal:

$$\hat{\mathbf{n}}_{\text{gal}} = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}$$

2. **Primary Solar-Type Star Initialization**:
   
   * Parameters: $M_* = 1.0 M_\odot$, age $t = 4.6\text{ Gyr}$, $B-V = 0.65$.
   
   * Gyrochronology calculation yields $P_{\text{rot}} \approx 25.4\text{ days}$ [R1429], corresponding to $v_{\text{eq}} \approx 2.0\text{ km/s}$.
   
   * Sample axial tilt: $\varepsilon = 7.2^\circ$ drawn from $\text{Rayleigh}(15^\circ)$.

3. **Close-In Terrestrial Exoplanet Initialization ($a = 0.03\text{ AU}$)**:
   
   * Parameters: $M_p = 1.0 M_\oplus$, $R_p = 1.0 R_\oplus$, $Q = 100$, $k_2 = 0.3$.
   
   * Synchronization timescale: $\tau_{\text{lock}} \approx 1.2\text{ Myr} < 4.6\text{ Gyr} \implies$ **Tidally Locked**.
   
   * Rotational state: $P_{\text{rot}} = P_{\text{orb}} = 2.1\text{ days}$, $\varepsilon = 0.0^\circ$, $\hat{\mathbf{\omega}}_p = \hat{\mathbf{L}}_{\text{orbit}}$.

4. **Outer Asteroid Belt Body Initialization ($D = 5\text{ km}, a = 2.5\text{ AU}$)**:
   
   * Cohesionless spin limit applies ($P_{\text{min}} = 2.2\text{ h}$) [R4623].
   
   * Period sampled: $P_{\text{rot}} = 6.4\text{ h}$. YORP tilt sampled: $\varepsilon = 172^\circ$ (retrograde rotator) [R4568].

5. Analytical Conclusions and Implementation Guidelines

Celestial rotational velocities are fundamentally constrained by mechanical breakup limits or despun by magnetic and tidal torques [R4544, R4623]. High-mass stars and monolithic asteroids rotate near physical breakup thresholds ($v_{\text{crit}}vcritv_{\text{crit}}vcrit​$ and $P_{\text{crit}} \approx 2.2\text{ h}Pcrit≈2.2 hP_{\text{crit}} \approx 2.2\text{ h}Pcrit​≈2.2 h$) (Roche, 1873; [R4623]), whereas cool dwarfs despin following Skumanich spin-down laws ($P_{\text{rot}} \propto t^{0.55}Prot∝t0.55P_{\text{rot}} \propto t^{0.55}Prot​∝t0.55$) [R1429, R4114]. Rotational vectors reflect persistent directional torques: close-in exoplanets align into zero-obliquity synchronized states, while thermal YORP torques align asteroid poles perpendicular to orbital planes [R4568]. Procedural engines incorporating these physical relationships generate maps that remain grounded in empirical astrophysics.

## References

* [R0135]

* [R1056]

* Roche, É. (1873). Critical rotation velocity and equilibrium shapes of self-gravitating fluid bodies. _Académie des Sciences et Lettres de Montpellier_.

* [R4114]

* [R4568]

* [R4623]
