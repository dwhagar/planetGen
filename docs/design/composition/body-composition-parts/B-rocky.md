# Part B: rocky and iron planets, mantle mineralogy, equations of state

Agent B, 2026-10-17. Read-only research; nothing in the repo was touched. Cross-checked against `/mnt/project-files/reports/bulk-composition-sources-check-2026-10-17.md` (cited here as G425 s.N, not re-verified). Labels: READ, 2ND-HAND, COMPUTED, SNIPPET, RECALLED, NOT REACHED as in the brief. Everything marked COMPUTED was run in this session with a scipy reference solver (kept in the scratchpad) and the standard-library solver printed in section 3.8; the check numbers in section 3.9 are the output of that code.

Hard-rule framing: Earth-calibrated laws ([R4882: Eq. 2], PREM) are calibrations. The generator uses physics (EOS plus hydrostatic balance plus stoichiometry); Earth, Venus, Mars, Mercury and the Moon appear only in section 4 as tests.

## 1. Scope and interfaces

**Covers.** Dry rock-and-iron bodies (about 0.1 to 20 Mearth), from coreless silicate through Earth-like to iron-rich (super-Mercury), plus carbon-rich (SiC and graphite) planets. Materials and EOS constants; mass-radius law and its inversion; how bulk composition follows the host star; two generator routes with a verified standard-library solver; hand-overs.

**Takes from other parts / the repo.**

| From | Quantity | Unit | Status |
|---|---|---|---|
| Stars (repo) | `[Fe/H]` | dex | stored (schema v95) |
| Stars (part A, `A-stars.md`, read after it landed) | `mg_fe`, `si_fe`, `o_fe`, `c_fe` (dex, new stored fields, A s.3.3); Fe/Si, Mg/Si, C/O, O/Si derived on demand as number and mass ratios (factors 1.9884, 0.8654, 0.7507, A s.2.3). A's section 2.7 also has hinge fits for `[Ca/Fe]` and `[Al/Fe]` (B14) that are not in its stored list | dex, ratios | new in A; B converts mass ratios to molar Fe/Mg, Si/Mg (3.3) |
| Planet draw (repo) | mass, radius (route 1), or mass only (route 2) | Mearth, Rearth | stored |
| Thermal lane | core temperature `T_core` (only a +/-3 percent radius effect, section 2.6) | K | derive or default 4000 K draw range 3000 to 6000 K |
| Part C | water mass fraction, base pressure of a volatile layer `P_base` (for the rock core radius under oceans) | fraction, Pa | C owns; B accepts `P_surf` |

**Hands over.**

| To | Quantity | Unit | Where |
|---|---|---|---|
| C | `CMF` (core mass / dry rock mass), `R_dry(M, CMF, x_fe, core_scale, P_surf)`, rock-core radius under a volatile layer | fraction, Rearth | 3.4, 3.8 (`structure(..., P_surf=P_base)`) |
| D | parent-body mantle oxide composition (source composition of moons, fragments) | mass fraction | 3.3 |
| E | `CMF`, mantle oxide mass fractions MgO, SiO2, FeO, CaO, Al2O3 (computed) and Na2O (parametrised) and K2O (gap), mantle FeO | mass fraction | 3.3; E's section 1 table asks for exactly these |
| E | bulk silicate (mantle plus crust) is not split here; crust types, melting and crust composition are E's (E s.2.2, 2.3) | | |
| Atmosphere lane | surface gravity `g`, escape speed `v_esc`, core mass fraction and mantle FeO (oxidation state), carbon-planet flag, and the statement that B does not derive C, N or H2O inventories (gap) | m/s2, m/s | 3.6, 3.9 |
| Tests | CMF identity, closure, monotonicity, round trip | | 3.10 |

## 2. Findings

### 2.1 Materials and EOS constants

All SLB rows below are 2ND-HAND: the numbers were parsed from the open-source file `burnman/minerals/SLB_2024.py` (burnman 2.1.0, from PyPI) and mean the [R4233] parameter set "Thermodynamics of mantle minerals III: the role of iron" (Stixrude and Lithgow-Bertelloni, GJI 237, 1699, [R4233]); the bibliography was confirmed via Crossref but the OUP full text was refused, so the original table was not read. Molar mass is computed from the formula with IUPAC standard atomic weights as used by burnman. EOS form for every SLB row: Mie-Grueneisen-Debye with 3rd-order Eulerian finite strain (burnman `slb3`), thermal expansion is not a tabulated parameter but follows from gamma, theta and the heat capacity; `q` is the volume exponent of gamma. ρ0 is the 1-bar, 300 K reference density (g/cm3); V0 is per formula unit as burnman defines it.

**Table 2.1a: silicates, oxides, alumina phases ([R4233] via burnman 2.1.0)**

| Phase | Formula (per burnman unit) | M (g/mol) | V0 (cm3/mol) | ρ0 (g/cm3) | K0 (GPa) | K0' | Theta0 (K) | gamma0 | q0 | Stable in mantle at (COMPUTED, 3.1) | Label |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Forsterite (olivine) | Mg2SiO4 | 140.691 | 43.603 | 3.227 | 127.96 | 4.218 | 809.2 | 0.993 | 2.107 | to about 10 GPa (Mg/Si > 1) | 2ND-HAND |
| Fayalite | Fe2SiO4 | 203.772 | 46.290 | 4.402 | 136.49 | 4.882 | 619.0 | 1.084 | 2.881 | | 2ND-HAND |
| Mg-wadsleyite | Mg2SiO4 | 140.691 | 40.515 | 3.473 | 168.70 | 4.123 | 849.1 | 1.205 | 2.208 | about 10 to 20 GPa | 2ND-HAND |
| Fe-wadsleyite | Fe2SiO4 | 203.772 | 42.800 | 4.761 | 168.57 | 4.123 | 636.8 | 1.205 | 2.208 | | 2ND-HAND |
| Mg-ringwoodite | Mg2SiO4 | 140.691 | 39.493 | 3.562 | 184.90 | 4.220 | 879.8 | 1.108 | 2.392 | about 17 to 25 GPa | 2ND-HAND |
| Fe-ringwoodite | Fe2SiO4 | 203.772 | 41.860 | 4.868 | 213.41 | 4.220 | 651.5 | 1.262 | 2.392 | | 2ND-HAND |
| Mg-bridgmanite | MgSiO3 | 100.388 | 24.445 | 4.107 | 250.57 | 4.134 | 893.0 | 1.545 | 0.834 | 25 to 127 GPa | 2ND-HAND |
| Fe-bridgmanite | FeSiO3 | 131.928 | 25.321 | 5.210 | 270.58 | 4.010 | 740.4 | 1.545 | 0.834 | | 2ND-HAND |
| Mg-post-perovskite | MgSiO3 | 100.388 | 24.111 | 4.164 | 247.74 | 4.000 | 931.0 | 1.640 | 1.765 | above about 127 GPa (T about 2530 K) | 2ND-HAND |
| Fe-post-perovskite | FeSiO3 | 131.928 | 25.296 | 5.215 | 247.74 | 4.000 | 769.3 | 1.640 | 1.765 | | 2ND-HAND |
| Ca-perovskite | CaSiO3 | 116.160 | 27.450 | 4.232 | 236.00 | 3.900 | 800.3 | 1.890 | 0.896 | high P | 2ND-HAND |
| Enstatite (orthopyroxene) | Mg2Si2O6 | 200.775 | 62.676 | 3.203 | 107.08 | 7.028 | 812.2 | 0.785 | 3.438 | to 9.4 GPa | 2ND-HAND |
| Ferrosilite | Fe2Si2O6 | 263.855 | 65.941 | 4.001 | 100.55 | 7.876 | 677.9 | 0.714 | 3.438 | | 2ND-HAND |
| HP-clinoenstatite | Mg2Si2O6 | 200.775 | 60.760 | 3.304 | 116.03 | 6.237 | 824.4 | 1.125 | 0.204 | about 9.4 to 17 GPa | 2ND-HAND |
| Clinoenstatite | Mg2Si2O6 | 200.775 | 62.500 | 3.212 | 113.76 | 4.806 | 805.6 | 1.009 | 0.601 | | 2ND-HAND |
| Diopside | CaMgSi2O6 | 216.548 | 66.039 | 3.279 | 113.76 | 4.806 | 782.6 | 1.009 | 0.601 | | 2ND-HAND |
| Mg-akimotoite | MgSiO3 | 100.388 | 26.354 | 3.809 | 210.69 | 5.215 | 929.0 | 1.193 | 2.225 | about 20 to 25 GPa | 2ND-HAND |
| Periclase | Mg4O4 (4 MgO) | 161.216 | 44.976 | 3.584 | 161.14 | 3.908 | 770.9 | 1.450 | 1.549 | with Mg/Si above 1 | 2ND-HAND |
| Wuestite | Fe4O4 (4 FeO) | 287.376 | 49.024 | 5.862 | 160.70 | 4.000 | 427.0 | 1.450 | 1.549 | | 2ND-HAND |
| Wuestite, low spin | Fe4O4 | 287.376 | 43.400 | 6.622 | 199.70 | 4.000 | 493.0 | 1.450 | 1.549 | | 2ND-HAND |
| Quartz | SiO2 | 60.084 | 22.421 | 2.680 | 61.43 | 19.780 | 883.5 | -0.040 | 1.000 | | 2ND-HAND |
| Coesite | SiO2 | 60.084 | 20.657 | 2.909 | 103.54 | 2.901 | 875.2 | 0.290 | 1.000 | | 2ND-HAND |
| Stishovite | SiO2 | 60.084 | 14.017 | 4.286 | 305.84 | 4.029 | 1096.1 | 1.557 | 2.211 | with Mg/Si below 1 | 2ND-HAND |
| alpha-PbO2 SiO2 (seifertite) | SiO2 | 60.084 | 13.670 | 4.395 | 327.16 | 4.017 | 1133.0 | 1.557 | 2.211 | very high P | 2ND-HAND |
| Corundum | Al2O3 | 101.960 | 25.577 | 3.986 | 252.59 | 3.887 | 932.2 | 1.308 | 1.712 | | 2ND-HAND |
| Al-perovskite | Al2O3 | 101.960 | 24.944 | 4.088 | 242.40 | 4.100 | 856.2 | 1.545 | 0.834 | | 2ND-HAND |
| Al-post-perovskite | Al2O3 | 101.960 | 24.500 | 4.162 | 247.74 | 4.000 | 752.0 | 1.865 | 1.765 | | 2ND-HAND |
| Anorthite | CaAl2Si2O8 | 278.204 | 100.610 | 2.765 | 84.09 | 6.734 | 754.5 | 0.385 | 1.000 | crust/shallow only | 2ND-HAND |
| Pyrope | Mg3Al2Si3O12 | 403.123 | 113.080 | 3.565 | 170.24 | 4.111 | 823.2 | 1.014 | 1.422 | | 2ND-HAND |
| Grossular | Ca3Al2Si3O12 | 450.442 | 125.120 | 3.600 | 167.06 | 3.915 | 822.8 | 1.054 | 1.889 | | 2ND-HAND |
| Mg-majorite | Mg4Si4O12 | 401.550 | 114.324 | 3.512 | 165.12 | 4.212 | 822.5 | 0.977 | 1.536 | | 2ND-HAND |
| Spinel | MgAl2O4 | 142.264 | 39.762 | 3.578 | 195.10 | 4.627 | 801.9 | 0.974 | 3.971 | | 2ND-HAND |
| Hematite | Fe2O3 | 159.687 | 30.287 | 5.272 | 204.25 | 4.100 | 653.8 | 1.589 | 2.225 | | 2ND-HAND |
| Magnetite | Fe3O4 | 231.531 | 44.528 | 5.200 | 183.88 | 5.228 | 529.5 | 1.358 | 0.419 | | 2ND-HAND |

Notes. (1) The "stable at" column is the output of my own Gibbs-energy minimisation on MgSiO3 stoichiometry (SLB end-members, burnman) along a 1600 K potential-temperature adiabat: orthoenstatite to HP-clinoenstatite 9.4 GPa; wadsleyite plus stishovite 17 GPa; ringwoodite/akimotoite 20 GPa; bridgmanite 25 GPa; post-perovskite 127 GPa (T about 2530 K), COMPUTED. (2) Periclase K0 161.14 GPa here vs 177 GPa in some older tables: G425 flagged this 177 vs 161.4 difference; the SLB value is carried. (3) Solid-state EOS extrapolates poorly above about 1 TPa for oxides (the code table stops at 10^12.45 Pa = 2.8 TPa).

**Table 2.1b: iron, iron alloys, light-element compounds, carbon**

| Material | Formula | M (g/mol) | ρ0 (g/cm3) | K0 (GPa) | K0' | Theta0 (K) | gamma0 (q) | EOS form | Source | Label |
|---|---|---|---|---|---|---|---|---|---|---|
| epsilon-Fe (hcp), [R4233] | Fe | 55.845 | 8.255 (V0 6.765 cm3/mol) | 165.00 | 4.970 | 379.7 | 2.131 (1.023) | Mie-Grueneisen-Debye, slb3 | [R4233] via burnman | 2ND-HAND |
| epsilon-Fe, Dewaele room-T Vinet | Fe | 55.845 | 8.270 (V0 6.75 cm3/mol) | 163.4 | 5.38 | 417 | 1.875 (q 3.289) | Vinet below 234.4 GPa | [R1752], SEOS | READ |
| epsilon-Fe, high-P | Fe | 55.845 | V0 4.28575 cm3/mol at 234.4 GPa | K 1145.7 at 234.4 GPa | | | | Holzapfel (c0 3.19, c2 -2.40) fitted to DFT | [R1752], SEOS | READ |
| epsilon-Fe, Seager | Fe | 55.845 | 8.300 | 156.2 | 6.08 | | | Vinet | [R1752] Table F.7 | READ |
| epsilon-Fe, Valencia | Fe | 55.845 | 8.300 | 160.2 | 5.82 | | | Vinet | [R1752] Table F.7 | READ |
| epsilon-Fe, Dorn | Fe | 55.845 | 8.341 | 173.98 | 5.297 | | | BM3 | [R1752] Table F.7 | READ |
| epsilon-Fe, Fei | Fe | 55.845 | 8.269 | 172.7 | 4.79 | | | BM3 | [R1752] Table F.7 | READ |
| **hcp-Fe used in the solver** | Fe | 55.845 | 8.2642 | 157.07 | 5.6167 | | | single Vinet fit to the Dewaele-plus-Hakim hybrid, error 0.3 % over 1 to 1000 GPa and 1.9 % to 20 TPa | this work | COMPUTED |
| FeS | FeS | 87.91 | 6.118 | 144.3 | 4.79 | | | Holzapfel ([R3923] refit) | [R1752] Table D.5 | READ |
| FeSi | FeSi | 83.93 | 6.522 | 219.5 | 4.29 | | | Holzapfel ([R3923] refit) | [R1752] Table D.5 | READ |
| Fe0.95O | Fe0.95O | | 5.764 | 152.7 | 4.13 | | | Holzapfel ([R3923] refit) | [R1752] Table D.5 | READ |
| Fe3C | Fe3C | | 7.980 | 289.4 | 3.81 | | | Holzapfel ([R3923] refit) | [R1752] Table D.5 | READ |
| Liquid Fe [R0111] | Fe | 55.845 | 7.019 | 109.7 | 4.66 (K'' -0.043 /GPa) | | | BM4 | `exoplex-1.0.1a1` sdist; original not reached | 2ND-HAND |
| Fe-Ni alloy | Fe(1-y)Ni(y) | | | | | | | no EOS reached; core Ni about 5 wt% (Earth, test target) changes density by under the light-element uncertainty | [R3626] (Earth core 87 wt% Fe, 5 Ni, 2 S, 4 to 7 Si) | READ for the wt%; EOS NOT REACHED |
| Diamond | C | 12.011 | 3.512 (V0 3.42 cm3/mol) | 446.5 | 1.61 | | alpha 4e-6 /K | Holland-Powell tmt; unsuitable above about 100 GPa; use [R1095] EOS and TFD above 1000 GPa | burnman HP_2011_ds62 (2ND-HAND); EOS route per [R2801, R2803] (READ) | 2ND-HAND |
| Graphite | C | 12.011 | (V0 5.3 cm3/mol) | 31.2 | 3.9 | | | Holland-Powell tmt | burnman HP_2011_ds62 | 2ND-HAND |
| Graphite, planet EOS | C | | 2.250 | 33.8 | 8.9 | | | BM3 | [R2801, R2803] ([R4615] fit) as used in my run | 2ND-HAND |
| Graphite to diamond transition | | | | | | | | at about 10 GPa | [R2801, R2803] | READ (paper's own statement) |
| SiC | SiC | 40.10 | 3.220 | 227 | 4.1 | | | BM3 | [R3992: Table 1] via G425 s.2.7 | 2ND-HAND |
| High-pressure ice VII/X, water | H2O | | | | | | | owned by part C | C-volatile-and-giants.md | not mine |

Molar masses of FeS, FeSi, SiC and diamond are sums of standard atomic weights (burnman values for Fe and Si; S 32.06 and C 12.011 RECALLED, trivially checkable). The hcp-Fe thermal behaviour is in 2.6.

**Table 2.1c: other ExoPlex numbers (2ND-HAND from the sdist, original papers not all reached)**

| Quantity | Value | Source | Label |
|---|---|---|---|
| Atomic masses used | Fe 55.845, Mg 24.306, Si 28.0867, O 15.9994, S 32.065, Ca 40.078, Al 26.981 | `exoplex-1.0.1a1` | 2ND-HAND |
| Solver | 10x10 linear mass balance (`get_percents`), potential temperature default 1600 K, mantle FeO 0 to 20 wt% grid, 0.1 <= Si/Mg <= 2, 0.02 <= Ca/Mg <= 0.09, 0.04 <= Al/Mg <= 0.11 | same | 2ND-HAND |
| Water EOS | BM4, V0 18.797 cm3/mol, K0 2.06 GPa, K' 6.29 | same | 2ND-HAND |

### 2.2 Mass-radius relations (COMPUTED with the reference and the stdlib solver, then cross-checked against READ literature)

Method: two layers (hcp-Fe core, mantle table from 2.1), hydrostatic balance in the mass coordinate s = m^(1/3), 4th-order Runge-Kutta, shooting on the central pressure to reach P_surf = 1 bar. The reference version uses scipy; the stdlib code (3.8) agrees to 0.01 percent on the check case. The mantle density is the Fe-free MgSiO3-stoichiometry assemblage on a 1600 K adiabat times (1 + 0.264 x_Fe), x_Fe = Fe/(Fe+Mg); the factor is pressure-independent to +/-0.1 percent (COMPUTED from SLB Fe and Mg end-members).

**Table 2.2a: radius limits (Rearth), x_Fe = 0.11 unless stated, core density factor 1 unless stated**

| M (Mearth) | Coreless, x=0 | Coreless, x=0.25 | CMF 0.325 (x 0.11, scale 0.9) | CMF 0.68 (x 0.05, scale 0.9) | Pure Fe | Label |
|---|---|---|---|---|---|---|
| 0.1 | 0.548 | 0.534 | 0.505 | 0.461 | 0.389 | COMPUTED |
| 0.3 | 0.757 | 0.738 | 0.700 | 0.644 | 0.543 | COMPUTED |
| 1 | 1.077 | 1.048 | 0.994 | 0.911 | 0.768 | COMPUTED |
| 3 | 1.466 | 1.424 | 1.347 | 1.231 | 1.028 | COMPUTED |
| 10 | 2.004 | 1.940 | 1.829 | 1.664 | 1.365 | COMPUTED |
| 20 | 2.363 | 2.281 | 2.144 | 1.941 | 1.569 | COMPUTED |

Radius against CMF at 1 Mearth (x 0.11, scale 0.9): CMF 0, 0.1, 0.2, 0.3, 0.4, 0.5 give R = 1.0644, 1.0431, 1.0216, 0.9996, 0.9770, 0.9534 (about -2.1 percent per 0.1 of CMF); at 5 Mearth CMF 0, 0.2, 0.3, 0.4 give 1.6572, 1.5857, 1.5491, 1.5118 (COMPUTED). Light-element and FeO sensitivities: pure silicate R(1 Mearth) = 1.0774, R(10) = 2.0039 for x = 0 and no core.

**Table 2.2b: independent checks of the mass-radius law**

| Check | Literature (READ) | This work (COMPUTED) | Comment |
|---|---|---|---|
| Pure silicate, [R3800] radius where a rocky planet stops being rocky: R_thresh 1.48 median, 1.59 (95 % bound) Rearth, which Rogers converts to 3.5 and 5.0 Mearth for pure silicate and 4.5 and 6.0 Mearth for Earth-like | [R3800] (READ) | coreless R at 3.5, 5.0 = 1.528, 1.680 (+3.2 %, +5.7 % vs 1.48, 1.59); Earth-like (CMF 0.325, x 0.11, scale 0.9) at 4.5, 6.0 Mearth = 1.499, 1.613 (+1.3 %, +1.4 %) | The Earth-like case agrees to 1.5 percent; the pure-silicate case is 3 to 6 percent larger (reason not isolated: different silicate EOS and temperature treatment) |
| [R3384] rocky fit | R = 1.03 M^0.29 (density above 3.3 g/cm3), M = 0.90 R^3.45 | Earth-like branch is within that scatter | data-driven fit, not a law |
| NRPZ envelopes, [R4467] Eqs. 23 to 28: g_SM = 4.3 + 6.6 R^2.4 m/s2, rho_SM = 4.6 + 1.5 R^2.1 g/cm3, M_SM = 0.25 + 0.87 R^4.2 Mearth (super-Mercury bound); MN-WW (mantle-only, water-world) g = 3.3 + 4.7 R^2.3, rho = 3.6 + 0.91 R^2.1, M = 0.19 + 0.64 R^4.1 (R in Rearth) | [R4467] | At R = 1.0 rho_MN-WW 4.51 vs mine pure silicate 4.25; rho_SM 6.10 vs mine CMF 0.5 = 6.73 (rho_SM corresponds to CMF about 0.42); at R = 1.5: 5.73 and 8.11 vs mine 5.32 (CMF 0), 7.20, 9.13, 14.5 (CMF 0.3, 0.5, 0.8); at R = 2.0: 7.50 and 11.03 vs 6.83, 9.80, 13.14 | Envelope is compatible with my law within the CMF range it covers |
| [R4882: Eq. 2] / PREM closed forms | G425 s.2 | "100 % Fe" column of [R4882] is not pure Fe (agrees with G425) | |
| [R4879] closed forms [R4879]: P_typ = GM^2/(4 pi R^4); P_cmb = P_typ ln(1/CMF); P0 = P_cmb + 1.5 P_typ; P0 about (200 + 600 CMF)(M/Mearth) GPa (Eq. 23); CRF about sqrt(CMF); MoI C = (1/3)(1 + CMF^2/5) | READ | vs my solver (M 0.3 to 10, CMF 0.1 to 0.4): P_cmb within 30 %, Eq. 22 P0 within 70 %, Eq. 23 within 22 %, CRF within 12 %, MoI within 12 % | Use the solver, not the closed forms (COMPUTED) |

**Rocky / non-rocky boundary and scatter.** [R1478] gap 1.5 to 2.0 Rearth (G425). [R3800] R_thresh median 1.48, 95 % bound 1.59 (READ). [R3384] density cut 2.8 to 3.3 g/cm3 separating rocky from volatile-rich (READ). [R4467]: of 85 planets only 9 are nominally rocky; about 20 percent are super-Mercuries and 30 percent volatile-rich; selection-biased sample (READ). [R1752]: rock-iron miscibility is possible above about 1.75 to 1.9 Rearth, which would blur the core boundary (READ). Rule for the generator: a body with radius above the dry limit R_dry(M, CMF=0, x=0) (for example 1.466 at 3 Mearth) cannot be dry rock; hand it to part C (water or envelope). Route 1 returns None in that case (check number 3.9: M 5, R 2.3).

**Compression effect.** Lowering the surface pressure to 1 bar vs a volatile layer base: P_surf = 5 GPa lowers the rock-core radius by 0.94 percent at 1 Mearth and 0.49 percent at 4 Mearth; 10 GPa gives -1.77 and -0.93 percent (COMPUTED). This is the "rock core under an ocean" correction for part C.

### 2.3 Density-to-composition inversion and the degeneracy (COMPUTED)

Mass and radius fix the *total iron fraction* well (when the core is close to pure Fe) but not the split between core and mantle FeO, and they cannot fix the core's light-element content. Ranges over x_Fe in {0, 0.11, 0.25} and core density scale in {1, 0.9, 0.8}, for a reference planet with CMF 0.30, x 0.11, scale 0.9:

| M (Mearth) | R (Rearth) | CMF range | Total-iron range (mass fraction) | Label |
|---|---|---|---|---|
| 0.3 | 0.7042 | 0.215 to 0.396 | 0.310 to 0.396 | COMPUTED |
| 1.0 | 0.9996 | 0.213 to 0.398 | 0.309 to 0.398 | COMPUTED |
| 3.0 | 1.3548 | 0.210 to 0.403 | 0.308 to 0.403 | COMPUTED |
| 8.0 | 1.7436 | 0.206 to 0.409 | 0.307 to 0.409 | COMPUTED |

At core scale 1.0 (pure-Fe core) the total iron fraction is 0.308 to 0.317, a spread of only 0.009. The remaining degeneracy is the core density (light elements, thermal): scale 1 to 0.8 moves total Fe 0.31 to 0.40. FeO-versus-core split moves R by under 1 percent (agrees with [R4467] sections 5 to 6 and [R3545], READ). Consequence: a measured (M, R) gives CMF to roughly +/-0.1 absolute unless the light-element content is known.

Single-parameter sensitivities of the inferred CMF at M = 1, R = 1.0 (COMPUTED, `cmf_from_mr`): x_Fe = 0, 0.11, 0.25 at core scale 0.9 give CMF 0.341, 0.299, 0.238 (about -0.04 per +0.1 in x_Fe); core scale 1.0, 0.9, 0.8 at x_Fe 0.11 give 0.263, 0.299, 0.347 (about +0.04 per -0.1 in scale).

Measurement scatter is separate: radius errors of 5 percent translate to about +/-0.2 in CMF at 1 Mearth (slope 2.1 percent in R per 0.1 CMF), which is why observed CMFs are very uncertain [R3545, R3978].

### 2.4 Stoichiometry: bulk ratios to layers

[R4467] ([R4467], ExoPlex paper), READ:

| Relation | Form | Eq. | Label |
|---|---|---|---|
| Mean mantle formula mass per Mg | mu_bar = sum over X in {Mg, Si, Ca, Al} of (X/Mg)(mu_X + n_X mu_O), n = 1 (Mg, Ca), 1.5 (Al), 2 (Si) | | READ |
| CMF without mantle FeO | CMF = (Fe/Mg) mu_Fe / [(Fe/Mg) mu_Fe + mu_bar] | Eq. 8 | READ |
| Mantle FeO and core Si | reaction 17: 3 Fe + SiO2 -> 2 FeO + FeSi, so core Si moles = FeO moles / 2 | Eqs. 12 to 14, 18 to 22 | READ |
| Core light element mass fraction w_LE | CMF = [1 + (1 - w_LE) mu_bar / mu_Fe / (Fe/Mg)]^-1; core density scales with (mean molar weight)/(mu_Fe), 10 wt% O gives 0.8 | | READ |

The coded version `layers_from_bulk` adds an explicit closed form for the coupled case: e = w K / (mu_FeO - w (mu_FeO - 0.5 mu_SiO2)), where K = mu_MgO + Ca/Mg mu_CaO + Al/Mg mu_AlO1.5 + Si/Mg mu_SiO2 and w is the FeO mass fraction of the mantle (COMPUTED). The ideal-mixing density rule (rho proportional to mean atomic mass of the alloy) is *approximate*: it gives 0.787 for 50 mol% S ([R1752]: 0.8) but 0.751 for 50 mol% Si (Hakim: 0.8) and 0.634 for 23 wt% O (Hakim: 0.8 at 51 mol% = 23 wt% O; the Fe0.95O end-member is 5.764/8.27 = 0.70). So the rule is within 2 to 6 percent for S and Si and too low by about 20 percent for large O contents; the nuisance range 0.8 to 1.0 of the core density factor covers realistic light-element loads (COMPUTED comparison; Hakim READ).

**Sensitivity of CMF to the stellar ratios (COMPUTED, no mantle FeO, Ca/Mg 0.07, Al/Mg 0.09):**

| Si/Mg \ Fe/Mg | 0.4 | 0.55 | 0.71 | 0.9 | 1.2 | 1.6 |
|---|---|---|---|---|---|---|
| 0.6 | 0.208 | 0.266 | 0.318 | 0.372 | 0.441 | 0.513 |
| 0.79 | 0.188 | 0.242 | 0.292 | 0.343 | 0.410 | 0.481 |
| 1.0 | 0.170 | 0.220 | 0.267 | 0.316 | 0.381 | 0.451 |

With 7.7 wt% mantle FeO (Putirka median) and Si/Mg 0.79: Fe/Mg 0.55, 0.71, 0.9 give CMF 0.206, 0.258, 0.312 and core Si 5.79, 4.32, 3.32 wt% (COMPUTED).

### 2.5 Mantle composition and mineralogy

| Item | Value | Source | Label |
|---|---|---|---|
| Mg/Si rule (mass thresholds 0.86 and 1.73 equal molar 1 and 2): below 1 pyroxene plus silica (stishovite/seifertite at depth); 1 to 2 pyroxene plus olivine; above 2 olivine plus MgO | | [R4335] via G425 s.2.10; [R0455] | 2ND-HAND |
| Deep mantle, my Gibbs minimisation | Mg/Si = 0.7: pyroxene then bridgmanite plus stishovite/seifertite; Mg/Si = 1.0: single MgSiO3 phases through the table; Mg/Si = 1.5: fo/wa/ri plus en, then pv plus periclase | this work | COMPUTED |
| Density law with Fe | rho(x) = rho(0)(1 + 0.264 x), x = Fe/(Fe+Mg), pressure-independent to +/-0.1 % | this work | COMPUTED |
| Fe effect factor | 1.029 (x 0.11), 1.066 (x 0.25) | this work | COMPUTED |
| Temperature | 300 K isotherm vs 1600 K adiabat: up to +4.4 % at about 1 GPa and +2.4 % at 127 GPa | this work | COMPUTED |
| Mantle composition difference in total mass: Mg2SiO4 vs MgSiO3/SiO2 mixtures differ by under 2 to 3 %; FeSiO3 gives +32 % / +48 % at Earth-like core size | [R1752: Table 4] | READ |
| Mantle FeO of exoplanets | 2 to 12 wt% (median 7.7) at alpha_Fe = 0.263; alpha_Fe = FeBSP/FeBP ranges: Mercury 0 to 0.12, Earth 0.263 to 0.494, Mars 0.54 to 0.58 | [R3626] | READ |
| Mantle mineralogy classes | (FeO+MgO)/SiO2; mantles dominated by olivine and/or orthopyroxene; median olivine fraction 38 % (alpha 0.263), 60 % (Mars-like), 19 % (Mercury-like) | same | READ |
| Alkalis | Na2O maximum 1.33 to 1.45 wt% if no volatile loss (Earth BSE 0.13 to 0.36) | same | READ |
| Oxygen budget | 97.8 % of 3266 stars have enough O to oxidise Si, Mg, Fe and C | same | READ |
| Planetesimal sources | Solar abundances mirror Earth's to about 10 mol% in Mg, Si, Fe | [R4470] ([R4470], abstract) | READ |
| Fragment populations | polluted white dwarfs show core-rich, mantle-rich and crust-rich fragments | Xu and [R0458, R4797] | READ (qualitative) |
| Coreless planets | water-oxidised iron into mantle FeO | [R1248] ([R1248], abstract) | READ |

### 2.6 Thermal effect on the core

Density ratio rho(P, T)/rho(P, 300 K) of hcp-Fe from [R4233] `fee` (COMPUTED with burnman):

| P (GPa) | 1000 K | 2500 K | 4000 K | 6000 K | 8000 K |
|---|---|---|---|---|---|
| 135 | 0.9923 | 0.9722 | 0.9475 | 0.9055 | 0.8465 |
| 330 | 0.9962 | 0.9862 | 0.9747 | 0.9568 | 0.9355 |
| 500 | 0.9973 | 0.9902 | 0.9821 | 0.9698 | 0.9556 |
| 1000 | 0.9986 | 0.9945 | 0.9900 | 0.9833 | 0.9759 |
| 2000 | 0.9992 | 0.9969 | 0.9944 | 0.9907 | 0.9867 |
| 3000 | 0.9995 | 0.9978 | 0.9960 | 0.9934 | 0.9906 |

Fit in code: deficit = 0.0525 ((T - 300)/3700)(135/P_GPa)^0.84, about +/-15 percent of the deficit. Values at 50 GPa or below with high T are not physical (solid-state extrapolation); the code uses it only above the CMB pressure.

hcp-Fe density cross-comparison (g/cm3, COMPUTED): at P = 50 / 135 / 234.4 / 330 / 500 / 1000 / 2000 / 5000 GPa: SLB 300 K 9.973 / 11.654 / 12.998 / 14.006 / 15.430 / 18.368 / 22.143 / 28.756; SLB 4000 K 8.840 / 11.042 / 12.562 / 13.651 / 15.154 / 18.184 / 22.019 / 28.681; Seager Vinet 9.987 / 11.552 / 12.801 / 13.744 / 15.093 / 17.951 / 21.803 / 29.107; Dewaele Vinet 9.977 / 11.653 / 13.027 / 14.081 / 15.611 / 18.927 / 23.530 / 32.607; Hakim SEOS 10.114 / 11.738 / 13.030 / 14.022 / 15.471 / 18.646 / 23.086 / 31.695. At 10 TPa: Seager 37.131, Dewaele 43.013, Hakim 41.030 (SLB not defined). Vinet extrapolation departs from SEOS by more than 5 % at 10 TPa and BM3 must not be used to extrapolate ([R1752], READ). The thermal correction for a 1 Mearth planet (4000 K) is about 3 to 5 percent in core density, 1 percent in R.

### 2.7 Carbon planets (COMPUTED radii, READ criteria)

| Item | Value | Source | Label |
|---|---|---|---|
| Equilibrium C/O (number) above which solid carbon (SiC plus graphite) replaces silicates | 0.8 (at 1e-4 bar) | [R0455] / Thiabaud via G425 s.2.10 | 2ND-HAND |
| Sequential condensation | carbon-rich planetesimals as low as C/O about 0.65 | no source | READ |
| Host-star distribution (849 FGK dwarfs) | C/O median 0.47, sd 0.07, max 0.66 (+/-0.068), Sun 0.55; C/O > 1 in under 0.13 % ; Mg/Si median 1.02, 40 % below 1, 60 % in 1 to 1.4; C/O rises with [Fe/H] and plateaus above about -0.1 dex | [R0535, R0536] | READ |
| Older data | mean Mg/Si 1.32, C/O 0.77 for hosts; 21/60 above 0.8; with error propagation only 5 of 61 (8 %) at 2 sigma and 7 (11 %) at 1 sigma | [R0455] | 2ND-HAND (via G425) |
| EOS | graphite below 10 GPa, diamond above [R1095], TFD above 1000 GPa; SiC/MgSiO3/Fe from [R3992] | [R2803] | READ |
| Radii (Rearth) at 1 / 5 Mearth, Seager-constant BM3 with Fe Vinet core | SiC coreless 1.146 / 1.844; SiC CMF 0.2: 1.088 / 1.741; CMF 0.325: 1.049 / 1.673. Graphite coreless 1.236 / 1.976; CMF 0.2: 1.167 / 1.857; CMF 0.325: 1.120 / 1.779. MgSiO3 pv BM3 coreless 1.047 / 1.662; CMF 0.2: 1.001 / 1.583; CMF 0.325: 0.971 / 1.533 | this work | COMPUTED |

Carbon planets with C/O above the cut are rare in the solar neighbourhood (Brewer and Fischer: nothing above 0.66 in 849 stars), so the carbon branch should almost never fire in a sample drawn from solar-neighbourhood statistics. The Seager MgSiO3 comparison row is stiffer than the SLB assemblage (0.971 vs 0.994 at CMF 0.325), which is the sort of EOS-choice spread (about 2 percent) to expect.

### 2.8 Star and planet composition link

| Study | Result | Label |
|---|---|---|
| [R0032] | f_iron(star, Sun) = 33.2 +/- 1.7 %; f_core(planet) = 6.8 +/- 1.3 f_star - 184 +/- 40 (full sample, 22 rocky planets, 5 super-Mercuries); super-Earths only 4.5 +/- 0.8 (-115); rho/rho_Earth-like vs f_star slope 0.074 +/- 0.012 (-1.4) | READ |
| [R3545] | stars Fe/Si 1.63 (+0.91/-0.85) vs planets 1.15 (+1.43/-0.76); CMF stars 0.32 (+0.14/-0.12) vs planets 0.24 (+0.33/-0.18); max enrichment rho0 about 6 g/cm3, CMF about 0.8 | READ |
| [R3978] | 10 of 11 planets (91 %) consistent with host at 2 sigma; mean CMF_star 0.28 (-0.06) | READ |
| [R0541] | 20 planets, 1 to 1.8 Rearth, homogeneous stellar data; OLS slope m = 1.3 +/- 1.0 (b = -0.05 +/- 0.28); ODR slope m = 5.6 +/- 1.6 (b = -1.4 +/- 0.5); 75 % of planets within 1 sigma of host CMF; no distinct super-Mercury population; "current uncertainties and modelling differences prevent definitive conclusions" | READ (abstract and sections 5, 6.1.1) |
| Plotnykov et al. ([R3546] v2, January 2026) | 21 planets, 20 stars, homogeneous data; best-fit line slope m = 0.94 (+1.02/-1.07), intercept c = -0.02 (+0.31/-0.29) for planet Fe mass fraction vs star-equivalent; consistent at 1 sigma with other homogeneous studies, not with heterogeneous samples (steeper); super-Mercuries reduced from four to one; zero slope also allowed at 1 sigma; outliers recovered in mock populations about 15 % of the time; 150 planets needed to settle; planet Fe-MF vs [Fe/H]: m = 0.24 +/- 0.16, c = 0.23 +/- 0.03 | READ (abstract and section 3.1) |
| Hypatia ([R4467] s.3) | molar Fe/Mg = 0.71 +/- 0.18, Si/Mg = 0.79 +/- 0.19 (correlated, correlation not given); mean CMF 0.29 (mu_bar 97.2) | READ |

Reading. The steep slope (above 4) of [R0032] depends on a heterogeneous sample and on the fitting method (OLS 1.3 vs ODR 5.6 on the same data, [R0541]); the cleaner homogeneous sample gives a slope consistent with 1. So the default must be a 1:1 relation with scatter, with the steep-slope super-Mercury branch as an optional low-probability draw. 

Brewer's Mg/Si median 1.02 (Si/Mg about 0.98) and A's BF16-based number median 1.07 both disagree with Hypatia's Si/Mg 0.79 (Mg/Si 1.27) by 15 to 20 percent; the effect on CMF is only about 0.02 (COMPUTED, table in 2.4).

## 3. How the generator should compute or draw it

### 3.1 Recommendation

Primary design: **route 2** (draw composition first, derive density and radius). Reason: it keeps mass, radius, density, gravity, CMB pressure and moment of inertia mutually consistent by construction, ties the planet to its star through Fe/Mg and Si/Mg, and gives the other lanes (E, atmosphere) a composition that was not back-fitted. **Route 1** stays available as the fill-in for already-drawn (M, R): it returns the composition with an honest uncertainty, because (M, R) alone fix the total iron fraction but not the split (2.3). The galaxy's existing class draws use density; route 1 is the migration path.

### 3.2 Stored quantities (per body)

| Name | Unit | Range / default | Per | Status |
|---|---|---|---|---|
| `mass` | Mearth | 0.1 to 20 | body | stored |
| `radius` | Rearth | from structure | body | stored (derived in route 2) |
| `cmf` | mass fraction of dry rock | 0 to 0.8 | body | **new** |
| `mantle_feo_wt` or `x_fe` | mass fraction / mole fraction | 0 to 0.20 (median 0.077, Putirka) | body | **new** |
| `core_w_si`, `core_w_o`, `core_w_s` | mass fraction of core | 0 to 0.10 each (parametrisation, see 3.4) | body | **new** |
| `core_scale` | density factor | 0.8 to 1.0 (computed from the above and T) | body | **new, derived** |
| `mg_si`, `si_mg`, `fe_mg`, `ca_mg`, `al_mg` | molar | from star via A's `mg_fe`, `si_fe`, [Fe/H] (Mg/Si 0.7 to 1.5 typical) | body (or derive from the star on demand) | **new in B** or derived |
| `c_o_flag` | bool | C/O above 0.8 (or 0.65) | body | **new** |
| `t_core` | K | 3000 to 6000 | body | derive/draw |
| `na2o_retention` | fraction of the 1.33 to 1.45 wt% cap | 0 to 1 | body | **new, parametrised** |
| `mantle_oxides` | mass fractions MgO, SiO2, FeO, CaO, Al2O3, Na2O, K2O | sums to 1 | body | **new, derived** |

### 3.3 Route 2 (composition to structure)

1. Star: take the molar ratios from part A's derived mass ratios: Fe/Mg = (Fe/Si_mass / 1.9884) / (Mg/Si_mass / 0.8654) and Si/Mg = 0.8654 / Mg/Si_mass; C/O (number) = C/O_mass / 0.7507 (conversion factors from A s.2.3, atomic weights Fe 55.845, Si 28.0855, Mg 24.305, C 12.011, O 15.999; the same weights as in `AMU`). A draws `mg_fe`, `si_fe`, `o_fe`, `c_fe` from [Fe/H] with branch-dependent scatter (A s.2.7), so B needs no draw of its own for the star. Ca/Mg and Al/Mg: if A adds `[Ca/Fe]` and `[Al/Fe]` to its draw (its 2.7 hinge fits exist), convert the same way; otherwise draw uniform inside the ExoPlex grid limits Ca/Mg 0.02 to 0.09 and Al/Mg 0.04 to 0.11 (2ND-HAND, PARAMETRISATION). The CMF this gives for A's medians is in the table below (COMPUTED); it is a weak function of [Fe/H] (0.20 to 0.30 over [Fe/H] -0.8 to +0.5).
   **CMF implied by the star, using A's median ratios by [Fe/H] bin (A s.2.8 table; no mantle FeO, Ca/Mg 0.07, Al/Mg 0.09; COMPUTED):**

   | [Fe/H] bin | Fe/Si mass | Mg/Si mass | Fe/Mg molar | Si/Mg molar | CMF (no FeO) | CMF (7.7 wt% FeO) | Label |
   |---|---|---|---|---|---|---|---|
   | -1.0 to -0.6 | 1.18 | 1.31 | 0.392 | 0.661 | 0.198 | 0.161 | COMPUTED |
   | -0.6 to -0.3 | 1.40 | 1.17 | 0.521 | 0.740 | 0.238 | 0.202 | COMPUTED |
   | -0.3 to 0 | 1.63 | 1.03 | 0.689 | 0.840 | 0.279 | 0.245 | COMPUTED |
   | 0 to 0.3 | 1.70 | 0.93 | 0.796 | 0.931 | 0.298 | 0.265 | COMPUTED |
   | 0.3 to 0.6 | 1.70 | 0.93 | 0.796 | 0.931 | 0.298 | 0.265 | COMPUTED |
   | thin-branch sample | 1.68 | 0.97 | 0.754 | 0.892 | 0.291 | 0.258 | COMPUTED |
   | BF16 local medians | 1.81 | 0.93 | 0.847 | 0.931 | 0.311 | 0.279 | COMPUTED |

   A's local Mg/Si (mass 0.93, number 1.07) agrees with Brewer and Fischer's number median 1.02 and differs from Hypatia's Si/Mg 0.79 (number Mg/Si 1.27); the CMF effect of that choice is about 0.02.
2. Planet vs star: take the planet Fe/Mg and Si/Mg equal to the star's (1:1, Thiabaud/Schulze/[R3546]), plus an optional intrinsic spread; with probability p_SM (PARAMETRISATION, 0.05 from [R3546], 1 of 21; up to 0.23 from [R0032], 5 of 22 heterogeneous) switch to a super-Mercury branch with Fe/Mg multiplied by a factor drawn so that CMF is uniform in 0.5 to 0.8 ([R3545]: max CMF about 0.8). A coreless branch (CMF 0 with high mantle FeO) is available for water-rich formation [R1248], tied to part C's water fraction.
3. Mantle FeO: draw `w_feo` from a distribution with median 0.077 and range 0.02 to 0.12 [R3626] as a PARAMETRISATION (a clipped lognormal). If `w_feo` would use more Fe than exists (e > Fe/Mg) redraw.
4. Light elements: core Si from reaction 17 (set by `w_feo`); `core_w_o`, `core_w_s` from 0 to 0.10 each (PARAMETRISATION; Earth's core 2 wt% S, 4 to 7 wt% Si and unknown O are test values only). `core_scale_le` from the ideal-mixing rule (2.4).
5. `layers_from_bulk` gives `cmf`, `x_fe`, core alloy, mantle oxide mass fractions, total Fe.
6. Thermal: `t_core` default 4000 K; two passes of `thermal_deficit` at the CMB pressure.
7. `structure(M, cmf, x_fe, core_scale)` gives radius, central and CMB pressure, g, v_esc, MoI factor.

### 3.4 Route 1 (keep drawn M and R, fill composition)

1. Check the body is dry rock: R must lie between the pure-iron and coreless limits for its mass (2.2a). Otherwise hand to part C.
2. Draw nuisance parameters: mantle FeO uniform 0 to 0.20 wt fraction, core density factor uniform 0.8 to 1.0 (PARAMETRISATIONS spanning 2.3).
3. For each draw invert `radius_fit` for CMF (bisection); keep draws inside [0, 0.9].
4. Report the median CMF and the 16 to 84 percent range; draw the stored value from that sample. Total iron fraction: median and range.
5. Optional anchoring: when the star's Fe/Mg and Si/Mg are known, weight the draws by a Gaussian centred on the star-implied CMF (2.4 table) with width 0.1 (PARAMETRISATION, of the order of the planet-minus-star spread in [R3545]); route 1 then becomes route 2 with a measured-radius likelihood.

### 3.5 Mg/Si and mantle mineralogy rule

| Mg/Si (molar) | Mantle assemblage at the surface | Deep mantle | Note |
|---|---|---|---|
| below 1 | pyroxene plus silica | bridgmanite plus stishovite/seifertite | silica excess |
| 1 to 2 | pyroxene plus olivine | wadsleyite/ringwoodite then bridgmanite plus periclase | |
| above 2 | olivine plus MgO | | rare |

Carbon-planet rule: if C/O (number) is above 0.8 (equilibrium; 0.65 for sequential condensation) set `c_o_flag`; the body uses a SiC and graphite mantle with no silicate (radius table 2.7).

### 3.5a Pseudo-code

```text
route2(star, M, u):                        # u() = seeded uniform in [0,1)
    fe_mg, si_mg, ca_mg, al_mg = star ratios (or drawn from Hypatia/Brewer, section 3.3)
    if u() < p_SM: fe_mg *= factor so that CMF ~ U(0.5, 0.8)   # optional branch
    w_feo   = clipped_lognormal(median 0.077, range 0.02..0.12)
    w_o, w_s = U(0, 0.10), U(0, 0.10)
    lay = layers_from_bulk(fe_mg, si_mg, ca_mg, al_mg, w_feo, w_o, w_s)
    if lay is None: redraw w_feo                       # not enough Fe
    if C/O > 0.8: use carbon mantle (SiC + graphite), section 2.7
    st  = structure(M, lay.cmf, lay.x_fe, lay.core_scale_le * (1 - thermal_deficit(P_cmb, T_core)))
    store cmf, x_fe, core alloy, mantle oxides, R, g, v_esc, MoI, P_cmb

route1(M, R, u):
    if R > R_coreless(M) or R < R_pureFe(M): return NOT_DRY_ROCK  # to part C
    for k in 1..200: draw w_feo ~ U(0, 0.20), core_scale ~ U(0.8, 1.0)
        cmf_k = invert radius_fit for CMF at (M, R, x_from_feo_wt(w_feo), core_scale)
    store median and 16..84 percent range; draw the stored CMF from the sample
```

### 3.6 Derived surface quantities

g = GM/R^2; v_esc = sqrt(2GM/R); MoI factor C/(MR^2) from the integrated profile; P_cmb from the same integration. Handover: g, v_esc and (`cmf`, `mantle_feo_wt`) to the atmosphere lane. The volatile inventory (water, carbon, nitrogen) is not derived here: see the gap list.

### 3.7 Interface to part C

C takes `CMF` as the iron share of the dry part. B's `structure(M, cmf, x_fe, core_scale, P_surf=P_base)` returns the rock-core radius under a volatile layer. C's dry-fit `R = (1.0688 - 0.2317 CMF (1 - x)) M^(1/3.7) (1 + 0.4492 x - 0.1597 x^2)` (x the water fraction) at CMF 0.325, x 0 gives 0.993; B's solver at CMF 0.325, x_Fe 0.11, scale 0.9 gives 0.994 at 1 Mearth. They agree to 0.1 percent at that point; at the limits (coreless x=0, 1 Mearth) C's fit gives 1.0688 vs my 1.077 (+0.8 %). Ocean planets are treated by B only up to the rock-core boundary.

### 3.8 Verified code (standard library only)

The listing below is the exact file that produced the check numbers in 3.9. Core: Vinet fit; mantle: 19-node log-log table with the Fe factor; solver: RK4 in s = m^(1/3) with Illinois shooting on the central pressure (60 core steps, 160 mantle steps). Validation, all COMPUTED: (a) against the scipy reference built from the full burnman mantle tables (90 cases: x_Fe 0, 0.11, 0.25; M 0.1 to 20 Mearth; CMF 0 to 0.8): max radius error 0.20 percent, max MoI error 0.28 percent, max P_cmb error 0.93 percent, max central-pressure error 1.9 percent; (b) the fast fit `radius_fit` against `structure` on 80 fresh random draws over M 0.1 to 20, CMF 0 to 0.8, x 0 to 0.3, scale 0.8 to 1: max error 0.40 percent, rms 0.10 percent.

```python
"""Two-layer rocky-planet structure, standard library only.  Units: SI.
Core  : single Vinet fit to hcp-Fe ([R1096] Vinet below 234.4 GPa + [R1752] SEOS above).
Mantle: Fe-free MgSiO3-composition (Mg/Si = 1) density along a 1600 K potential-temperature adiabat,
        from a Gibbs-energy minimisation with [R4233] end-members; (1 + 0.264 x_Fe) scales it
        for Fe/(Fe+Mg) = x_Fe.  Tables are COMPUTED (see section 3 of the part file)."""
import math
from bisect import bisect_right

G = 6.6743e-11            # [R4361], m3 kg-1 s-2
M_EARTH = 5.9722e24       # kg  (unit conversion only)
R_EARTH = 6.371e6         # m   (unit conversion only)

# --- core: Vinet, rho0 [kg/m3], K0 [Pa], K0' ; fit range 1-1000 GPa, max error 0.3 % in density
CORE = (8264.2, 157.07e9, 5.6167)

def core_density(P, scale=1.0):
    """Vinet inverse by bisection on x = (rho0/rho)^(1/3); scale = light-element density factor."""
    r0, K0, Kp = CORE
    if P <= 0.0:
        return r0 * scale
    lo, hi = 0.05, 1.0
    for _ in range(60):
        x = 0.5 * (lo + hi)
        p = 3.0 * K0 * (1.0 - x) / (x * x) * math.exp(1.5 * (Kp - 1.0) * (1.0 - x))
        if p > P:
            lo = x
        else:
            hi = x
    return r0 / (0.5 * (lo + hi)) ** 3 * scale

# --- mantle table: log10 P [Pa] -> density [kg/m3] at x_Fe = 0, Mg/Si = 1, 1600 K adiabat
_LP = [5.5, 8.0, 9.0, 9.4, 9.7, 9.95, 10.1, 10.25, 10.35, 10.45, 10.6, 10.8, 11.0, 11.15, 11.3,
       11.6, 11.9, 12.2, 12.45]
_RHO = [3055, 3060, 3097, 3154, 3234, 3361, 3457, 3880, 4049, 4351, 4494, 4748, 5091, 5534, 5951,
        7042, 8560, 10631, 12902]
_LR = [math.log(r) for r in _RHO]
_SLOPE_TOP = (_LR[-1] - _LR[-2]) / (_LP[-1] - _LP[-2])

def mantle_density(P, x_fe=0.0):
    lp = math.log10(max(P, 1.0))
    if lp <= _LP[0]:
        lr = _LR[0]
    elif lp >= _LP[-1]:
        lr = _LR[-1] + _SLOPE_TOP * (lp - _LP[-1])
    else:
        i = bisect_right(_LP, lp) - 1
        t = (lp - _LP[i]) / (_LP[i + 1] - _LP[i])
        lr = _LR[i] + t * (_LR[i + 1] - _LR[i])
    return math.exp(lr) * (1.0 + 0.264 * x_fe)

def _rk4_layer(y, s0, s1, nsteps, rho_fn):
    """y = [r, P, C]; independent variable s = m**(1/3) so that dm = 3 s^2 ds."""
    h = (s1 - s0) / nsteps
    def f(s, r, P):
        rho = rho_fn(max(P, 1.0))
        return (3.0 * s * s / (4.0 * math.pi * r * r * rho),
                -3.0 * G * s ** 5 / (4.0 * math.pi * r ** 4),
                2.0 * r * r * s * s)
    s = s0
    r, P, C = y
    for _ in range(nsteps):
        k1 = f(s, r, P)
        k2 = f(s + h / 2, r + h / 2 * k1[0], P + h / 2 * k1[1])
        k3 = f(s + h / 2, r + h / 2 * k2[0], P + h / 2 * k2[1])
        k4 = f(s + h, r + h * k3[0], P + h * k3[1])
        r += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        P += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        C += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        s += h
    return [r, P, C]

def _run(Pc, M, cmf, x_fe, core_scale, n_core, n_mant):
    rho_c = core_density(Pc, core_scale)
    s0 = (1e-6 * M) ** (1.0 / 3.0)
    r0 = s0 * (3.0 / (4.0 * math.pi * rho_c)) ** (1.0 / 3.0)
    P0 = Pc - (3.0 * G / (8.0 * math.pi)) * (4.0 * math.pi * rho_c / 3.0) ** (4.0 / 3.0) * s0 * s0
    y = [r0, P0, 0.0]
    sc = (cmf * M) ** (1.0 / 3.0)
    Pcmb = rcore = None
    if cmf > 1e-9:
        y = _rk4_layer(y, s0, sc, n_core, lambda p: core_density(p, core_scale))
        Pcmb, rcore = y[1], y[0]
        s_from = sc
    else:
        s_from = s0
    if cmf < 1.0 - 1e-9:
        y = _rk4_layer(y, s_from, M ** (1.0 / 3.0), n_mant, lambda p: mantle_density(p, x_fe))
    return y, Pcmb, rcore

def structure(M_earth, cmf, x_fe=0.11, core_scale=1.0, P_surf=1e5, n_core=60, n_mant=160):
    """Radius, central/CMB pressure, g, v_esc and C/(M R^2) for a core (mass fraction cmf) + mantle planet.
    Shooting on log central pressure by the Illinois (regula falsi) method."""
    M = M_earth * M_EARTH
    res = lambda lpc: _run(math.exp(lpc), M, cmf, x_fe, core_scale, n_core, n_mant)[0][1] - P_surf
    a, b = math.log(1e8), math.log(2e15)
    fa, fb = res(a), res(b)
    side = 0
    for _ in range(80):
        c = (a * fb - b * fa) / (fb - fa)
        fc = res(c)
        if abs(fc) < 2e3 or abs(b - a) < 1e-9:
            break
        if fc * fb > 0:
            b, fb = c, fc
            if side == -1:
                fa *= 0.5
            side = -1
        else:
            a, fa = c, fc
            if side == 1:
                fb *= 0.5
            side = 1
    (r, P, C), Pcmb, rcore = _run(math.exp(c), M, cmf, x_fe, core_scale, n_core, n_mant)
    return dict(R=r / R_EARTH, Pc=math.exp(c), Pcmb=Pcmb, Rcore=rcore, g=G * M / r ** 2,
                v_esc=math.sqrt(2.0 * G * M / r), moi=C / (M * r * r), rho_bar=M / (4.0 / 3.0 * math.pi * r ** 3))


# ---------------------------------------------------------------------------
# Fast radius: ln(R/R_earth) = sum_ij a_ij u^i c^j + x(b0 + b1 c + b2 u + b3 u c)
#                              + d c (e0 + e1 c + e2 u + e3 u c),  u = ln(M/M_earth), c = CMF,
# x = Fe/(Fe+Mg) of the mantle, d = 1 - core_scale.  Fitted to structure() (1500 random draws);
# valid 0.1 <= M <= 20, 0 <= CMF <= 0.8, 0 <= x <= 0.3, 0.8 <= core_scale <= 1.
_A = ((0.075176, -0.23494, -0.033795, -0.044414),
      (0.28771, -0.0099019, 0.0021907, -0.0049603),
      (-0.0052534, -0.0041485, 0.0081412, -0.009653),
      (-0.0011712, 0.0006952, -0.0026943, 0.0033964))
_B = (-0.11559, 0.086052, -0.0067087, 0.0040513)
_E = (0.25003, 0.18349, 0.022991, 0.013879)

def radius_fit(M_earth, cmf, x_fe=0.11, core_scale=0.9):
    u = math.log(M_earth)
    s = sum(_A[i][j] * u ** i * cmf ** j for i in range(4) for j in range(4))
    s += x_fe * (_B[0] + _B[1] * cmf + _B[2] * u + _B[3] * u * cmf)
    s += (1.0 - core_scale) * cmf * (_E[0] + _E[1] * cmf + _E[2] * u + _E[3] * u * cmf)
    return math.exp(s)

def cmf_from_mr(M_earth, R_earth, x_fe=0.11, core_scale=0.9, lo=0.0, hi=0.9):
    """Invert radius_fit for the core mass fraction by bisection (R falls as CMF rises).
    Returns None when (M, R) is outside [lo, hi] (denser than iron-rich or lighter than coreless)."""
    f = lambda c: radius_fit(M_earth, c, x_fe, core_scale) - R_earth
    if f(lo) < 0.0 or f(hi) > 0.0:
        return None
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

# ---------------------------------------------------------------------------
# Stoichiometry: bulk molar ratios -> core mass fraction, mantle oxides, core alloy.
AMU = {"Fe": 55.845, "Mg": 24.305, "Si": 28.0855, "O": 15.999, "Ca": 40.078, "Al": 26.9815, "S": 32.06}
MU_MGO = AMU["Mg"] + AMU["O"]                 # 40.304
MU_SIO2 = AMU["Si"] + 2 * AMU["O"]            # 60.085
MU_CAO = AMU["Ca"] + AMU["O"]
MU_ALO15 = AMU["Al"] + 1.5 * AMU["O"]
MU_FEO = AMU["Fe"] + AMU["O"]                 # 71.844

def layers_from_bulk(fe_mg, si_mg, ca_mg=0.07, al_mg=0.09, w_feo=0.0, w_o_core=0.0, w_s_core=0.0,
                     si_into_core=True):
    """Per mole of Mg.  w_feo = FeO mass fraction of the mantle.  If si_into_core, the mantle FeO is made
    by 3 Fe + SiO2 -> 2 FeO + FeSi ([R4467] reaction 17), so core Si moles = FeO moles / 2.
    w_o_core, w_s_core = O and S mass fractions of the core (extra light elements, not tied to the mantle)."""
    a0 = MU_MGO + ca_mg * MU_CAO + al_mg * MU_ALO15
    k = a0 + si_mg * MU_SIO2
    if si_into_core:
        e = w_feo * k / (MU_FEO - w_feo * (MU_FEO - 0.5 * MU_SIO2))
        s_core = 0.5 * e
    else:
        e = w_feo * k / (MU_FEO * (1.0 - w_feo))
        s_core = 0.0
    if e > fe_mg:
        return None                               # not enough iron to make that much FeO
    m_mantle = k - s_core * MU_SIO2 + e * MU_FEO
    m_fe, m_si = (fe_mg - e) * AMU["Fe"], s_core * AMU["Si"]
    m_core = (m_fe + m_si) / (1.0 - w_o_core - w_s_core)
    cmf = m_core / (m_core + m_mantle)
    w = {"Fe": m_fe / m_core, "Si": m_si / m_core, "O": w_o_core, "S": w_s_core}
    mu_core = 1.0 / sum(w[el] / AMU[el] for el in w)           # ideal-mixing mean atomic mass
    return {
        "cmf": cmf, "x_fe": e / (1.0 + e),                      # Fe/(Fe+Mg) in the mantle
        "core_wt": w, "core_scale_le": mu_core / AMU["Fe"],     # density factor from light elements
        "mantle_wt": {"MgO": MU_MGO / m_mantle, "SiO2": (si_mg - s_core) * MU_SIO2 / m_mantle,
                      "FeO": e * MU_FEO / m_mantle, "CaO": ca_mg * MU_CAO / m_mantle,
                      "Al2O3": 0.5 * al_mg * (2 * AMU["Al"] + 3 * AMU["O"]) / m_mantle},
        "mantle_si_mg": si_mg - s_core, "fe_total_wt": fe_mg * AMU["Fe"] / (m_core + m_mantle),
    }

def thermal_deficit(p_gpa, t_core_k):
    """Fractional density loss of hcp-Fe at t_core_k relative to 300 K at p_gpa (fit to [R4233], 135-3000 GPa,
    2500-6000 K, about +/-15 percent of the deficit)."""
    return 0.0525 * ((t_core_k - 300.0) / 3700.0) * (135.0 / p_gpa) ** 0.84

def planet_from_composition(M_earth, fe_mg, si_mg, ca_mg=0.07, al_mg=0.09, w_feo=0.0, w_o_core=0.0,
                            w_s_core=0.0, t_core_k=4000.0):
    """Route 2: composition -> structure.  Two passes: the thermal deficit needs P_cmb."""
    lay = layers_from_bulk(fe_mg, si_mg, ca_mg, al_mg, w_feo, w_o_core, w_s_core)
    if lay is None:
        return None
    scale = lay["core_scale_le"]
    for _ in range(2):
        st = structure(M_earth, lay["cmf"], lay["x_fe"], scale)
        p = (st["Pcmb"] or 3e11) / 1e9
        scale = lay["core_scale_le"] * (1.0 - thermal_deficit(p, t_core_k))
    st["layers"] = lay
    st["core_scale"] = scale
    return st

def route1(M_earth, R_earth, u, n=200, feo_wt_range=(0.0, 0.20), core_scale_range=(0.80, 1.0)):
    """Route 1: drawn (M, R) -> CMF with uncertainty from nuisance draws.  u() returns a uniform in [0,1)
    from the seeded wrapper.  Returns None if the planet is too large for a dry rock (hand to the water/gas lane)."""
    cm, fe = [], []
    for _ in range(n):
        w = feo_wt_range[0] + u() * (feo_wt_range[1] - feo_wt_range[0])
        x = x_from_feo_wt(w)
        s = core_scale_range[0] + u() * (core_scale_range[1] - core_scale_range[0])
        c = cmf_from_mr(M_earth, R_earth, x, s)
        if c is not None:
            cm.append(c)
            fe.append(c + (1.0 - c) * x * AMU["Fe"] / (x * (AMU["Fe"] + AMU["Si"] + 3 * AMU["O"]) +
                                                      (1 - x) * (AMU["Mg"] + AMU["Si"] + 3 * AMU["O"])))
    if len(cm) < n // 4:
        return None
    cm.sort(); fe.sort()
    q = lambda a, p: a[min(len(a) - 1, int(p * len(a)))]
    return {"cmf_median": q(cm, 0.5), "cmf_16_84": (q(cm, 0.16), q(cm, 0.84)), "cmf_min_max": (cm[0], cm[-1]),
            "fe_total_median": q(fe, 0.5), "fe_total_16_84": (q(fe, 0.16), q(fe, 0.84)), "feasible": len(cm) / n}

def x_from_feo_wt(w, si_per_cation=1.0):
    """Mantle Fe/(Fe+Mg) from FeO mass fraction for a mantle with Si/(Fe+Mg) = si_per_cation (bisection)."""
    lo, hi = 0.0, 0.99
    for _ in range(60):
        x = 0.5 * (lo + hi)
        f = MU_FEO * x / (MU_FEO * x + MU_MGO * (1 - x) + si_per_cation * MU_SIO2)
        if f < w:
            lo = x
        else:
            hi = x
    return 0.5 * (lo + hi)
```

Run time in CPython (COMPUTED): about 0.09 s per `structure` call, 4 microseconds per `radius_fit`, 0.04 s for `route1` with 200 draws; call `structure` once per body at generation and use `radius_fit` for anything repeated.

### 3.9 Check numbers (output of the code above, COMPUTED)

```text
1. [R4467] Eq.8: Fe/Mg 0.71, Si/Mg 0.79 (Hypatia means), no FeO, no light elements
   CMF 0.2917 (paper: mean 0.29, mu_bar 97.2)
2. Earth bulk (Fe/Mg 0.9, Si/Mg 0.9, Ca/Mg 0.07, Al/Mg 0.09), no FeO: CMF 0.3282
   + 8 wt% mantle FeO via reaction 17: CMF 0.2954, core Si 3.73 wt%, mantle Si/Mg 0.840, x_Fe 0.1073, core scale (Si only) 0.9644
   mantle oxides wt%: {'MgO': 37.35, 'SiO2': 46.76, 'FeO': 8.0, 'CaO': 3.64, 'Al2O3': 4.25} sum 1.0000
3. Solver, M=1 CMF 0.3 x=0.11 scale 0.9:  {'R': 0.9996, 'Pc': 370491451548.3137, 'Pcmb': 145072429958.4354, 'Rcore': 3329388.4363, 'g': 9.8278, 'v_esc': 11188.2896, 'moi': 0.331, 'rho_bar': 5519.729}
   fit R: 0.9997
4. Zeng&Sasselov-type pure silicate (CMF 0, x 0): R(1 M_E) 1.0774 R(10) 2.0039
5. Route 2, Earth bulk, 8 wt% FeO: R 0.9993 CMF 0.2954 core_scale 0.9176 g 9.835 m/s2 vesc 11.190 km/s MoI 0.3302 Pcmb 148.0 GPa Pc 379 GPa
6. Route 1, M=1, R=1.0: {'cmf_median': 0.295, 'cmf_16_84': (0.242, 0.343), 'cmf_min_max': (0.194, 0.379), 'fe_total_median': 0.34, 'fe_total_16_84': (0.315, 0.366), 'feasible': 1.0}
   Route 1, M=5, R=1.55: {'cmf_median': 0.286, 'cmf_16_84': (0.231, 0.339), 'cmf_min_max': (0.188, 0.393), 'fe_total_median': 0.345, 'fe_total_16_84': (0.315, 0.367), 'feasible': 1.0}
   Route 1, M=5, R=2.3 (too big for dry rock): None
7. Fit inversion round trip: CMF from (M=3, R=fit(3,0.37,.15,.85)) = 0.37
```

Reading the check numbers: item 1 reproduces [R4467]'s mean CMF (0.29, mu_bar 97.2) with the Hypatia means (0.2917). Item 2 with 8 wt% mantle FeO gives core Si 3.73 wt% (paper 3.7), mantle Si/Mg 0.840 (paper text 0.83). Items 3 and 5 reproduce an Earth-sized planet within 0.1 percent in radius with MoI 0.3302 and a CMB pressure 148 GPa (see 4 for why this is 9 percent above PREM). Item 6: for a planet of 1 Mearth and 1.0 Rearth, CMF 0.295 (0.242 to 0.343 at 16 to 84 percent, 0.194 to 0.379 extremes) and total Fe 0.34 (0.315 to 0.366).

### 3.10 Tests and invariants a test can assert

1. CMF identity: CMF from `layers_from_bulk` equals core mass over total mass (recompute from the returned wt%).
2. Mass-fraction closure: mantle oxide mass fractions sum to 1 within 1e-9; core wt% sum to 1.
3. Monotonic: R decreases with CMF at fixed M, increases with M at fixed composition, decreases with x_Fe.
4. Round trip: `cmf_from_mr(M, radius_fit(M, c, x, s), x, s)` returns c within 1e-3 (check: 0.37 returned for 0.37).
5. Fit vs solver: `radius_fit` within 0.5 percent of `structure` on a fixed grid.
6. Limits: coreless x=0 R(1 Mearth) in [1.07, 1.085]; pure Fe R(1 Mearth) in [0.76, 0.78].
7. Route 1 returns None for R above the coreless limit (M 5, R 2.3).
8. Seed determinism: the same `u()` stream gives the same composition.

## 4. Calibration and test targets

These are tests only; no default uses them.

| Target | Data (NASA fact sheets, via the scratchpad HTML) | Outcome | Label |
|---|---|---|---|
| Earth | M 5.9722e24 kg, R 6371 km, MoI 0.3308, core radius 3485 km | Route 2 from Earth-like bulk (Fe/Mg 0.9, Si/Mg 0.9, values as quoted by Unterborn; calibration only) with 8 wt% FeO: R = 0.9993, MoI 0.3302 (NASA 0.3308), g 9.835, v_esc 11.190 km/s, P_cmb 148.0 GPa (PREM about 135.8 GPa per G425, +9 %), Pc 379 GPa. Fitting CMF 0.326: core scale 0.835, MoI 0.3335, Rcore/R 0.554 (0.547), P_cmb 132.9 | COMPUTED |
| Earth, (R, MoI) fit | | CMF 0.286, scale 0.926, Rcore/R 0.508 | COMPUTED |
| Venus | M 4.8673e24, R 6051.8 km, MoI 0.337 | CMF 0.283, scale 0.850 | COMPUTED |
| Mars | M 0.64169e24, R 3389.5 km, MoI 0.366 | CMF 0.181, scale 0.705, Rcore/R 0.476 (NASA 0.54): a pure-Fe core cannot match; needs an S-rich core (Hakim: 0.8 rho_Fe at 50 mol% S) and 17 wt% FeO | COMPUTED |
| Mercury | M 0.33010e24, R 2439.7 km, MoI 0.35 | CMF 0.770, scale 0.769, Rcore/R 0.847 | COMPUTED |
| Moon | M 0.07346e24, R 1737.4 km, MoI 0.394 | CMF 0.046, scale 0.55; residual 0.74 % in R | COMPUTED |
| x from FeO | 8 wt% -> 0.116; 17.3 -> 0.262; 3.5 -> 0.050 | | COMPUTED |

Findings to flag: the Mars calibration needs a core density factor near 0.7; the Moon radius matches only to 0.74 percent; the solver P_cmb for an Earth-like composition is 9 percent above PREM, meaning an Earth-like composition drawn from the star does not reproduce PREM's exact core size; this is acceptable inside the stated +/-0.1 CMF uncertainty.

## 5. Gaps: Boss must research, and to verify

### Boss must research (still unreachable or not obtained)

Ranked.

1. **Stellar Ca and Al, and the planet-versus-star scatter**: `A-stars.md` now supplies Mg, Si, O, C (stored `mg_fe`, `si_fe`, `o_fe`, `c_fe`) but its stored list has no Ca or Al, which `layers_from_bulk` uses (small effect, a few percent on CMF); the planet-minus-star scatter law is still open (only a width of about 0.1 CMF is parametrised). Gives: Ca/Mg, Al/Mg per star and a published scatter.
2. **[R0535] full tables and Hypatia bulk table with the Fe/Mg-Si/Mg correlation**: VizieR refused; the Hypatia API needs star names. Gives the joint distribution to draw from.
3. **Mantle FeO / oxygen fugacity rule versus star and water** ([R1248] full text; Schaefer and Fegley): only an abstract and Putirka's 2 to 12 wt% range are available. Gives a proper `w_feo` law.
4. **Crust mass fraction and crust composition law** (part E owns the types); B hands over mantle oxides only.
5. **Volatile (C, N, H) inventory of the mantle and core**: no source gives ranges per planet; needed by the atmosphere lane. Candidates: [R3961], Bergin et al.
6. **K2O and Na2O mantle contents**: only the Na2O ceiling (Putirka) was read; K2O is part E's radiogenic table.
7. **Primary EOS papers**: [R4232] and 2024 originals (OUP refused), [R1096] [R0111, R1095, R3923] diamond, Occelli diamond, Sotin/Noack/[R4335, R4481, R4583]/2016. Gives original tabulated EOS to replace the 2ND-HAND copy.
8. **Low-pressure Fe-alloy EOS (Fe-Ni, Fe-S, Fe-Si, Fe-O) and liquid-core density** with temperature; only Hakim's 300-K-type refits were read. The 0.8 to 1.0 `core_scale` range is the stand-in.
9. **IOP-hosted full texts**: a bot check blocked access; not bypassed.

### To verify (RECALLED or weakly sourced)

1. Earth CMF 0.325 (used via G425 s.2.3 from [R4882]) and PREM CMB pressure 135.8 GPa (via G425): only test-side numbers.
2. Thermal-deficit fit coefficient (0.0525, 0.84): my own fit to [R4233] output.
3. Ideal-mixing light-element density rule for O (20 percent too low versus Hakim at 23 wt% O).
4. The probability p_SM for the super-Mercury branch (parametrisation).
5. The Moon and Mars calibration pairs (only as test hints).

## 6. Proposed work items for the TODO tree

| Title | Scope | Depends on | Size |
|---|---|---|---|
| B.1 Interface: A's ratios to molar B inputs | Implement the mass-to-molar conversion of 3.3 step 1 with A's factors; decide whether Ca/Mg, Al/Mg come from A or from the ExoPlex ranges | A's `mg_fe`, `si_fe` stored | small |
| B.2 Table: mantle density table | Store `_LP/_RHO` and the Fe factor with source rows (2.1) in a data file | B.1 | small |
| B.3 Table: EOS materials | Data module of 2.1 rows with source and label | none | small |
| B.4 Draw: composition from star | Implement 3.3 steps 1 to 4 with seeded draws | B.1, B.3 | medium |
| B.5 Draw: super-Mercury and coreless branches | p_SM and CMF range as parameters | B.4 | small |
| B.6 Derive: `layers_from_bulk` and oxides | Stoichiometry and closure test | B.3 | small |
| B.7 Derive: structure solver | `structure`, `radius_fit`, `cmf_from_mr` | B.2 | medium |
| B.8 Derive: route 1 fill-in | `route1` with nuisance priors, hand-off to C if too large | B.7 | small |
| B.9 Storage: new fields (3.2) in DB schema and save/load | cmf, mantle_feo_wt, core_w_*, core_scale, mg_si, mantle_oxides | B.4 | medium |
| B.10 Hand-over to E, C, atmosphere | Expose oxides, CMF, g, v_esc, P_surf API | B.6, B.7 | small |
| B.11 Tests | The eight invariants (3.10) and the calibration table | B.7 | small |
| B.12 Docs | Section for `docs/design` with the tables and the labelled gaps | all | small |

## 7. Verification summary

Counts are of table rows whose label cell carries the label (the prose lines that quote a READ paper are not counted again).

| Label | Rows | What |
|---|---|---|
| READ | 33 | [R0032, R0535, R0541, R1752, R3545, R3546, R3626, R3978, R4467]/2026, [R3384, R3800, R4879], [R2801, R2803], abstracts of [R1248, R4470] |
| 2ND-HAND | 46 | 32 [R4233] rows (via burnman 2.1.0), SLB hcp-Fe, liquid Fe and ExoPlex constants (sdist), diamond and graphite (burnman HP), SiC and Mg/Si thresholds (via G425) |
| COMPUTED | 31 | mass-radius tables, degeneracy, limits, carbon radii, thermal factor, stoichiometry check numbers, calibration fits, solver validation |
| SNIPPET | 0 | the earlier snippet-level slope claims were replaced by the full-text reads of [R0541, R3546] |
| RECALLED | 2 table items | S 32.06 and C 12.011 atomic weights; plus the test-side Earth CMF 0.325 and PREM P_cmb 135.8 GPa, taken from G425 (so 2ND-HAND there) and listed under "to verify" |
| NOT REACHED | 1 table row, 9 gap items | Fe-Ni EOS; the gap list in section 5 |

Agreement and differences with G425: agree on Seager constants, the [R4882] PREM calibration, the MgO K0 spread (177 vs 161.4 GPa) and the Mg/Si thresholds; differ in that the "100 percent Fe" Zeng column is not pure Fe (G425 already notes this) and that I recommend the integrated two-layer solver over the closed forms (P_cmb within 30 %, Eq. 22 within 70 %, MoI within 12 %).

