# GEN.533 bar and bulge sources (research)

Status: research read-out, 2026-10-17, written by a read-only research worker for the GEN.533 build (Core 1). Its proposed draws are implemented in `galaxy/structure_draw.py` (`draw_bar`). The rule applied to every conflict is Boss's: verified beats second-hand, the newest of two verified wins, the widest range wins where both are in question.
Status words: READ = number is in the paper text or a table printed in the PDF. FIG-LABEL = number printed inside a
figure panel, extracted as text. FIG-ONLY = visible only as a plot (not used for numbers). ABSTRACT = abstract only.
COMPUTED = I computed it from a table in the paper or its VizieR copy (method named).
Galaxy 'stage' in this file: early = Sa-Sb, late = Sc and later (papers use other cuts; each row says which).

## 0. Corrections to the lead list (check before citing)
- arXiv:1706.04625 is NOT Erwin & Debattista; it is an unrelated maths paper. The B/P paper is Erwin & Debattista 2017, arXiv:1703.01602.
- Gadotti 2009 (arXiv:0810.1953) is a bulge/bar/disc structure paper, not a bar-fraction paper. Its only bar fraction: 42% of disc galaxies (mass > 1e10, face-on SDSS), conclusions section, READ.
- Masters et al. 2011 is arXiv:1003.0449 (MNRAS 411, 2026; arXiv v2 Oct 2010). Sheth 2008 is arXiv:0710.4552. Nair & Abraham 2010 are 1004.0684 (bar fraction) and 1001.2401 (catalogue). Aguerri 2009 is 0901.2346. Diaz-Garcia 2016 is 1509.06743. Erwin 2018 is 1711.04867. Gadotti 2011 is 1003.1719. Erwin 2005 is astro-ph/0508590 (title: How Large Are the Bars in Barred Galaxies).
- Laurikainen & Salo 2017 = arXiv:1609.01936 (Barlenses and X-shape features compared). Not found as a separate 2017 item.
- Dwek 1995 and Martinez-Valpuesta 2006 aside, all PDFs below were fetched; Dwek 1995 is not on arXiv and was not read (G2 values are as given in the task).

## A. Bar fraction (any bar unless 'strong' is said)
| Quantity | Value (scatter / counts) | Paper | Place | Status |
|---|---|---|---|---|
| Any bar, S4G spirals, all masses | f=0.563+-0.019 (Sample 1, D<25 Mpc); 0.618+-0.020 (logM>=8.5); 0.623+-0.019 (Sample 2m, logM>=9) | Erwin 2018, 1711.04867 | Table 2 | READ |
| Strong (SB) / weak (SAB) split, same samples | f_SB 0.234 to 0.290; f_SAB 0.329 to 0.351 | Erwin 2018 | Table 2 | READ |
| Bar fraction vs stellar mass (S4G spirals) | quadratic logistic f=1/(1+exp(-(a+b1*m+b2*m^2))), m=log10(M/Msun), a=-82.2+-22.1, b1=17.1+-4.62, b2=-0.88+-0.24; peak ~0.70 at m~9.7; binned max 0.76 at m~9.8; declines both sides | Erwin 2018 | Sect. 3 text, Eq. 6, Fig. 5 | READ (coefficients), fit range m=8 to 12 |
| Same fit evaluated by me | m=9.48 (3e9): 0.69; 10.0: 0.69; 10.5: 0.58; 10.8: 0.46; 11.0: 0.36; 11.3 (2e11): 0.21 | Erwin 2018 Eq. 6 | COMPUTED | high-mass end is the thin end of S4G data; coverage not stated in text |
| Bar fraction vs g-r colour and gas fraction | roughly constant, g-r 0.1 to 0.8; logMHI/M* -2.5 to 1 | Erwin 2018 | abstract, Sect. 3 | READ |
| Bar fraction vs Hubble stage, visual (S4G, i<65deg) | T [0,3) (Sa,Sab): 69.8+-3.2% (N=212); T [3,5) (Sb,Sbc): 55.8+-3.6% (190); T [5,7] (Sc,Scd): 75.1+-2.4% (329); T>7: 60.6+-2.2% (498); S0 (T<0): 49.1+-4.6% (116) | Diaz-Garcia 2016, 1509.06743 | Table 1 (visual column), Sect. 5.1 | READ; errors binomial |
| Same, pooled by me to the Galaxy stages | early (T 0 to <5, Sa to Sbc): 254/402 = 0.632+-0.024; late (T>=5): 549/827 = 0.664+-0.016; difference 0.032+-0.029 (1.1 sigma) | Diaz-Garcia 2016 Table 1 | COMPUTED | Sbc (T=4) sits in the early bin; table cannot split T=3 from 4 |
| Ellipse-fit / Fourier columns of same table | T[0,3) 61.8%; T[3,5) 48.9%/46.8%; T[5,7] 60.2%/50.8%; T>7 37.1%/25.5% | Diaz-Garcia 2016 | Table 1 | READ; lower for late types by method |
| Abstract headline | ~75% for Sab, local minimum ~40% at T=4, drops below M*~1e9.5-1e10 | Diaz-Garcia 2016 | abstract, Fig. 10, Fig. 19 | ABSTRACT/FIG (does not match Table 1 bins; Table 1 preferred) |
| Bar fraction vs mass (visual) | no obvious trend above ~1e9; drop for faintest; ellipse/Fourier methods rise to ~50% (60% peak at 10^10.5) | Diaz-Garcia 2016 | Sect. 6.1, Fig. 19 | READ text; curves FIG-ONLY |
| SDSS strong bars, counts by type (no inclination cut) | Sa 287/1322 (22%); Sab 243/893 (27%); Sb 649/1890 (34%); Sbc 368/1207 (30%); Sc 379/1562 (24%); Scd 188/541 (35%); Sd 72/173 (42%); Sdm 13/52; Sm 18/68; all 2612/14034 (26%) | Nair & Abraham 2010a, 1001.2401 | Table 3 | READ; lower limits (text: +5% if unsure ones added) |
| Same pooled by me | Sa-Sb 1179/4105 = 0.287+-0.007; Sc-Sm 670/2396 = 0.280+-0.009; (Sa-Sbc 0.291; Sbc-Sm 0.288) | Nair & Abraham 2010a Table 3 | COMPUTED | strong bars only; no stage difference |
| SDSS strong bars vs stellar mass | ~40% at logM~9, falls to ~24% at 10.2, rises and plateaus ~30% to logM~11.4; N barred: S0+Sa 804, Sb 908, Sc+Sd 600 of 2312 | Nair & Abraham 2010b, 1004.0684 | Sect. 2 text, Fig. 1 legend | READ (text values); per-type curves FIG-ONLY |
| SDSS GZ2 strong bars, overall | 29.4+-0.5% (4020/13665), 0.01<z<0.06, logM~9 to 11 | Masters 2011, 1003.0449 | Sect. 3 | READ |
| GZ2 vs colour | blue cloud (g-r<0.6) constant ~20%; red sequence rises, >half of reddest; max ~70% at low L and g-r 0.7 to 0.8; 26+-1% for Mr>-20 to 37+-2% for Mr<-22 | Masters 2011 | Sect. 3.1, 3.1.1, 4.1 | READ text; curves FIG-ONLY |
| GZ2 vs bulge share (fracdeV) | rises monotonically; half of bulge-dominated discs barred; late-type (low fracdeV) low except reddest | Masters 2011 | Sect. 3.2, Figs 6-7 | abstract + text; no numbers per bin (FIG-ONLY) |
| GZ2 vs Hubble type (literature quoted) | RC3 strong bars ~60% (S0/a) falling to 30% (Sc), rising again later; Virgo 45-50% early vs 22-36% late spirals | Masters 2011 | Sect. 4.1 | second-hand |
| Local any-bar (SB+SAB) and strong | f_bar ~0.65, f_SB ~0.3 (consensus); SDSS local sample fSB 0.27+-0.05 (ElPa) / 0.35+-0.05 (visual); fbar independent of SED type at low z | Sheth 2008, 0710.4552 | Sect. 1, 4.1, 4.3 | READ |
| Ellipse-fit SDSS, by morphology | all 45%; S0 29%, early spirals (S0/a-Sb) 55%, late (Sbc-Sm) 54%; Fourier 26%/29%/33%/17%; visual re-check 38% (22/52/48%) | Aguerri 2009, 0901.2346 | abstract, Sect. 5 | READ; bars in less concentrated, bluer galaxies |
| NIR 2MASS spirals | 59% (67% with candidates), N=151; 36% for bars >4 kpc | Menendez-Delmestre 2007, astro-ph/0611540 | abstract, Sect. 3.1 | READ |
| Massive SDSS face-on | bars in 42% of disc galaxies, M>1e10 | Gadotti 2009, 0810.1953 | conclusions | READ |

Conflicts (A):
1. Mass trend. Erwin 2018 (S4G, f falls above 10^10) vs Nair & Abraham 2010b and Masters 2011 (SDSS, f flat or rising above 10^10.2; Masters rising with colour). Erwin 2018 attributes the gap to SDSS resolution (bars under twice the PSF are missed). Rule: both verified, newest wins: Erwin 2018 (2018) over Masters 2011 and N&A 2010. Widest mass range also Erwin (logM 8 to 12).
2. Stage trend. Diaz-Garcia 2016 (early 63% vs late 66%, no real difference), Aguerri 2009 (55% vs 54%), N&A 2010a strong bars (28.7% vs 28.0%) agree on no stage difference. Masters 2011 / RC3 say early spirals more barred, with a dip at Sc to Sbc. Rule picks Diaz-Garcia 2016 (newest verified with explicit stage bins).
3. Colour. Erwin 2018 flat vs Masters 2011 steep rise. Rule picks Erwin 2018.
4. Overall level: 0.29 (Masters, strong), 0.45 (Aguerri), 0.56 to 0.62 (Erwin), 0.65 (Sheth consensus), 0.59 to 0.67 (Menendez-Delmestre). These count different bar strengths; Erwin 2018 Table 2 is the one with SB/SAB split and the widest mass range.

## B. Bar shape, bar size, vertical thickness
| Quantity | Value (scatter) | Paper | Place | Status |
|---|---|---|---|---|
| Bar max isophotal ellipticity e (deprojected), S4G barred | mean 0.52+-0.01 (N=654); by T: S0 0.40+-0.02 (47); T[0,3) 0.51+-0.01 (131); [3,5) 0.52+-0.02 (93); [5,7] 0.53+-0.01 (198); >7 0.54+-0.01 (185); +- is error of the mean. By bar family (four classes, symbols lost in text extraction) 0.59, 0.55, 0.46, 0.36 | Diaz-Garcia 2016 | Table 3 (third column), Fig. 15; abstract 'spirals e~0.5, independent of T' | READ; abstract says S0 e~0.2 but Table 3 gives 0.40 (internal mismatch, spirals unaffected) |
| Same, full distribution | N=504 with e: mean 0.534, SD 0.148, median 0.546, 5th 0.26, 16th 0.38, 84th 0.69, 95th 0.74, min 0.10, max 0.92; b/a=1-e: mean 0.466, SD 0.148; mass 3e9-2e11 subset (407): mean 0.533, SD 0.149 | Diaz-Garcia 2016 VizieR J/A+A/587/A160 table a3 (Ell) joined to table a1 (M*) | COMPUTED by me (python3 -I) | includes S0 bars; no T column in the tables |
| e independent of host mass | stated | Diaz-Garcia 2016 | Sect. 6 text, Fig. 20 | READ |
| Bar e from 2D bulge/bar/disc fits (SDSS, M>1e10, face-on, N=291) | median 0.63, SD 0.10 (error 0.02); boxiness c median 2.7, SD 0.3; text: distribution peaks at ~0.6, about 20% higher than ellipse-fit peaks | Gadotti 2011, 1003.1719 | Fig. 1 panel labels; Sect. 3.1 text | FIG-LABEL (median, SD) + READ (peak 0.6) |
| Deprojected bar e, early (S0-Sb) / late (Sc-Sd) | early SB 0.49+-0.13, SAB 0.36+-0.15; late SB 0.64+-0.18, SAB 0.40+-0.20 | Erwin 2005, astro-ph/0508590 | Sect. 4.3 (Bar Size and Bar Strength) | READ |
| Projected bar axis ratio, 2MASS | typical b/a ~0.5 (e~0.5), skewed to high e, sharp drop near e~0.75; e not correlated with Hubble type | Menendez-Delmestre 2007 | abstract, Sect. 3.2.2 | READ (no SD) |
| Detection floor | smallest bar e 0.16 (Aguerri); e>~0.2 required (M-D07) | Aguerri 2009 bar-strength subsection; M-D07 Sect. 3.2.2 | READ |
| Strong-bar ideal | 'A strong bar has an axial ratio of ~1/5'; ovals ~0.85 | Kormendy & Kennicutt 2004, astro-ph/0407343 | Sect. 3.2 | READ (rule of thumb, no sample) |
| Bar e in simulations | maximum e tracked, values only in Fig. 5 | Martinez-Valpuesta 2006, astro-ph/0507219 | Fig. 5 | FIG-ONLY |
| B/P in barred galaxies, by mass | fraction = 1/(1+exp(-(a+b*m))), a=-42.79, b=4.13 (50% at m=10.37); 79% of bars at m>=10.4 vs 12% lower; overall B/P 0.524+-0.054 of GoodPA (84) barred galaxies; bar fraction in parent 63.4% (+3.5-3.6) | Erwin & Debattista 2017, 1703.01602 | abstract, Sect. 4.1 Table 2, Sect. 2 | READ; my evaluation: m=9.48 0.03, 10.0 0.18, 10.7 0.80, 11.0 0.93, 11.3 0.98 (COMPUTED) |
| B/P in barred galaxies, by stage | 77% of S0-Sbc bars, 15% of Sc-Sd bars (trend is a mass side effect); ~80% flat for S0-Sbc | Erwin & Debattista 2017 | abstract, Sect. 4 (Fig. 7 text) | READ |
| B/P size (in plane) vs bar | Rbox/Lbar = 0.42+-0.09 (range 0.25 to 0.76, median 0.43); Rbox/a = 0.53+-0.12; earlier 0.43+-0.10 (24 galaxies, EB13 1301.0638); inverted 0.38+-0.06 | Erwin & Debattista 2017 | text and Fig. 12 (section number not checked) | READ |
| X-shape size vs bar | a/rbar 0.2 to 0.5, mean ~0.35; barlens axis ratio b/a 0.5 to 1.0 (sky peak ~0.75, disc plane ~0.85) | Laurikainen & Salo 2017, 1609.01936 | Sect. 4, Fig. 9, conclusions | READ |
| Thin bar: length / thickness (edge-on, NIR) | 14+-4 (so thickness/length 0.071, range 0.056 to 0.10); literature range 4 to 15 (Sellwood & Wilkinson 1993) | astro-ph/0008478 | abstract, Sect. 6.4 | READ (the 4 to 15 is second-hand) |
| Peanut/X height vs its length (edge-on) | zmax/R ~0.5 to 0.6 for X/P with R<~2.5 kpc (11 galaxies); height 0.1h to 0.5h of disc | Ciambur & Graham 2016, 1603.00019 | Sect. 4.3, Fig. 11 | READ (no SD) |
| Derived peanut height / bar length | 0.42 x 0.55 = 0.23 (range 0.14 to 0.30 from 0.25-0.76 x 0.5-0.6 ) | Erwin & Debattista 2017 x Ciambur & Graham 2016 | COMPUTED, two samples multiplied | STAND-IN basis, not a measured ratio |
| Simulated bar flatness c/a | ~0.1 at formation, grows to ~0.2 over a Hubble time | Martinez-Valpuesta, Shlosman | Sect. 3 and Sect. 4 text | READ (simulation, not data) |
| Edge-on B/P frequency | 45% of ~1350 RC3 edge-on discs; >40% for all types S0 to Sd | astro-ph/0006359 | abstract | READ; E&D17 says this is high (Sect. 5) |
| Bureau & Athanassoula 2005 (astro-ph/0403226) | N-body kinematic diagnostics; no axis-ratio table found | | grep of full text | not found |
| MILKY WAY check only: Wegg & Gerhard 2013, 1308.0593 | exponential scale lengths 0.70:0.44:0.18 kpc = axis ratios 10:6.3:2.6 (y/x 0.63, z/x 0.26); projected 1:2.1; bar angle 27+-2 deg; test model x0=0.68, y0=0.28, z0=0.26 kpc (y0/x0 0.41, z0/x0 0.38) at 26.5 deg | | abstract; Sect. 7; model-test section (Eq. 16 context) | READ |

Conflicts (B):
1. Bar e. Diaz-Garcia 2016 mean 0.53 (ellipse-fit isophote, 2016) vs Gadotti 2011 median 0.63 (2D fit, 2011) vs Erwin 2005 early SB 0.49, late SB 0.64 vs M-D07 ~0.5. Rule: newest verified = Diaz-Garcia 2016 (also widest range: e 0.10 to 0.92, N=504). Caution recorded: Gadotti says ellipse-fit e is ~20% lower than the bar itself; the model's y0/x0 is a bar-body ratio, so the true centre may sit nearer 0.6 (b/a 0.4). Milky Way (e=0.61) lies 0.5 sigma above the Diaz-Garcia centre and 0.3 sigma below Gadotti.
2. Milky Way in-plane ratio: Dwek G2 0.392 (task) vs Wegg & Gerhard 0.63 (2013, exponential law) vs 0.41 (their Cao-model test). Definitions differ; no law taken from MW.
3. Vertical ratio, thin bar: 0.071 vs Martinez-Valpuesta 0.1 to 0.2 (sim). Rule: data over simulation, so 0.071. B/P part: derived 0.23 vs MW 0.26 (W&G) vs 0.38 (Cao model) vs 0.272 (Dwek). No direct external-galaxy B/P axis ratio was found: need edge-on 3D fits (z/x of the B/P body vs bar length) for a sample.

## C. Bulge axis ratios without a bar (B/A in plane, C/A vertical; A longest)
| Quantity | Value (scatter) | Paper | Place | Status |
|---|---|---|---|---|
| CALIFA unbarred bulges, most probable B/A, C/A | N=41 (clean) from Table C.1: B/A mean 0.83, SD 0.21, median 0.91; C/A mean 0.54, SD 0.19, median 0.54, 16th-84th 0.31-0.74 (43 rows incl. 2 with C/A>1: C/A mean 0.57) | Costantin et al. 2018, 1710.05222 | Table C.1 | COMPUTED from printed table (python3 -I) |
| Same, mass 3e9-2e11 (logM 9.48-11.3) | N=36: B/A 0.85+-0.20; C/A 0.55+-0.20 | Costantin 2018 Table C.1 | COMPUTED | |
| Barred-galaxy bulges, same survey (for contrast) | N=40: B/A 0.85+-0.13, C/A 0.58+-0.22 | Costantin 2018 Table C.2 | COMPUTED | |
| Distribution peaks | <B/A>=0.85, <C/A>=0.55; per-bulge errors sigma_B/A=0.15, sigma_C/A=0.25; 66% oblate, 19% prolate, 15% triaxial (79 bulges: 41 unbarred, 38 barred); 75% of triaxial in barred; low-n and late-type bulges are heterogeneous, high-n and early oblate | Costantin 2018 | Sect. 5 (text next to Fig. 7), abstract | READ |
| Unbarred S0-Sb bulges, equatorial | <E>=0.16, <B/A>=0.85; ~80% triaxial; PDF falls for B/A>0.93; no bulge with B/A<0.5; no difference by S0 vs spiral, n, luminosity | Mendez-Abreu et al. 2008, 0710.5466 | abstract, Sect. 5.6, Fig. 11-12 | READ |
| Flattening F=C^2/R^2 (R^2=(A^2+B^2)/2), unbarred S0-Sb, 115 clean bulges | mean of <F>=0.88 and of F1/2 =0.71; 33/55/43% have B/A<0.9 (max/mean/median); triaxiality bimodal (65% oblate-ish); C/R = sqrt(F) = 0.84 to 0.94 so C/A ~0.78 to 0.87 if B/A=0.85 | Mendez-Abreu et al. 2010, 1006.5449 | Sect. 4 (flattening), Eq. 41, Sect. 5 | READ (C/A conversion is mine) |
| Older statistics | <B/A>=0.86, <C/A>=0.65 (32 S0-Sb); Fathi & Peletier 03 <B/A>=0.79 early, 0.71 late; Noordermeer 07 <C/A>=0.55 (oblate assumed); Mosenkov 10 median C/A 0.63 (edge-on) | quoted in Mendez-Abreu 2010 (Intro) and Costantin 2018 (Intro) | | second-hand |
| Fathi & Peletier 2003 (astro-ph/0305326) | late-type bulges more elongated, KS 98.4%; ellipticity histograms in 0.1 bins only | | Sect. 3, Figs 2-4 | text gives no mean; FIG-ONLY |
| Pseudobulge flatness | e_bulge/e_disc>0.9 in 36% of S0-Sb and 51% of Sbc-Sm (35 each); Moellenhoff & Heidt: 10% | Kormendy & Kennicutt 2004 quoting Fathi & Peletier | Sect. 4.3 | READ (qualitative: many pseudobulges as flat as discs) |

Conflicts (C):
1. C/A. Costantin 2018 0.54 (SD 0.19) vs Mendez-Abreu 2010 ~0.78 to 0.87 (via F) vs Bertola 0.65, Mosenkov 0.63, Noordermeer 0.55 (second-hand). Rule: verified beats second-hand; of the two verified, newest is Costantin 2018 (also lower bound 0.19 to 1.0 range is wider). Older value kept: MA10.
2. B/A. Costantin 0.83 to 0.85, MA08 0.85, Bertola 0.86 agree; FP03 0.79/0.71 (second-hand) differ by stage. No stage split in Costantin's tables (no T column), so no stage dependence is drawn.

## D. Bar angle
| Quantity | Statement | Paper | Place | Status |
|---|---|---|---|---|
| Orientation of bar in its disc | 'bars are expected to be randomly oriented with respect to the underlying disk', dPA in 0 to 90 deg; ~10% within 10 deg of disc PA on geometry alone | Menendez-Delmestre 2007 | Sect. 3.1 | READ (expectation, not a measurement) |
| Mock-observation setup | bar angle drawn uniformly in (0, 90 deg) from the line of nodes; mock PAs uniform 0 to 180 deg | Erwin 2018 Sect. 5 (simulations); Aguerri 2009 Sect. 4 | READ (assumption) |
| Empirical check | 84 of 118 barred galaxies have dPA<=60 deg (0.71); uniform expects 0.667 (60/90); +1.0 sigma | Erwin & Debattista 2017 Sect. 2 | COMPUTED | consistent, weak test; selection favours small dPA |
No paper in this set states a measured uniform distribution; every statement is an assumption or expectation. By rotational symmetry the bar angle against any fixed reference direction is uniform on [0,180) unless the host has an external alignment (not found). Need for a measured test: bar PA relative to a fixed sky axis for a face-on, large-sample catalogue (e.g. Diaz-Garcia/Herrera-Endoqui PA columns, not tested here). The Milky Way 27 deg is one draw.

## Proposed draws (m = log10 stellar mass in Msun, range 9.48 to 11.30)
1. Bar present: p = 1/(1+exp(-(-82.2+17.1*m-0.88*m^2))) for BOTH early and late (0.69, 0.69, 0.58, 0.46, 0.36, 0.21 at m = 9.48, 10.0, 10.5, 10.8, 11.0, 11.3). Basis Erwin 2018 Eq. 6 (VERIFIED, read, any bar incl. weak). No stage factor because Diaz-Garcia 2016 gives early 0.632+-0.024 vs late 0.664+-0.016 (1.1 sigma). Flags: (a) the quadratic is least safe above m~11 (coverage not stated); an alternative from the same S4G data, Diaz-Garcia 2016 visual, is flat ~0.6 to 0.65 above 1e9, so p=0.36 to 0.21 at the top end is the rule-picked but widest-uncertainty value; (b) 'strong bar only' would be 0.23 to 0.29 (VERIFIED, Erwin Table 2) with no stage split (N&A Table 3: 0.287 vs 0.280).
2. Bar in-plane axis ratio y0/x0: draw e from a normal, centre 0.53, SD 0.15, clamp e to 0.15 to 0.90, then y0/x0 = 1-e (so mean 0.47, clamp 0.10 to 0.85). Numbers VERIFIED (Diaz-Garcia 2016 Table 3 centre 0.52 to 0.54; SD COMPUTED from the VizieR table, N=504; clamp from table range 0.10 to 0.92 and Aguerri e_min 0.16). Mapping isophote e to the model's y0/x0 is a STAND-IN; Gadotti 2011 (centre 0.63, SD 0.10) suggests the bar body is ~0.1 more elongated, so a build that wants the body ratio would shift to 0.63/0.10. Same draw for early and late (Diaz-Garcia Table 3 differences 0.51 vs 0.53; M-D07 finds no T trend).
3. Bar vertical axis ratio z0/x0, two branches:
   - Boxy/peanut branch with probability P_BP = 1/(1+exp(-(-42.79+4.13*m))) (0.03 at 3e9, 0.51 at 2.3e10, 0.98 at 2e11; Erwin & Debattista 2017, VERIFIED read): normal, centre 0.23, SD 0.06, clamp 0.12 to 0.38. Centre is a STAND-IN derived from Rbox/Lbar 0.42 x Ciambur & Graham zmax/R 0.55; SD set to 25% of the centre (E&D17's Rbox/Lbar scatter is 21%); upper clamp 0.38 from the Cao-model MW test fit (check only). Milky Way 0.272 lies +0.8 sigma.
   - Thin-bar branch (probability 1-P_BP): normal, centre 0.071, SD 0.02, clamp 0.05 to 0.15. Centre VERIFIED (length/thickness 14+-4, edge-on, galaxies that show peanuts); SD from the quoted +-4; clamp upper from simulated c/a 0.1 to 0.2 (Martinez-Valpuesta 2006, a simulation).
   Enforce z0 <= y0 after drawing (design choice, STAND-IN).
4. Unbarred bulge axis ratios (y0/x0 = B/A, z0/x0 = C/A): B/A normal centre 0.83, SD 0.21, clamp 0.30 to 1.00 (Costantin 2018 Table C.1, N=41, COMPUTED from the printed table; MA08 centre 0.85 and MA08 floor B/A>=0.5 agree on the centre but not the floor; clamp uses the wider table range). C/A normal centre 0.54, SD 0.19, clamp 0.19 to 1.00, then C/A = min(C/A, B/A) (last step STAND-IN). VERIFIED numbers, STAND-IN truncation and ordering. Same for early and late (no stage split in any verified source; Mendez-Abreu 2010 would put C/A near 0.8, recorded as the older, conflicting value).
5. Bar angle: uniform on [0, 180) deg, independent of everything else. STAND-IN (symmetry plus assumed-random statements in M-D07, Erwin 2018; no measured test found). Milky Way 27 deg is only a check draw.

Not found, with what is needed: (i) a direct external-galaxy B/P vertical axis ratio relative to bar length (needs an edge-on sample with 3D or deprojected B/P height and bar length, e.g. S4G edge-on); (ii) bar fraction by stage AND mass jointly in numbers (only figures in N&A 2010b Fig. 1c and Diaz-Garcia 2016 Fig. 19); (iii) bulge share (B/T) bins as numbers (Masters Fig. 6-7, Sheth Fig. 5 are plots); (iv) a measured test of bar-angle uniformity; (v) Dwek 1995 G2 paper text (not on arXiv).
Scratch: PDFs and scripts in /tmp/claude-0/gen533-research/ (set1 PDFs, set2 VizieR tables); nothing left in /mnt/project-files except this file.
