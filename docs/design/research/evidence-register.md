# The evidence register: whether the thing a stand-in stands for has been seen (GEN.1428)

Status: hand-kept, checked by `src/tests/test_stand_in_values.py` and shown as the evidence-status column of `stand-in-values.md`. The tier of a row says how good its source is; the status says whether the thing itself has been seen (Boss, 10:55 AM EST on 2026-11-08, from the unseen-compounds report, `reports/unseen-compounds-ranked-2026-11-08.md`).

There are four statuses: `observed in a planet`, `observed in a lab`, `predicted only` and `disputed` (the list is in `stand-in-values.md`). **A row is `predicted only` until a line here says otherwise**: observed in a planet or a lab only where the source says so. A line is `| row | status | note |`:

- `6.229:CV1` names one registry row (the section number, a colon, the fill ID) and sets its status.
- `6.196:*` sets the status of every row of that section at tier T1 or T2 (the rows whose value is a source's own); rows of the section at a lower tier stay `predicted only` unless a line names them.

The note says what the source saw, so the next reader can check it. A row's owner confirms its line, as for the tier. A later finding that promotes a predicted material (GEN.1429 keeps the watch lines) changes its line, not the row.

| row | status | note |
|---|---|---|
| 6.196:* | observed in a lab | melting curves and rheology constants of water ice are laboratory measurements |
| 6.129:* | observed in a lab | the strengths of the ice types are laboratory measurements |
| 6.228:SV1 | observed in a lab | the activity-model parameters are fits to laboratory vapor-pressure data |
| 6.228:SV2 | observed in a lab | the pure-acid vapor pressure is measured between 338 and 445 K; its cold extension is a calculation |
| 6.228:SV3 | observed in a lab | the water pressure over a solution is measured in the laboratory |
| 6.228:SV4 | observed in a lab | the supercooled water curve is a fit to laboratory data |
| 6.229:CV1 | observed in a planet | sulfuric acid droplets of about 75 percent acid make the Venus cloud (descent probes) |
| 6.229:CV2 | observed in a planet | sulfur dioxide is measured in the Venus air |
| 6.229:CV3 | observed in a planet | hydrogen chloride, hydrogen fluoride, carbonyl sulfide and carbon monoxide are measured in the Venus air |
| 6.229:CV4 | predicted only | sulfur monoxide is observed on Venus, disulfur monoxide is predicted only; the row holds both, so it takes the weaker status |
| 6.237:NM1 | observed in a lab | one group, a preprint of 2026-11-04 |
| 6.237:NM2 | observed in a lab | diamond-cell diffraction; the shock experiment disagrees on where it starts |
| 6.237:NM3 | observed in a lab | synthesised above 94 GPa, peroxide ion seen by Raman |
| 6.237:NM4 | observed in a lab | predicted, then synthesised |
| 6.237:NM5 | observed in a lab | made, survives release at room temperature |
| 6.237:NM6 | observed in a lab | made in a diamond cell |
| 6.237:NM7 | observed in a lab | two 2024 syntheses at ambient pressure |
| 6.238:UC1 | disputed | a solid-phase claim is made and doubted; a comment rebuts the 2020 measurement |
| 6.238:UC2 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC3 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC4 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC5 | observed in a lab | made and simulated in the laboratory, not shown on Titan |
| 6.238:UC6 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC7 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC8 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC9 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC10 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC11 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC12 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC13 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC14 | disputed | the detection is claimed and doubted |
| 6.238:UC15 | disputed | the detection is claimed and doubted |
| 6.238:UC16 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC17 | predicted only | predicted or computed only; nobody has seen it |
| 6.238:UC18 | predicted only | predicted or computed only; nobody has seen it |
| 6.287:AF1 | observed in a lab | the freezing points the curve is read from are laboratory measurements (Gable, Betz and Maron; Kunzler and Giauque) |
