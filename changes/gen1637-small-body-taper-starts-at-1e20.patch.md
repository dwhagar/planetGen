**Headline:** Small bodies of 1e19 to 1e21 kg keep more of their pores, so 300 to 450 km Kuiper belt bodies are less dense.

### Changed
- GEN.1637: the compaction taper (`physics/compaction.py`, a stand-in of decision S9) now closes the pores from 1e20 kg to 1e21 kg, not from 1e19 kg. Data audit D17 found a 70 percent rock body of 300 to 450 km at 1 657 to 1 836 kg/m3 where 2002 UX25, Varda, Salacia and Orcus measure 820 to 1 530 and keep their pores to about 5e20 kg. At a rock fraction of 0.6 Varda, Salacia and Orcus now come out within 15 percent (2002 UX25 stays 44 percent over). The compaction table is rebuilt (`scripts/build_compaction_table.py`); a stored belt member under 1e19 kg still takes its family's drawn porosity (`compaction.DRAWN_POROSITY_BELOW_KG`). The registry row MC1 moves to tweak level L2 (tuned to measured bodies).

generation-output: changed
generation-affects: fill
