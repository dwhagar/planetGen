# tests/test_datasets_hyg.py

"""VIEW.30: the HYG catalogue as a test oracle (`datasets/hyg_stars.py`): it is downloaded with the other datasets, read for its counts, and never required by bulk generation."""

import pytest

from planetgen.datasets import catalogue, fetch, hyg_stars

HEADER = '"id","hip","hd","hr","gl","bf","proper","ra","dec","dist","pmra","pmdec","rv","mag","absmag"\n'
ROWS = [
    '0,,,,"","",Sol,0.0,0.0,0.0,0.0,0.0,0.0,-26.7,4.85\n',                 # the Sun is not in anyone's sky
    '1,1,,,"","","",0.1,1.0,219.7,0.0,0.0,0.0,6.4,2.39\n',
    '2,2,,,"","","",0.2,1.0,100000.0,0.0,0.0,0.0,6.5,2.39\n',               # a distance of 100,000 pc is the catalogue's "unknown"
    '3,3,,,"","","",0.3,1.0,12.0,0.0,0.0,0.0,6.6,2.39\n',
    '4,4,,,"","","",0.4,1.0,,0.0,0.0,0.0,7.4,2.39\n',                       # no distance at all
    '5,5,,,"","","",0.5,1.0,30.0,0.0,0.0,0.0,,2.39\n',                      # no magnitude: left out
]


def test_the_reader_counts_the_stars_at_or_brighter_than_a_limit_and_leaves_out_the_sun(tmp_path, monkeypatch):
    path = tmp_path / hyg_stars.FILE
    path.write_text(HEADER + "".join(ROWS), encoding="utf-8")
    monkeypatch.setattr(hyg_stars, "stars", lambda: hyg_stars._magnitudes(str(path)))
    assert hyg_stars.count_at_or_brighter(6.5) == 2
    assert hyg_stars.count_at_or_brighter(6.5, known_distance=True) == 1
    assert hyg_stars.count_at_or_brighter(7.5) == 4
    assert hyg_stars.count_at_or_brighter(7.5, known_distance=True) == 2


def test_the_catalogue_is_fetched_by_setup_but_only_a_test_reads_it():
    assert hyg_stars.DATASET in catalogue.ALL and not hyg_stars.DATASET.required
    assert all(dataset.required for dataset in catalogue.ALL if dataset is not hyg_stars.DATASET)
    assert hyg_stars.DATASET.pieces[0].url.startswith("https://") and hyg_stars.DATASET.pieces[0].size == 33932548


def test_the_real_catalogue_holds_the_counts_the_plan_recalls():
    if fetch.locate(hyg_stars.DATASET).missing:
        pytest.skip("the HYG catalogue is not downloaded (python3 -m planetgen.cli.datasets fetch)")
    assert 8900 <= hyg_stars.count_at_or_brighter(6.5) <= 8940           # sky-view.md: about 8,920 at V <= 6.5
    assert 25700 <= hyg_stars.count_at_or_brighter(7.5) <= 25900


@pytest.mark.real_datasets
def test_bulk_generation_does_not_ask_for_the_catalogue(monkeypatch):
    """`catalogue.missing` (the gate in front of bulk generation) leaves out a dataset that is not required, however absent it is."""
    looked = []
    monkeypatch.setattr(fetch, "locate", lambda dataset: looked.append(dataset.key) or fetch.Located(dataset, "/here"))
    assert catalogue.missing() == []
    assert "hyg_stars" not in looked and "isochrones" in looked
