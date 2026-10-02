# tests/test_scatter_progress.py

"""
PERF.9 and PERF.4: the bright-star scatter's progress. The main bar
counts each layer's expected work (`brightStars.layer_weight`), credited
star by star as layers report (`scatter_layer`'s `on_progress`), so the
near-empty edge layers no longer count as much as the dense middle; and
while layers finish slowly a second bar shows the stars of the layers
being drawn (`generate._LayerTracker`), mirrored to the web job's
progress file as its `detail` and shown on the Generate page.
"""

import json

import pytest

import generate
from stellarObjects import brightStars, progressFile

from tests.test_bright_star_scatter import EDGE_PC, E_VALUE, EXTENTS, SHAPE, THRESHOLD


class _Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


class _FakeProgress:
    """Just the `Progress` calls `_LayerTracker` makes."""

    def __init__(self):
        self.tasks = {}
        self.next_id = 0
        self.main_task = self.detail_task = None

    def add_task(self, description, total=None, completed=0, **fields):
        self.next_id += 1
        self.tasks[self.next_id] = {"description": description, "total": total, "completed": completed, **fields}
        return self.next_id

    def update(self, task_id, **kwargs):
        self.tasks[task_id].update(kwargs)

    def remove_task(self, task_id):
        del self.tasks[task_id]


def _fractions():
    return brightStars.band_fractions(THRESHOLD)


def test_the_expected_count_follows_what_the_scatter_places():
    expected = sum(brightStars.layer_expected_stars(SHAPE, layer, outer, EDGE_PC, E_VALUE, _fractions())
                   for layer, outer in EXTENTS)
    placed = [len(list(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, THRESHOLD, seed)))
              for seed in range(6)]
    mean = sum(placed) / len(placed)
    assert mean > 20
    # A coarse sample of the same means the draw uses: close, not exact.
    assert expected == pytest.approx(mean, rel=0.5)


def test_dense_layers_weigh_more_than_edge_layers_and_rings_count_too():
    fractions = _fractions()
    middle, middle_stars = brightStars.layer_weight(SHAPE, 0, 8, EDGE_PC, E_VALUE, fractions)
    edge, edge_stars = brightStars.layer_weight(SHAPE, 1, 6, EDGE_PC, E_VALUE, fractions)
    assert middle > edge
    assert middle_stars > edge_stars
    assert middle == pytest.approx(middle_stars + brightStars.RING_WEIGHT_STARS * 9)
    empty, empty_stars = brightStars.layer_weight(SHAPE, 40, 3, EDGE_PC, E_VALUE, fractions)
    assert empty_stars == 0.0
    assert empty == brightStars.RING_WEIGHT_STARS * 4


def test_a_band_weighs_less_than_the_whole_scatter():
    whole = brightStars.layer_expected_stars(SHAPE, 0, 8, EDGE_PC, E_VALUE, brightStars.band_fractions(THRESHOLD))
    band = brightStars.layer_expected_stars(SHAPE, 0, 8, EDGE_PC, E_VALUE,
                                            brightStars.band_fractions(THRESHOLD, THRESHOLD * 4))
    assert 0 < band < whole


def test_a_layer_reports_its_stars_and_ends_at_its_count():
    reports = []
    rows = list(brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, 7,
                                          on_progress=lambda done, estimate: reports.append((done, estimate))))
    assert len(rows) > 10
    dones = [done for done, _ in reports]
    assert dones == sorted(dones)
    assert reports[-1] == (len(rows), len(rows))
    assert all(done <= estimate + 1e-9 for done, estimate in reports)


def test_reporting_progress_leaves_the_stars_unchanged():
    plain = list(brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, 7))
    reported = list(brightStars.scatter_layer(SHAPE, 0, 8, EDGE_PC, E_VALUE, THRESHOLD, 7,
                                              on_progress=lambda done, estimate: None))
    assert plain == reported


def test_the_main_bar_credits_layers_in_progress_and_never_goes_back():
    progress, clock = _FakeProgress(), _Clock()
    tracker = generate._LayerTracker(progress, "Bright stars", {0: 100.0, 1: 20.0, -1: 20.0}, clock=clock)
    main = progress.tasks[tracker.task]
    assert main["total"] == 140.0 and main["percent"] is True
    assert progress.main_task == tracker.task
    tracker.layer_progress(0, 25, 50)
    assert main["completed"] == pytest.approx(50.0)
    tracker.layer_progress(0, 30, 120)  # the estimate grew: no credit taken back
    assert main["completed"] == pytest.approx(50.0)
    tracker.layer_done(0)
    assert main["completed"] == pytest.approx(100.0)
    assert main["description"] == "Bright stars (1 of 3 layers)"
    tracker.layer_progress(7, 1, 1)  # not a layer of this run
    assert main["completed"] == pytest.approx(100.0)


def test_a_late_report_for_a_finished_layer_is_ignored():
    """PERF.23: a report that reaches the tracker after its layer's
    `layer_done` (the channel's thread lagging behind) used to put the
    layer back in flight, counting it twice and ending the bar at 101%."""
    progress, clock = _FakeProgress(), _Clock()
    tracker = generate._LayerTracker(progress, "Bright stars", {0: 100.0, 1: 20.0}, clock=clock)
    main = progress.tasks[tracker.task]
    tracker.layer_progress(0, 40, 50)
    tracker.layer_done(0)
    tracker.layer_progress(0, 50, 50)  # late
    assert main["completed"] == pytest.approx(100.0)
    tracker.layer_progress(1, 10, 10)
    tracker.layer_done(1)
    tracker.layer_progress(1, 10, 10)  # late
    tracker.layer_done(1)  # and counted once only
    assert main["completed"] == pytest.approx(main["total"])
    assert main["description"] == "Bright stars (2 of 2 layers)"
    assert tracker.in_flight == {}


def test_a_weighted_bar_never_shows_more_than_100_percent(tmp_path, monkeypatch, generate_page):
    """PERF.23: the terminal, the progress file and the Generate page all
    cap a share at 100%, whatever the count says."""
    path = tmp_path / "progress.json"
    monkeypatch.setenv(progressFile.ENV_VAR, str(path))
    progress = generate._generation_progress()
    main = progress.add_task("Bright stars (41 of 41 layers)", total=200.0, percent=True)
    progress.main_task = main
    progress.update(main, completed=202.6, total=200.0)  # a new total writes the file at once
    assert generate._CountColumn().render(progress.tasks[main]).plain == "100%"
    data = json.loads(path.read_text())
    assert data["completed"] == 200.0
    job = {"id": "abc", "created_at": 0, "finished": False, "progress": dict(data, completed=202.6)}
    assert generate_page._job_view(job)["progress_text"] == "Bright stars (41 of 41 layers): 100%"


def test_slow_layers_add_a_second_bar_that_goes_once_they_speed_up():
    progress, clock = _FakeProgress(), _Clock()
    weights = {layer: 10.0 for layer in range(-5, 6)}
    tracker = generate._LayerTracker(progress, "Bright stars", weights, clock=clock)
    tracker.layer_progress(0, 5, 50)
    assert tracker.detail is None and len(progress.tasks) == 1
    clock.now += 31
    tracker.layer_progress(0, 10, 50)
    tracker.layer_progress(1, 4, 40)
    detail = progress.tasks[tracker.detail]
    assert progress.detail_task == tracker.detail
    assert detail["completed"] == 14 and detail["total"] == 90
    assert detail["description"].strip() == "Layers 0, 1: stars"
    # Still slow by the stricter line once the bar shows.
    clock.now += 25
    tracker.layer_done(1)
    assert tracker.detail is not None
    assert progress.tasks[tracker.detail]["description"].strip() == "Layer 0: stars"
    # Layers come quickly again: the bar goes.
    for layer in (2, 3, 4, 5, -1, -2):
        clock.now += 2
        tracker.layer_progress(layer, 1, 2)
        tracker.layer_done(layer)
    tracker.layer_progress(0, 20, 50)
    assert tracker.detail is None and progress.detail_task is None
    assert list(progress.tasks) == [tracker.task]


def test_the_progress_file_carries_the_second_bar(tmp_path, monkeypatch):
    path = tmp_path / "progress.json"
    monkeypatch.setenv(progressFile.ENV_VAR, str(path))
    progress = generate._generation_progress()
    main = progress.add_task("Bright stars (0 of 3 layers)", total=200.0, percent=True)
    progress.main_task = main
    progress.update(main, completed=50.0)
    detail = progress.add_task("  Layer 0: stars", total=80, completed=10)
    progress.detail_task = detail
    progress.update(detail, completed=20, total=80)
    data = json.loads(path.read_text())
    assert data["description"] == "Bright stars (0 of 3 layers)"
    assert data["completed"] == 50.0 and data["percent"] is True
    assert data["detail"]["description"] == "Layer 0: stars"
    assert data["detail"]["completed"] == 20 and data["detail"]["total"] == 80
    progress.remove_task(detail)
    progress.detail_task = None
    data = json.loads(path.read_text())
    assert data["detail"] is None


def test_the_count_column_shows_a_share_for_a_weighted_bar():
    progress = generate._generation_progress()
    weighted = progress.add_task("Bright stars", total=200.0, completed=50.0, percent=True)
    counted = progress.add_task("Sectors", total=40, completed=12)
    column = generate._CountColumn()
    assert column.render(progress.tasks[weighted]).plain == "25%"
    assert column.render(progress.tasks[counted]).plain == "12/40"


@pytest.fixture
def generate_page(monkeypatch):
    import web.generate_page as generate_page

    monkeypatch.setattr(generate_page, "url_for", lambda *a, **k: "/job")
    return generate_page


def test_the_job_view_shows_the_share_and_the_second_bar(generate_page):
    job = {"id": "abc", "created_at": 0, "finished": False, "progress": {
        "description": "Bright stars (3 of 81 layers)", "completed": 30.0, "total": 120.0, "percent": True,
        "detail": {"description": "Layer 0: stars", "completed": 1200, "total": 4000, "eta_s": 75.0},
    }}
    view = generate_page._job_view(job)
    assert view["progress_text"] == "Bright stars (3 of 81 layers): 25%"
    assert view["progress_detail_text"] == "Layer 0: stars: 1,200 of 4,000, about 1 m 15 s left"


def test_the_job_view_counts_plain_bars_and_drops_the_detail_when_finished(generate_page):
    job = {"id": "abc", "created_at": 0, "finished": True, "progress": {
        "description": "Sectors", "completed": 12, "total": 40,
        "detail": {"description": "Layer 0: stars", "completed": 1, "total": 2, "eta_s": None},
    }}
    view = generate_page._job_view(job)
    assert view["progress_text"] == "Sectors: 12 of 40"
    assert view["progress_detail_text"] == ""
