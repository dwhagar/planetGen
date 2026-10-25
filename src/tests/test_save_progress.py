"""PERF.50: a worker reports the progress of one sector's save to the run's bar."""


from planetgen import tuning
from planetgen.generation import run_galaxy, steps
from planetgen.galaxy.spherical import sector_position_pc

from tests.test_phenomenon_scatter import EDGE_PC


class _Channel:
    def __init__(self):
        self.items = []

    def put(self, item):
        self.items.append(item)


def test_a_sector_save_reports_each_system_it_inserts(mysql_config):
    # No stored galaxy shape on purpose (TEST.375): with one, the sector is filled from the galaxy's density and the
    # 3 systems asked for do not limit it (hundreds of systems, over 15 minutes); the save's report is the same.
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 3
    address = (1, 0, 2)
    channel = _Channel()
    _sid, _name, sector = run_galaxy.generate_and_save_sector_at(
        args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC, channel=channel)
    kinds = [item[0] for item in channel.items]
    assert kinds[0] == "start" and kinds[-1] == "end"
    _tag, _key, name, stats_kind, total, _workers = channel.items[0]
    assert stats_kind == "save" and name.startswith("Saving sector")
    assert total == len(sector.entries) + len(sector.phenomena)
    assert channel.items[-1][2] is True


def test_a_save_over_the_threshold_draws_its_bar_in_the_parent(monkeypatch):
    """The relay turns a worker's save into a bar drawn at once when it is predicted long."""
    monkeypatch.setattr(tuning, "PROGRESS_BAR_SECONDS", 0.0)
    drawn = []

    class Display:
        def add_task(self, description, total=None, **kwargs):
            drawn.append(description)
            return 1

        def update(self, *_a, **_k):
            pass

        def remove_task(self, *_a):
            pass

    relay = steps.Relay(None, Display())
    relay.put(("start", "k", "Saving sector 1", "save", 10, 1))
    relay.put(("progress", "k", 4.0, None, None))
    relay.put(("end", "k", True))
    assert drawn == ["Saving sector 1"]
    assert not relay.steps
