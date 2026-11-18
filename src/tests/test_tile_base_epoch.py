"""MAP.431: the map's tile cache key (`query.galaxy_content_state`'s `tile_base`) follows what the tiles are drawn from: the epoch a stored shape or scatter was made under counts as the last epoch that changed
that stage (`generator_epoch.last_change`), so an epoch that left the stage alone does not make every tile again."""

import pytest

from planetgen.db import query
from planetgen.db import store as _db
from planetgen.util import generator_epoch
from tests.test_galaxy_gen import _seed_skeleton


@pytest.fixture
def history(monkeypatch):
    """Five epochs: 1 first, 2 plan and scatter, 3 fill only, 4 derived only, 5 scatter only."""
    entries = [{"epoch": number, "first_release": None, "reason": "t", "affects": list(stages)}
               for number, stages in enumerate([[], ["plan", "scatter", "fill"], ["fill"], ["derived"], ["scatter"]], start=1)]
    monkeypatch.setattr(generator_epoch, "epoch_history", lambda: tuple(entries))
    return entries


def test_the_last_change_of_a_stage_is_the_newest_epoch_up_to_the_recorded_one_that_altered_it(history):
    assert [generator_epoch.last_change("plan", epoch) for epoch in (1, 2, 3, 4, 5)] == [1, 2, 2, 2, 2]
    assert [generator_epoch.last_change("scatter", epoch) for epoch in (1, 2, 3, 4, 5)] == [1, 2, 2, 2, 5]
    assert [generator_epoch.last_change("derived", epoch) for epoch in (1, 2, 3, 4, 5)] == [1, 1, 1, 4, 4]
    assert generator_epoch.last_change("plan", None) is None
    assert generator_epoch.last_change("plan", 99) == 2, "an epoch this code does not know is read as the newest it knows"


@pytest.mark.fresh_database
def test_a_stamp_at_an_epoch_that_left_the_plan_and_the_scatter_alone_keeps_the_tile_base(mysql_config, history):
    _seed_skeleton(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute("UPDATE galaxy_shape SET generator_epoch = 2 WHERE id = 1")
        conn.execute("INSERT INTO scatter_epochs (kind, generator_epoch) VALUES ('mass', 2) ON DUPLICATE KEY UPDATE generator_epoch = 2")
        conn.commit()
        first = query.galaxy_content_state(conn)
        conn.execute("UPDATE galaxy_shape SET generator_epoch = 4 WHERE id = 1")
        conn.execute("UPDATE scatter_epochs SET generator_epoch = 4 WHERE kind = 'mass'")
        conn.commit()
        after_derived = query.galaxy_content_state(conn)
        conn.execute("UPDATE scatter_epochs SET generator_epoch = 5 WHERE kind = 'mass'")
        conn.commit()
        after_scatter = query.galaxy_content_state(conn)
    finally:
        conn.close()
    assert after_derived["tile_base"] == first["tile_base"], "epochs 3 and 4 changed neither the plan nor the scatter"
    assert after_scatter["tile_base"] != first["tile_base"], "epoch 5 changed the scatter: the scatter drawn under it is another one"
    assert after_derived["base"] != first["base"], "the release-keyed base still moves with the stored epoch"
