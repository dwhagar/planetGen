"""
The bulk composition accessor (GEN.428's first slice, `planetgen.physics.bulk_composition`).
"""

from types import SimpleNamespace

import pytest

from planetgen.physics import bulk_composition as bc


def _body(**overrides):
    values = dict(metal_mass_fraction=0.3, rock_mass_fraction=0.65, water_mass_fraction=0.05, ice_mass_fraction=0.0,
                  envelope_mass_fraction=0.0)
    values.update(overrides)
    return SimpleNamespace(**values)


def test_a_body_without_the_attributes_has_no_composition():
    assert bc.bulk_fractions(SimpleNamespace()) is None
    assert bc.bulk_fractions(_body(metal_mass_fraction=None)) is None
    # the hydrosphere's water share alone is not a composition
    assert bc.bulk_fractions(SimpleNamespace(water_mass_fraction=0.01)) is None


def test_the_accessor_reads_the_fractions_and_the_descriptors():
    body = _body(rock_feo_fraction=0.08, rock_al2o3_fraction=0.045, rock_mg_si_ratio=1.1, carbon_rock=True,
                 composition_origin="inverse")
    bulk = bc.bulk_fractions(body)
    assert bulk.shares() == {"metal": 0.3, "rock": 0.65, "water": 0.05, "ice": 0.0, "envelope": 0.0}
    assert (bulk.rock_feo, bulk.rock_al2o3, bulk.rock_mg_si) == (0.08, 0.045, 1.1)
    assert bulk.envelope_helium is None and bulk.porosity is None
    assert bulk.carbon_rock is True and bulk.origin == "inverse"
    assert bulk.dry_core_fraction == pytest.approx(0.3 / 0.95)


def test_fractions_must_add_up_to_one():
    with pytest.raises(bc.CompositionError, match="add up"):
        bc.bulk_fractions(_body(rock_mass_fraction=0.6))
    with pytest.raises(bc.CompositionError, match="negative"):
        bc.bulk_fractions(_body(rock_mass_fraction=0.7, water_mass_fraction=-0.05))
    with pytest.raises(bc.CompositionError, match="origin"):
        bc.bulk_fractions(_body(composition_origin="guess"))


def test_a_body_with_no_dry_part():
    bulk = bc.bulk_fractions(_body(metal_mass_fraction=0.0, rock_mass_fraction=0.0, water_mass_fraction=0.0,
                                   envelope_mass_fraction=1.0))
    assert bulk.dry_core_fraction is None


def test_to_attributes_round_trips():
    bulk = bc.BulkFractions(0.25, 0.7, 0.05, 0.0, 0.0, rock_feo=0.08, porosity=None, origin="drawn")
    attributes = bulk.to_attributes()
    assert set(attributes) == set(bc.BODY_FIELDS)
    assert bc.bulk_fractions(SimpleNamespace(**attributes)) == bulk
