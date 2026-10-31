# tests/test_variability.py

"""GEN.944: the kind of variable a star is, from its state, class, temperature and mass (`physics/variability.py`)."""

from planetgen.physics import variability as v


def test_a_cepheid_is_a_classical_cepheid():
    assert v.kind_of("cepheid", "II", 6000.0, 5.0) == v.CEPHEID


def test_only_a_cool_asymptotic_giant_is_a_mira():
    assert v.kind_of("asymptotic-giant", "III", 3200.0, 1.5) == v.MIRA
    assert v.kind_of("asymptotic-giant", "III", 3800.0, 1.5) is None
    assert v.kind_of("red-giant", "III", 3200.0, 1.5) is None


def test_a_pre_main_sequence_star_is_a_t_tauri_or_a_herbig_star_by_mass():
    assert v.kind_of("pre-main-seq", "PMS", 4000.0, 0.8) == v.T_TAURI
    assert v.kind_of("pre-main-seq", "PMS", 9000.0, 2.0) == v.HERBIG
    assert v.kind_of("pre-main-seq", "PMS", 9000.0, 7.9) == v.HERBIG
    assert v.kind_of("pre-main-seq", "PMS", 9000.0, 8.0) is None
    assert v.kind_of("protostar", "protostar", 3000.0, 0.5) is None


def test_an_ordinary_star_has_no_kind():
    for state, letter in (("main-sequence", "V"), ("subgiant", "IV"), ("white-dwarf", "VII"), ("supergiant", "IB")):
        assert v.kind_of(state, letter, 5800.0, 1.0) is None
    assert v.kind_of(None, None, None, None) is None
