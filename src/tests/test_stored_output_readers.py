# tests/test_stored_output_readers.py

"""GEN.800: every magnetosphere and stellar activity output a generator stores has a reader that is not the sheet (the audit
reports/unused-stored-values-audit-2026-10-24.md named these)."""

import pathlib
import re

import pytest

SOURCE = pathlib.Path(__file__).resolve().parents[1] / "planetgen"
READERS = ("physics", "generation")

# stored name -> files (under planetgen/) that write it and must not count as a reader
WRITERS = {
    "xuv_flux_earth": {"physics/magnetism.py", "physics/activity.py"},
    "xuv_exposure_index": {"physics/magnetism.py", "physics/activity.py"},
    "flare_irradiation_index": {"physics/magnetism.py", "physics/activity.py"},
    "magnetopause_rp": {"physics/magnetism.py"},
    "dipole_class": {"physics/magnetism.py"},
    "rotation_period_hours": {"physics/spin.py"},
    "flare_alpha": {"physics/activity.py"},
    "xuv_saturated": {"physics/activity.py"},
    "ozone_loss_flag": set(),
}

NEVER_READ = {
    # Asked of TODO (GEN.800): a star's spin axis and tilt have no physical reader until the stellar obliquity is modelled.
    "star spin axis and axial tilt": "no model reads the star's obliquity yet",
}
"""Stored outputs allowed no reader, with why; the audit's remainder."""


def _readers(name, writers):
    """Files under physics and generation that read `name` as an attribute, a getattr string or a subject key, other than the file that sets it."""
    pattern = re.compile(rf"(?<![\w.])(?:\w+\.)?{name}\b(?!\s*=[^=])|[\"']{name}[\"']")
    found = []
    for folder in READERS:
        for path in sorted((SOURCE / folder).rglob("*.py")):
            relative = path.relative_to(SOURCE).as_posix()
            if relative in writers:
                continue
            code = [line for line in path.read_text(encoding="utf-8").splitlines() if not line.lstrip().startswith("#")]
            if any(pattern.search(line) for line in code):
                found.append(relative)
    return found


@pytest.mark.parametrize("name", sorted(WRITERS))
def test_every_stored_output_has_a_reader_besides_its_writer(name):
    assert _readers(name, WRITERS[name]), f"{name} is stored and nothing in physics or generation reads it"


def test_the_duplicate_exposure_index_is_gone_from_the_mixture_subject():
    text = (SOURCE / "physics" / "mixture_subject.py").read_text(encoding="utf-8")
    assert "xuv_exposure_index" in text and "activity.exposure_index(star, distance_au)" in text
    assert text.count("exposure_index(") == 1      # only as the fallback for a distance the body did not store


def test_the_allowed_unread_list_names_only_what_the_audit_left():
    assert set(NEVER_READ) == {"star spin axis and axial tilt"}


def test_the_runaway_history_follows_the_stars_own_saturation_age():
    import types
    from planetgen.physics import activity, constants, runaway_history
    mass_sol, age = 0.3, 3.0
    star = types.SimpleNamespace(mass=mass_sol * constants.SOLAR_MASS_TO_KG, age=age, xuv_saturated=True, log_lx_lbol=None)
    table_age = activity.saturation_age_gy(mass_sol)
    assert runaway_history._saturation_age_gy(star, mass_sol) == max(age, table_age)
    star.xuv_saturated, star.log_lx_lbol = False, -4.0
    assert runaway_history._saturation_age_gy(star, mass_sol) == pytest.approx(activity.effective_saturation_age_gy(mass_sol, age, -4.0))
    star.log_lx_lbol = None
    assert runaway_history._saturation_age_gy(star, mass_sol) == table_age
