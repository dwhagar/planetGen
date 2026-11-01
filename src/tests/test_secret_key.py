import multiprocessing
import os
import stat

from planetgen.util import secret_key


def _grab(path, queue):
    queue.put(secret_key.ensure("", path=path)[0])


def test_configured_key_wins_and_nothing_is_written(tmp_path):
    path = str(tmp_path / "secret_key")
    assert secret_key.ensure("explicit", path=path) == ("explicit", "configured")
    assert not os.path.exists(path)


def test_missing_key_is_made_once_mode_0600_and_reused(tmp_path):
    path = str(tmp_path / "state" / "secret_key")
    said = []
    key, kind = secret_key.ensure("", announce=said.append, path=path)
    assert kind == "created" and len(key) == 64
    assert stat.S_IMODE(os.stat(path).st_mode) == 0o600
    assert secret_key.ensure("", announce=said.append, path=path) == (key, "existing")
    assert len(said) == 1 and "missing" in said[0].lower() or "no secret_key" in said[0].lower()
    assert os.listdir(tmp_path / "state") == ["secret_key"]


def test_simultaneous_starts_end_with_one_key(tmp_path):
    path = str(tmp_path / "secret_key")
    queue = multiprocessing.Queue()
    procs = [multiprocessing.Process(target=_grab, args=(path, queue)) for _ in range(8)]
    for p in procs:
        p.start()
    keys = {queue.get(timeout=30) for _ in procs}
    for p in procs:
        p.join()
    assert len(keys) == 1 and keys == {secret_key._read(path)}


def test_unwritable_state_folder_falls_back_to_random(tmp_path):
    blocker = tmp_path / "file"
    blocker.write_text("x")
    key, kind = secret_key.ensure("", path=str(blocker / "secret_key"))
    assert kind == "ephemeral" and len(key) == 64


def test_the_key_never_appears_in_output(tmp_path, monkeypatch, capsys, caplog):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    monkeypatch.delenv("PLANETGEN_SECRET_KEY", raising=False)
    from planetgen.util import settings
    monkeypatch.setattr(settings, "get_settings", lambda *a, **k: type("S", (), {"secret_key": ""})())
    monkeypatch.setattr(secret_key, "get_settings", settings.get_settings)
    with caplog.at_level("DEBUG"):
        assert secret_key.main() == 0
    key = secret_key._read(str(tmp_path / "secret_key"))
    out = capsys.readouterr()
    assert key and key not in out.out + out.err and key not in caplog.text
    assert "setting one" in out.err
    assert secret_key.main() == 0
    again = capsys.readouterr()
    assert key not in again.out + again.err and again.err == ""
