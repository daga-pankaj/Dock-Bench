from backend import config


def test_find_obrms_uses_valid_environment_override(monkeypatch, tmp_path):
    executable = tmp_path / "obrms"
    executable.touch()
    monkeypatch.setenv("OBRMS_BIN", str(executable))
    monkeypatch.setattr(config, "_tool_candidates", lambda _tool: [])

    assert config.find_obrms() == str(executable)


def test_find_obrms_checks_bundled_tool_before_path(monkeypatch, tmp_path):
    executable = tmp_path / "obrms"
    executable.touch()
    monkeypatch.delenv("OBRMS_BIN", raising=False)
    monkeypatch.setattr(config, "_tool_candidates", lambda _tool: [executable])
    monkeypatch.setattr(config.shutil, "which", lambda _tool: None)

    assert config.find_obrms() == str(executable)


def test_find_obrms_falls_back_to_path(monkeypatch):
    monkeypatch.delenv("OBRMS_BIN", raising=False)
    monkeypatch.setattr(config, "_tool_candidates", lambda _tool: [])
    monkeypatch.setattr(config.shutil, "which", lambda tool: "/usr/bin/obrms" if tool == "obrms" else None)

    assert config.find_obrms() == "/usr/bin/obrms"
