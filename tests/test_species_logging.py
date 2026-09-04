"""Species fallback must leave enough context to diagnose failed calls."""

import pytest

from eca_pp import agent
from eca_pp.standardize.species import _llm_infer


@pytest.mark.parametrize("failure", [True, False])
def test_unresolved_call_logs_its_cause(monkeypatch, caplog, failure):
    monkeypatch.setattr(agent, "check_available", lambda: None)

    def ask(**kwargs):
        if failure:
            raise RuntimeError("test transport failure")
        return {"species": "human", "confidence": 0.2}, None, [], {}

    monkeypatch.setattr(agent, "ask_json", ask)
    assert _llm_infer(["ACTB"], {}) is None
    assert ("test transport failure" if failure else "confidence 0.200") in caplog.text
