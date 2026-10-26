"""M5 benchmark prompts: every reference spec exists and validates."""

from pathlib import Path

import pytest
import yaml
from validate import validate

ROOT = Path(__file__).resolve().parents[1]
PROMPTS = yaml.safe_load((ROOT / "examples" / "benchmarks" / "prompts.yaml").read_text())
REFS = [(b["id"], r) for b in PROMPTS for r in b["reference"]]


def test_prompts_are_complete():
    ids = [b["id"] for b in PROMPTS]
    assert len(ids) == len(set(ids))
    assert "nitration" in ids
    for b in PROMPTS:
        assert {"id", "prompt", "preset", "reference"} <= set(b), b
        assert b["prompt"].strip().endswith("."), b["id"]


@pytest.mark.parametrize("bench,ref", REFS, ids=[f"{b}:{Path(r).stem}" for b, r in REFS])
def test_reference_validates(bench, ref):
    path = ROOT / ref
    assert path.exists(), ref
    assert validate(path) == []
