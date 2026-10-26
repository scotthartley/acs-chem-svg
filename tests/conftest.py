import sys
from pathlib import Path

# scripts/ holds measure_refs.py, compare.py and cdxml_to_molfile.py (dev tools the tests reuse)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
