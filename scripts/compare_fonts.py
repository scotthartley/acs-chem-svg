"""Compare fallback-font advance widths with system Helvetica (design.md §4.3 [VERIFY]).

Usage: uv run python scripts/compare_fonts.py
"""

from acs_chem_svg.fontmetrics import list_system_helvetica_faces, load_fallback, load_system_helvetica

CHARS = list("CHNOSPFlBrI0123456789+−()")
SIZE = 10.0


def main() -> None:
    print("System Helvetica.ttc faces:")
    for i, fam, sub in list_system_helvetica_faces():
        print(f"  #{i}: {fam} / {sub}")
    for sub in ("Regular", "Bold"):
        ref, alt = load_system_helvetica(sub), load_fallback(sub)
        print(f"\n{sub}: {ref.describe()} vs {alt.describe()}  (advance at {SIZE:g} pt)")
        worst = (0.0, "")
        for ch in CHARS:
            missing = [f.path.name for f in (ref, alt) if not f.has_char(ch)]
            if missing:
                print(f"  {ch!r}: MISSING in {missing}")
                continue
            a, b = ref.advance(ch, SIZE), alt.advance(ch, SIZE)
            d = b - a
            if abs(d) > abs(worst[0]):
                worst = (d, ch)
            print(f"  {ch!r} U+{ord(ch):04X}: {a:7.3f} {b:7.3f}  Δ {d:+.3f}")
        print(f"  max |Δ advance| = {abs(worst[0]):.3f} pt at {worst[1]!r}")
        # kerning pairs within the label set (resvg applies both fonts' kerning)
        diffs = []
        for a in CHARS:
            for b in CHARS:
                ka, kb = ref.kerning(a, b, SIZE), alt.kerning(a, b, SIZE)
                if ka or kb:
                    diffs.append((kb - ka, a + b, ka, kb))
        diffs.sort(key=lambda d: -abs(d[0]))
        for d, pair, ka, kb in diffs:
            print(f"  kern {pair!r}: Helvetica {ka:+.3f}  fallback {kb:+.3f}  Δ {d:+.3f}")
        print(f"  max |Δ kern| = {abs(diffs[0][0]) if diffs else 0:.3f} pt"
              + (f" at {diffs[0][1]!r}" if diffs else ""))


if __name__ == "__main__":
    main()
