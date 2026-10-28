"""Write the CDXML header templates presets/cdxml/<preset>.xml from the ChemDraw references
(acs-style.md §6, CDXML export).

Each template is copied verbatim from reference_renders/chemdraw/m1_set_x<k>.cdxml: the XML
declaration and DOCTYPE, the root <CDXML> attributes, <colortable>, <fonttable> and the <page>
attributes. Captions (phase 2) need Helvetica Bold (compound numbers) and Helvetica Oblique (*n*-BuLi),
which m1_set_x* don't list: their <font> entries are copied verbatim from FONT_SOURCE
(scheme_multi_x1.0), ids included, after the m1_set fonts. Left out are what changes per file, which cdxml_writer fills in: the root
CreationProgram (ours, not ChemDraw's), Name and BoundingBox, the page id, size and contents. The templates are generated data, so the
style numbers in them (BondLength, LineWidth, LabelSize, ...) trace to the references;
tests/test_cdxml.py checks them against the presets and the other references at each scale.

Also presets/cdxml/<preset>.nicknames.xml from abbrev_x<k> (phase 4, acs-style.md §13.3): per
chem.ABBREVIATIONS key, the inner <fragment> of ChemDraw's nickname node in its left-justified
form (bond leaving to the left), its NodeType (Fragment or Nickname) and IsNickname, moved so the
attachment atom (the inner node at the outer node's p) is at the origin; the
ExternalConnectionPoint then lies on -x. Ids and attributes are ChemDraw's (cdxml_writer renumbers
them); the right-justified form is the same fragment turned 180° (to 0.01 pt in abbrev_x*).

Usage: uv run python scripts/cdxml_headers.py
"""

from __future__ import annotations

import sys
from pathlib import Path

from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_refs as mr  # noqa: E402

from acs_chem_svg.cdxml_writer import HEADER_DIR  # noqa: E402

SOURCES = {"acs1996": "m1_set_x1.0", "acs1996_x0.8": "m1_set_x0.8", "acs1996_x0.6": "m1_set_x0.6"}
NICKNAME_SOURCES = {"acs1996": "abbrev_x1.0", "acs1996_x0.8": "abbrev_x0.8", "acs1996_x0.6": "abbrev_x0.6"}
PER_FILE_ROOT = ("CreationProgram", "Name", "BoundingBox")
FONT_SOURCE = "scheme_multi_x1.0"
CAPTION_FONTS = ("Helvetica Bold", "Helvetica Oblique")


def _parse(ref: str):
    parser = etree.XMLParser(load_dtd=False, no_network=True, remove_blank_text=True)
    return etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), parser)


def header(ref: str) -> bytes:
    doc = _parse(ref)
    src = doc.getroot()
    root = etree.Element("CDXML", {k: v for k, v in src.attrib.items() if k not in PER_FILE_ROOT})
    for tag in ("colortable", "fonttable"):
        el = src.find(tag)
        if el is not None:
            el = etree.fromstring(etree.tostring(el))
            root.append(el)
    table = root.find("fonttable")
    have = {f.get("name") for f in table}
    for f in _parse(FONT_SOURCE).getroot().find("fonttable"):
        if f.get("name") in CAPTION_FONTS and f.get("name") not in have:
            table.append(etree.fromstring(etree.tostring(f)))
    page = src.find("page")
    etree.SubElement(root, "page", {k: v for k, v in page.attrib.items() if k != "id"})
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", doctype=doc.docinfo.doctype,
                          pretty_print=True)


def _xy(n) -> tuple[float, float]:
    x, y = n.get("p").split()[:2]
    return float(x), float(y)


def nicknames(ref: str) -> bytes:
    """The <nicknames> template: one <nickname key NodeType [IsNickname]> per ABBREVIATIONS key,
    holding the inner <fragment> of its first left-justified node in `ref`, relative to the
    attachment atom."""
    from acs_chem_svg.chem import ABBREVIATIONS
    from acs_chem_svg.cdxml_writer import pts
    found = {}
    for n in _parse(ref).getroot().iter("n"):
        if n.get("NodeType") not in ("Fragment", "Nickname"):
            continue
        key = "".join(s.text or "" for s in n.find("t").iter("s"))
        ax, ay = _xy(n)
        inner = n.find("fragment")
        ecp = next(m for m in inner.findall("n") if m.get("NodeType") == "ExternalConnectionPoint")
        ex, ey = _xy(ecp)
        if key in found or not (ex < ax and abs(ey - ay) < 0.005):  # the left-justified form
            continue
        frag = etree.fromstring(etree.tostring(inner))
        frag.attrib.clear()
        for m in frag.findall("n"):
            x, y = _xy(m)
            m.set("p", pts(x - ax, y - ay))
        el = etree.Element("nickname", {"key": key, "NodeType": n.get("NodeType")})
        if n.find("t").get("IsNickname"):
            el.set("IsNickname", n.find("t").get("IsNickname"))
        el.append(frag)
        found[key] = el
    missing = [k for k in ABBREVIATIONS if k not in found]
    if missing:
        raise ValueError(f"{ref}: no left-justified nickname node for {', '.join(missing)}")
    root = etree.Element("nicknames", {"source": f"{ref}.cdxml"})
    root.extend(found[k] for k in ABBREVIATIONS)
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", pretty_print=True)


def main() -> int:
    HEADER_DIR.mkdir(parents=True, exist_ok=True)
    for preset, ref in SOURCES.items():
        path = HEADER_DIR / f"{preset}.xml"
        path.write_bytes(header(ref))
        print(f"wrote {path.relative_to(mr.REPO)} from {ref}.cdxml")
    for preset, ref in NICKNAME_SOURCES.items():
        path = HEADER_DIR / f"{preset}.nicknames.xml"
        path.write_bytes(nicknames(ref))
        print(f"wrote {path.relative_to(mr.REPO)} from {ref}.cdxml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
