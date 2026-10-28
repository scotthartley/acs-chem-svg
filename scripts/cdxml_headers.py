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


def main() -> int:
    HEADER_DIR.mkdir(parents=True, exist_ok=True)
    for preset, ref in SOURCES.items():
        path = HEADER_DIR / f"{preset}.xml"
        path.write_bytes(header(ref))
        print(f"wrote {path.relative_to(mr.REPO)} from {ref}.cdxml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
