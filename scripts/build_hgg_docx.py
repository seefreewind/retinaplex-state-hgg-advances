from pathlib import Path
import csv
import re

from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION_START
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Mm, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[1]
MD_PATH = ROOT / "submission" / "HGG_ADVANCES_MANUSCRIPT_v2.md"
OUT_PATH = ROOT / "submission" / "HGG_ADVANCES_MANUSCRIPT_v2.docx"
FIG_DIR = ROOT / "figures" / "final"


def set_cell_shading(cell, fill):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_borders(cell, color="D9D9D9", size="4"):
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    borders = tc_pr.first_child_found_in("w:tcBorders")
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tc_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = "w:" + edge
        element = borders.find(qn(tag))
        if element is None:
            element = OxmlElement(tag)
            borders.append(element)
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), size)
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), color)


def set_cell_margins(cell, top=90, start=100, bottom=90, end=100):
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for m, v in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn("w:" + m))
        if node is None:
            node = OxmlElement("w:" + m)
            tc_mar.append(node)
        node.set(qn("w:w"), str(v))
        node.set(qn("w:type"), "dxa")


def set_repeat_table_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def add_hyperlink(paragraph, text, url):
    part = paragraph.part
    r_id = part.relate_to(url, "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink", is_external=True)
    hyperlink = OxmlElement("w:hyperlink")
    hyperlink.set(qn("r:id"), r_id)
    run = OxmlElement("w:r")
    r_pr = OxmlElement("w:rPr")
    color = OxmlElement("w:color")
    color.set(qn("w:val"), "1F4E79")
    r_pr.append(color)
    underline = OxmlElement("w:u")
    underline.set(qn("w:val"), "single")
    r_pr.append(underline)
    run.append(r_pr)
    text_node = OxmlElement("w:t")
    text_node.text = text
    run.append(text_node)
    hyperlink.append(run)
    paragraph._p.append(hyperlink)


def add_markdown_text(paragraph, text):
    # Minimal inline Markdown support for readable manuscript typography.
    pattern = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`|https?://\S+)")
    pos = 0
    for match in pattern.finditer(text):
        if match.start() > pos:
            paragraph.add_run(text[pos:match.start()])
        token = match.group(0)
        if token.startswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.bold = True
        elif token.startswith("*"):
            run = paragraph.add_run(token[1:-1])
            run.italic = True
        elif token.startswith("`"):
            run = paragraph.add_run(token[1:-1])
            run.font.name = "Courier New"
        else:
            add_hyperlink(paragraph, token.rstrip(".,"), token.rstrip(".,"))
            trailing = token[len(token.rstrip(".,")):]
            if trailing:
                paragraph.add_run(trailing)
        pos = match.end()
    if pos < len(text):
        paragraph.add_run(text[pos:])


def add_body_paragraph(doc, text, style="Body Text"):
    p = doc.add_paragraph(style=style)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.line_spacing = 1.12
    add_markdown_text(p, text)
    return p


def add_bullet(doc, text):
    p = doc.add_paragraph(style="List Bullet")
    p.paragraph_format.space_after = Pt(3)
    p.paragraph_format.line_spacing = 1.05
    add_markdown_text(p, text)
    return p


def parse_pipe_row(line):
    return [x.strip() for x in line.strip().strip("|").split("|")]


def is_table_separator(line):
    """Recognize the Markdown delimiter row so it is not emitted as data."""
    cells = parse_pipe_row(line)
    return bool(cells) and all(re.fullmatch(r":?-{3,}:?", cell) for cell in cells)


def add_table(doc, rows, landscape=False):
    if not rows:
        return
    if landscape:
        sec = doc.add_section(WD_SECTION_START.NEW_PAGE)
        sec.orientation = WD_ORIENT.LANDSCAPE
        sec.page_width, sec.page_height = sec.page_height, sec.page_width
        sec.top_margin = Inches(0.55)
        sec.bottom_margin = Inches(0.55)
        sec.left_margin = Inches(0.55)
        sec.right_margin = Inches(0.55)
    table = doc.add_table(rows=1, cols=len(rows[0]))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    table.style = "Table Grid"
    hdr = table.rows[0]
    set_repeat_table_header(hdr)
    for i, val in enumerate(rows[0]):
        cell = hdr.cells[i]
        cell.text = val
        set_cell_shading(cell, "1F4E79")
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        set_cell_borders(cell)
        set_cell_margins(cell)
        for p in cell.paragraphs:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(7)
                run.font.name = "Arial"
    for row_idx, row_values in enumerate(rows[1:], start=1):
        cells = table.add_row().cells
        for i, val in enumerate(row_values):
            cell = cells[i]
            cell.text = val
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            if row_idx % 2 == 0:
                set_cell_shading(cell, "F2F6FA")
            set_cell_borders(cell)
            set_cell_margins(cell)
            for p in cell.paragraphs:
                p.paragraph_format.space_after = Pt(0)
                for run in p.runs:
                    run.font.size = Pt(6.5)
                    run.font.name = "Arial"
    doc.add_paragraph().paragraph_format.space_after = Pt(2)


def setup_document(doc):
    sec = doc.sections[0]
    sec.page_width = Inches(8.5)
    sec.page_height = Inches(11)
    sec.top_margin = Inches(0.85)
    sec.bottom_margin = Inches(0.85)
    sec.left_margin = Inches(0.9)
    sec.right_margin = Inches(0.9)
    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Arial"
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = RGBColor(0, 0, 0)
    normal.paragraph_format.line_spacing = 1.12
    normal.paragraph_format.space_after = Pt(6)
    body = styles["Body Text"]
    body.font.name = "Arial"
    body.font.size = Pt(10.5)
    body.font.color.rgb = RGBColor(0, 0, 0)
    for name, size, before, after in (("Heading 1", 15, 16, 8), ("Heading 2", 12, 12, 5), ("Heading 3", 11, 8, 4)):
        st = styles[name]
        st.font.name = "Arial"
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = RGBColor(0, 0, 0)
        st.paragraph_format.space_before = Pt(before)
        st.paragraph_format.space_after = Pt(after)
        st.paragraph_format.keep_with_next = True
    if "Caption Custom" not in [s.name for s in styles]:
        st = styles.add_style("Caption Custom", WD_STYLE_TYPE.PARAGRAPH)
    else:
        st = styles["Caption Custom"]
    st.font.name = "Arial"
    st.font.size = Pt(9)
    st.font.color.rgb = RGBColor(0, 0, 0)
    st.paragraph_format.space_after = Pt(6)


def add_figure(doc, number):
    path = FIG_DIR / f"Figure{number}.png"
    if not path.exists():
        return
    doc.add_picture(str(path), width=Inches(6.5))
    alt_text = {
        1: "Study architecture and diffuse genetic landscape across AMD, POAG, and RE",
        2: "Donor-aware broad retinal cellular effects across AMD, POAG, and RE",
        3: "Contrasting genome-wide genetic direction and cellular-context similarity",
        4: "Orthogonal gene-level and robustness evidence for selected broad-cell associations",
        5: "Boundary of cellular convergence across RNA, state, and open-chromatin evidence",
    }[number]
    inline = doc.inline_shapes[-1]._inline
    inline.docPr.set("descr", alt_text)
    inline.docPr.set("title", f"Figure {number}")
    p = doc.paragraphs[-1]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(8)


def main():
    doc = Document()
    setup_document(doc)
    lines = MD_PATH.read_text(encoding="utf-8").splitlines()
    current_heading = ""
    pending_table = []
    figure_for_heading = {
        "Ocular traits show heterogeneous and partly opposing genome-wide architecture": 1,
        "Donor-aware retinal mapping identifies a shared broad Müller-glial context": 2,
        "Cellular-context convergence can oppose genome-wide genetic direction": 3,
        "Orthogonal gene-level analysis supports selected broad-cell associations": 4,
        "Pooled retinal open-chromatin annotations do not provide independent class-level support": 5,
    }
    in_table_data = False
    table_landscape = False

    def flush_table():
        nonlocal pending_table, in_table_data
        if pending_table:
            add_table(doc, pending_table, landscape=table_landscape)
            pending_table = []
        in_table_data = False

    def flush_figure():
        if current_heading in figure_for_heading:
            add_figure(doc, figure_for_heading[current_heading])

    for raw in lines:
        line = raw.rstrip()
        if line.startswith("|"):
            if is_table_separator(line):
                continue
            if not pending_table and line.count("|") > 2:
                pending_table.append(parse_pipe_row(line))
            elif pending_table:
                pending_table.append(parse_pipe_row(line))
            continue
        if pending_table:
            flush_table()
        if line.startswith("# "):
            # Use a plain paragraph here: Word's built-in Title style can carry
            # a decorative blue bottom border in the user's default template.
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_after = Pt(18)
            run = p.add_run(line[2:].strip())
            run.font.name = "Arial"
            run.font.size = Pt(16)
            run.font.bold = True
            run.font.color.rgb = RGBColor(0, 0, 0)
            continue
        if line.startswith("## "):
            flush_figure()
            current_heading = line[3:].strip()
            doc.add_paragraph(current_heading, style="Heading 1")
            continue
        if line.startswith("### "):
            current_heading = line[4:].strip()
            doc.add_paragraph(current_heading, style="Heading 2")
            if current_heading.startswith("Table "):
                table_landscape = True
            continue
        if line.startswith("#### "):
            doc.add_paragraph(line[5:].strip(), style="Heading 3")
            continue
        if line.startswith("- "):
            add_bullet(doc, line[2:].strip())
            continue
        if not line.strip():
            continue
        if line.startswith("**") and line.endswith("**"):
            p = doc.add_paragraph()
            p.paragraph_format.space_after = Pt(4)
            run = p.add_run(line.strip("*"))
            run.bold = True
            run.font.name = "Arial"
            run.font.size = Pt(10.5)
            continue
        add_body_paragraph(doc, line)

    flush_table()
    flush_figure()

    # Restore portrait orientation for any trailing section created by tables.
    doc.core_properties.title = "Cross-trait genetic architecture reveals shared retinal cellular contexts despite opposing genetic direction"
    doc.core_properties.subject = "HGG Advances Article submission manuscript"
    doc.core_properties.author = "Da Lin; Ying Chen; Yue Liu; Yu Zhang"
    doc.save(OUT_PATH)
    print(OUT_PATH)


if __name__ == "__main__":
    main()
