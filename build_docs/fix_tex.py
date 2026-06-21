"""
Post-process pandoc LaTeX output:
1. Replace floating \begin{figure} environments with inline images + bold captions.
2. Fix image paths for Overleaf (strip output/ prefix from media paths).
3. Fix bottom margin so page numbers are not clipped.
4. Insert clearpage before each lab/appendix section.
"""

import re
import sys


def fix_tex(filepath):
    with open(filepath, "r") as f:
        content = f.read()

    # 1. Strip output/ prefix from media paths (for Overleaf)
    content = content.replace("output/media_B/", "media_B/")
    content = content.replace("output/media_A/", "media_A/")
    content = content.replace("output/media_C/", "media_C/")

    # 2. Fix bottom margin for page numbers
    content = content.replace(
        r"\usepackage[top=2cm, bottom=1cm, left=2cm, right=1cm]{geometry}",
        r"\usepackage[top=2cm, bottom=2cm, left=2cm, right=1cm]{geometry}",
    )

    # 3. Replace figure floats with inline images
    pattern = re.compile(
        r'\\begin\{figure\}\s*\n'
        r'\\centering\s*\n'
        r'(\\pandocbounded\{\\includegraphics\[.*?\]\{.*?\}\})\s*\n'
        r'\\caption\{(.*?)\}\s*\n'
        r'\\end\{figure\}',
        re.DOTALL
    )

    def replace_figure(match):
        image_line = match.group(1)
        caption = match.group(2)
        return (
            f"\\begin{{center}}\n"
            f"{image_line}\n\n"
            f"\\textbf{{{caption}}}\n"
            f"\\end{{center}}"
        )

    content = pattern.sub(replace_figure, content)

    # 4. Insert \clearpage before each lab/section, appendix, and Shared Infrastructure
    # Match both \subsection{\texorpdfstring{X.\d and \subsection{X.\d
    for prefix in ['A', 'B', 'C']:
        content = re.sub(
            rf'(?=\\subsection\{{(?:\\texorpdfstring\{{)?{prefix}\.\d+)',
            '\\\\clearpage\n',
            content,
        )

    # Appendices: \section{...Appendix...
    content = re.sub(
        r'(?=\\section\{(?:\\texorpdfstring\{)?Appendix)',
        '\\\\clearpage\n',
        content,
    )

    # Shared Infrastructure section
    content = re.sub(
        r'(?=\\subsection\{Shared Infrastructure)',
        '\\\\clearpage\n',
        content,
    )

    # 5. Shrink the report title so it fits on one line
    content = re.sub(
        r'\\section\{From the DFT to the SPWVD: Time-Frequency Analysis Applied to\s+Neonatal\s+EEG\}\\label\{[^}]*\}',
        r'\\begin{center}{\\Large\\bfseries From the DFT to the SPWVD:\\\\ Time-Frequency Analysis Applied to Neonatal EEG}\\end{center}',
        content,
    )

    # 6. Fix Table A.2 column widths (35/25/20/20 instead of 25/25/25/25)
    if "Library / Function" in content:
        col = r'>{\raggedright\arraybackslash}p{(\linewidth - 6\tabcolsep) * \real{0.2500}}'
        old_block = (col + '\n  ' + col + '\n  ' + col + '\n  ' + col)
        c1 = col.replace('0.2500', '0.3500')
        c2 = col
        c3 = col.replace('0.2500', '0.2000')
        c4 = col.replace('0.2500', '0.2000')
        new_block = (c1 + '\n  ' + c2 + '\n  ' + c3 + '\n  ' + c4)
        content = content.replace(old_block, new_block, 1)

    with open(filepath, "w") as f:
        f.write(content)

    remaining_figs = content.count(r"\begin{figure}")
    clearpage_count = content.count(r"\clearpage")
    print(f"Fixed: {filepath}")
    print(f"  Remaining \\begin{{figure}}: {remaining_figs}")
    print(f"  \\clearpage inserted: {clearpage_count}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fix_tex.py <file.tex>")
        sys.exit(1)
    fix_tex(sys.argv[1])
