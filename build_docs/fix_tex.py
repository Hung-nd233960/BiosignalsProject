"""
Post-process pandoc LaTeX output:
1. Replace floating \begin{figure} environments with inline images + bold captions.
2. Fix image paths for Overleaf (strip output/ prefix from media paths).
3. Fix bottom margin so page numbers are not clipped.
4. Insert \clearpage before each lab/appendix section.
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

    # 4. Insert \clearpage before each lab section and appendix
    # Labs: \subsection{...B.1 - Lab..., \subsection{...B.2 - Lab..., etc.
    content = re.sub(
        r'(?=\\subsection\{\\texorpdfstring\{B\.\d+ - Lab)',
        r'\\clearpage\n',
        content,
    )
    # Also before: B.7, B.8 (which use an em-dash not hyphen — check both)
    # Already covered by the pattern above since they match B.\d+

    # Appendices: \section{...Appendix...
    content = re.sub(
        r'(?=\\section\{\\texorpdfstring\{Appendix|\\section\{Appendix)',
        r'\\clearpage\n',
        content,
    )

    # Volume C sections: \subsection{...C.1, C.2, etc.
    content = re.sub(
        r'(?=\\subsection\{\\texorpdfstring\{C\.\d+)',
        r'\\clearpage\n',
        content,
    )

    # Volume A sections: \subsection{...A.1, A.2, etc.
    content = re.sub(
        r'(?=\\subsection\{\\texorpdfstring\{A\.\d+)',
        r'\\clearpage\n',
        content,
    )

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
