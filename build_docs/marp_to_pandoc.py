"""
Convert Marp slide markdown to pandoc-compatible slide markdown.

Changes:
1. Remove Marp frontmatter (--- marp: true ... ---)
2. Add pandoc title block (% Title / % Author / % Date)
3. Convert ![w:900](path) to ![](path){ width=90% }
4. Remove HTML comments (<!-- ... -->)
5. Keep --- slide separators (pandoc uses them too with --slide-level=2)
"""

import re
import sys


def convert(inpath, outpath):
    with open(inpath, "r") as f:
        content = f.read()

    # 1. Remove Marp frontmatter
    content = re.sub(
        r'^---\nmarp:.*?^---\n',
        '',
        content,
        count=1,
        flags=re.MULTILINE | re.DOTALL,
    )

    # 2. Convert Marp image syntax: ![w:900](...) -> ![](...){ width=90% }
    content = re.sub(
        r'!\[w:\d+\]\(([^)]+)\)',
        r'![](\1){ width=90% }',
        content,
    )

    # 3. Remove HTML comments
    content = re.sub(r'<!--.*?-->', '', content, flags=re.DOTALL)

    # 4. Clean up multiple blank lines
    content = re.sub(r'\n{3,}', '\n\n', content)

    with open(outpath, "w") as f:
        f.write(content)

    slide_count = content.count('\n---\n') + 1
    print(f"Converted: {inpath} -> {outpath} (~{slide_count} slides)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python marp_to_pandoc.py input.md output.md")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
