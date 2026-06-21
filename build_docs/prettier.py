"""
prettier.py — Touch up markdown files before pandoc conversion.
Run after writing reports, before building docx/pdf.

Usage: python prettier.py docs/volume_B.md
       python prettier.py docs/*.md
"""

import sys
import re


def prettify(text):
    """Apply all touch-ups to markdown text."""

    # 1. Remove horizontal rules (--- or *** on their own line)
    text = re.sub(r"^\s*[-*]{3,}\s*$\n?", "", text, flags=re.MULTILINE)

    # 2. Replace em dashes (—) with normal dashes (-)
    text = text.replace("—", "-")

    # 3. Replace en dashes (–) with normal dashes (-)
    text = text.replace("–", "-")

    # 4. Remove blockquote lines (> text) — internal editorial notes
    text = re.sub(r"^\s*>.*$\n?", "", text, flags=re.MULTILINE)

    # 5. Remove consecutive blank lines (keep at most one)
    text = re.sub(r"\n{3,}", "\n\n", text)

    return text


def main():
    if len(sys.argv) < 2:
        print("Usage: python prettier.py <file1.md> [file2.md ...]")
        sys.exit(1)

    for filepath in sys.argv[1:]:
        with open(filepath, "r") as f:
            original = f.read()

        result = prettify(original)

        with open(filepath, "w") as f:
            f.write(result)

        n_dashes = original.count("—") + original.count("–") - result.count("—") - result.count("–")
        n_rules = len(re.findall(r"^\s*[-*]{3,}\s*$", original, re.MULTILINE))
        n_quotes = len(re.findall(r"^\s*>", original, re.MULTILINE))
        print(f"{filepath}: removed {n_rules} horizontal rules, {n_quotes} blockquotes, replaced {n_dashes} fancy dashes")


if __name__ == "__main__":
    main()
