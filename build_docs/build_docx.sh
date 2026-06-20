#!/bin/bash
# Build all volumes to docx.
# Margins: top=2cm, bottom=1cm, left=2cm, right=1cm
# Fonts: main=Roboto, code=Consolas (set via reference.docx template)
# Page numbers: bottom center (set via reference.docx template)
#
# Usage: bash build_docx.sh

PANDOC=~/miniforge3/bin/pandoc
OUTDIR=output
DOCS=docs

mkdir -p "$OUTDIR"

# Build each volume from its docs/ directory so image paths resolve
for vol in volume_A volume_B; do
    if [ -f "$DOCS/$vol.md" ]; then
        echo "Building $vol.docx ..."
        (cd "$DOCS" && $PANDOC "$vol.md" \
            -o "../$OUTDIR/$vol.docx" \
            --from=markdown+tex_math_dollars+pipe_tables \
            --mathml \
            --highlight-style=tango \
            --reference-doc="../template/reference.docx" \
            2>&1 | grep -v "^$")
        echo "  -> $OUTDIR/$vol.docx"
    fi
done

echo "Done."
