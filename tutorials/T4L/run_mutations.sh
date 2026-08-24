#!/bin/bash
# Set up every mutation listed in a file, both legs of the thermodynamic cycle.
#
#   bash run_mutations.sh mutations_neutral.txt [cluster]
#
# The water sphere is rebuilt for each mutation, centred on the residue being
# mutated. That is not a detail: a residue far from the centre gets neutralised by
# `qprep_prot` as an out-of-sphere charge, and a neutralised residue is no longer
# the residue the mutation names.
#
# Expects, in this directory:
#   2LZM_prep.pdb            the prepared protein
#   <MUT><POSITION>.pdb      one per mutation, from `pymol -c mutagenesis.pml`
set -euo pipefail

MUTATIONS=${1:-mutations_neutral.txt}
CLUSTER=${2:-SNELLIUS}
STRUCTURE=2LZM_prep.pdb
FF=OPLSAAM
RADIUS=25
CHAIN=A

here=$(pwd)
mkdir -p protein tripeptide

# `|| [ -n "$mutation" ]` so a final line without a trailing newline is still read.
while read -r mutation || [ -n "$mutation" ]; do
    mutation=$(echo "$mutation" | tr -d '[:space:]')
    [ -z "$mutation" ] && continue
    position=$(echo "$mutation" | grep -oE '[0-9]+')
    mutant=$(echo "$mutation" | sed -E 's/^[A-Z]+[0-9]+//')
    echo "=== $mutation ==="

    if [ ! -f "${mutant}${position}.pdb" ]; then
        echo "  skipping: ${mutant}${position}.pdb not found (run the PyMOL mutagenesis first)"
        continue
    fi

    work=$here/work/$mutation
    rm -rf "$work"; mkdir -p "$work"
    cp "$STRUCTURE" "$work/protein.pdb"
    cp "${mutant}${position}.pdb" "$work/"
    cd "$work"

    # Centre the sphere on the mutated side chain (CB; glycine has no CB, so HA3).
    cog=$(python - protein.pdb "$position" "$CHAIN" <<'PY'
import sys
pdb, position, chain = sys.argv[1], int(sys.argv[2]), sys.argv[3]
for line in open(pdb):
    if (line.startswith("ATOM") and int(line[22:26]) == position
            and line[21] == chain and line[12:16].strip() in ("CB", "HA3")):
        print(line[30:38].strip(), line[38:46].strip(), line[46:54].strip())
        break
PY
)
    if [ -z "$cog" ]; then
        echo "  skipping: no CB/HA3 found for residue $position of chain $CHAIN"
        cd "$here"; continue
    fi

    qprep_prot -i protein.pdb -FF "$FF" -r "$RADIUS" -cog $cog

    for leg in protein tripeptide; do
        qresfep -m "$mutation" -mc "$CHAIN" -S "$leg" -t A -FF "$FF" -c "$CLUSTER" \
                -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10 -clean dcd
        # Both legs are written as FEP_<mutation>; move each out before the next.
        mv "FEP_$mutation" "$here/$leg/"
    done

    cd "$here"
done < "$MUTATIONS"

echo
echo "Done. Submit with:"
echo "  for d in protein/FEP_* tripeptide/FEP_*; do (cd \$d && bash FEP_submit.sh); done"
echo "Then analyse with:"
echo "  qresfep_analyze -p protein -t tripeptide -T 298"
