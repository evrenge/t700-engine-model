#!/usr/bin/env bash
# Rerun every digitizer and prove the committed data is what they produce today.
#
# A digitized CSV is only as trustworthy as its generator. If a tool has been improved
# since its output was written, the file is stale and every downstream number inherits it.
# This reruns all of them and leaves `git status` as the verdict: clean means current.
#
#   bash tools/reproduce_all.sh
#
# Needs the rendered page images under validation/out/digitize/, and takes a few minutes.
set -u
cd "$(dirname "$0")/.."

if [ -n "$(git status --porcelain data/)" ]; then
    echo "data/ is already dirty -- commit or stash first, or the verdict is meaningless."
    exit 2
fi

for t in digitize_a1 digitize_a2 digitize_a6810 digitize_a7 digitize_a9 \
         digitize_fig678 digitize_fig910; do
    printf '%-22s ' "$t"
    if timeout 1800 python3 "tools/$t.py" >"/tmp/reproduce_$t.log" 2>&1; then
        echo "ok"
    else
        echo "FAILED (see /tmp/reproduce_$t.log)"
    fi
done

echo
if [ -z "$(git status --porcelain data/)" ]; then
    echo "VERDICT: every data file reproduces byte for byte. The data is current."
else
    echo "VERDICT: the following files differ from what their tool produces today --"
    git status --porcelain data/
    echo
    echo "Inspect with 'git diff data/'. If the tool is right, commit the new data and say"
    echo "in the message what moved and why. Never hand-edit a generated file: an"
    echo "annotation added to the output does not survive the next rerun, which is how the"
    echo "Figure 9 defect notes were lost once already."
fi
