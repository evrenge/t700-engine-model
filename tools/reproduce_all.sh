#!/usr/bin/env bash
# Rerun every digitizer and prove the committed data is what they produce today.
#
# A digitized CSV is only as trustworthy as its generator. If a tool has been improved
# since its output was written, the file is stale and every downstream number inherits it.
# This reruns all of them and leaves `git status` as the verdict: clean means current.
#
#   bash tools/reproduce_all.sh              # from inside the t700-dev container
#   distrobox enter t700-dev -- bash tools/reproduce_all.sh
#
# Needs poppler (pdftoppm, pdfimages) and takes a few minutes. It pulls the page rasters it
# needs from the PDF itself -- until 2026-09-13 digitize_a1 did not, and its one input was
# a gitignored intermediate that no tool in this list produced, so the gate could not be
# run from a fresh clone at all.
#
# ## Why this script checks its own interpreter before doing anything
#
# It used to run `python3` and compute its verdict from `git status data/` alone. Run from
# the Bazzite host, which has no NumPy, every digitizer crashed on import -- and `data/`
# was then clean *because nothing had written to it*, so a total failure to execute
# reported as "every data file reproduces byte for byte." A gate that cannot fail is not a
# gate, and CLAUDE.md lists this script as one of only two `live` mechanisms behind the
# rule that generated files are never hand-edited. Both holes are closed below: the
# interpreter is checked up front, and a failed tool now poisons the verdict.
set -u
cd "$(dirname "$0")/.."

if ! python3 -c 'import numpy' 2>/dev/null; then
    echo "this python3 has no NumPy, so every digitizer would crash on import and the"
    echo "verdict would be vacuous. Run inside the dev container:"
    echo
    echo "    distrobox enter t700-dev -- bash tools/reproduce_all.sh"
    exit 2
fi

# poppler is the prerequisite a Python dependency list cannot express: every digitizer
# shells out to pdftoppm or pdfimages to raster the source scan. Checked here for the same
# reason NumPy is -- without it the tools fail one by one with a FileNotFoundError naming
# a binary, which reads as a broken digitizer rather than a missing package.
for bin in pdftoppm pdfimages; do
    if ! command -v "$bin" >/dev/null 2>&1; then
        echo "$bin is not on PATH. The digitizers raster the source scan with poppler:"
        echo
        echo "    dnf install poppler-utils        # or the equivalent for your distro"
        echo
        echo "See docs/notes/environment.md."
        exit 2
    fi
done

if [ -n "$(git status --porcelain data/)" ]; then
    echo "data/ is already dirty -- commit or stash first, or the verdict is meaningless."
    exit 2
fi

failed=0
# Every tool that writes a committed file under data/. `digitize_multi_curve` was missing
# from this list until 2026-09-12, which left the seven Appendix C schedules under
# data/schedules/ with no reproducibility check at all; tests/test_reproducibility_gate.py
# now fails if a data-writing tool is absent here. `digitize.py` and `digitize_native.py`
# are shared libraries and write nothing, so they are correctly not listed.
# Entries are "<tool> [args]". `digitize_multi_curve` takes the figure key as an argument
# and only c23 has committed output -- c30 is open question #50 and the tool refuses it.
for entry in "digitize_a1" "digitize_a2" "digitize_single_curve" "digitize_a7" "digitize_a9" \
             "digitize_multi_curve c23" "digitize_c30" "digitize_fig678" \
             "digitize_fig910"; do
    set -- $entry
    t=$1
    shift
    printf '%-26s ' "$entry"
    if timeout 1800 python3 "tools/$t.py" "$@" >"/tmp/reproduce_$t.log" 2>&1; then
        echo "ok"
    else
        echo "FAILED (see /tmp/reproduce_$t.log)"
        failed=$((failed + 1))
    fi
done

echo
if [ "$failed" -ne 0 ]; then
    echo "VERDICT: $failed digitizer(s) did not run, so NOTHING is proven. A clean data/"
    echo "here means the tools never wrote to it, not that the files are current. Fix the"
    echo "failures above and rerun before trusting any number downstream of data/."
    exit 1
elif [ -z "$(git status --porcelain data/)" ]; then
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
