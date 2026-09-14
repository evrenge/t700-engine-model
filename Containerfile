# The environment every number in this repository was produced under.
#
# Until 2026-09-13 there was no container definition, no lock file and no recorded
# interpreter version: the only description of the environment was the sentence "work
# inside the t700-dev distrobox" in CLAUDE.md, and that container existed on exactly one
# machine. A replication whose environment cannot be reproduced is a replication with an
# unrecorded input.
#
#   podman build -t t700-dev -f Containerfile .
#   distrobox create --name t700-dev --image localhost/t700-dev
#   distrobox enter t700-dev
#
# `poppler-utils` is the prerequisite a Python dependency list cannot express: every
# digitizer under tools/ shells out to `pdftoppm`, `pdfimages` or `pdfinfo` to raster the
# source scan. `src/t700/` needs none of it -- the model reads only the committed CSVs.
FROM registry.fedoraproject.org/fedora:44

RUN dnf install -y \
        python3 python3-pip \
        poppler-utils \
        git \
        ibm-plex-sans-fonts ibm-plex-mono-fonts ibm-plex-serif-fonts \
    && fc-cache -f \
    && dnf clean all

# The three IBM Plex families are what `validation/plotstyle.py` sets Matplotlib to, so
# the figures in `site/` are typeset in the same faces as the page around them. Without
# them Matplotlib falls back to DejaVu Sans silently and the plots still render -- the
# style module says so rather than failing, because a missing font is not a reason to
# stop a validation run.

COPY pyproject.toml /tmp/t700/pyproject.toml
COPY requirements-recorded.txt /tmp/t700/requirements-recorded.txt

# Pinned to what produced the published numbers. `requirements-recorded.txt` is the
# record; loosening it is a decision, not a maintenance step -- see
# docs/notes/environment.md for what was measured across versions.
RUN python3 -m pip install --no-cache-dir -r /tmp/t700/requirements-recorded.txt
