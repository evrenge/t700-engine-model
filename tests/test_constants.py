"""Table A.1 [TM-100991 pdf p.55], read from the page image, against `t700.constants`.

`TABLE_A1` is transcribed independently of the module, so a mistyped constant fails here
rather than moving both sides of a check.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from t700 import constants as c

TABLE_A1 = {
    "HVF": 18300.0,
    "J_GT": 0.0445,
    "J_PT": 0.062,
    "K_BL": 0.7826,
    "K_B3": 0.0025,
    "K_DAMP": 0.06854,
    "K_DPB": 0.03045,
    "K_H2": 0.239,
    "K_H3_1": 0.2496,
    "K_H3_2": -8.4,
    "K_H41_1": 0.3010,
    "K_H41_2": -86.905,
    "K_H45": 0.9623,
    "K_PS3": 0.956,
    "K_QC_1": 0.71,
    "K_QC_2": 0.29,
    "K_T41_1": 3.322,
    "K_T41_2": 288.7,
    "K_T45_1": 3.519,
    "K_T45_2": 179.1,
    "K_T49_1": 3.516,
    "K_T49_2": 172.3,
    "K_TH41_1": 0.0018326,
    "K_TH41_2": 0.0856,
    "K_TH45_1": 0.0018326,
    "K_TH45_2": 0.0856,
    "K_V3": 0.97,
    "K_V41": 6.17,
    "K_V45": 13.63,
    "K_WGT": 0.0876,
    "NG_DES": 44700.0,
    "NP_DES": 20900.0,
    "TC_T41": 0.29,
}
"""All 33 rows of Table A.1 in printed order [pdf p.55]. `_1`/`_2` are the printed
subscripts on subscripts (slope and intercept of one fit; see CLAUDE.md)."""

DERIVED = {"J_LOAD_UH60A"}
"""Module constants that are not Table A.1 rows; each cites its own derivation."""


@pytest.mark.parametrize("name", list(TABLE_A1))
def test_the_constant_is_the_printed_value(name):
    assert getattr(c, name) == TABLE_A1[name], (name, getattr(c, name), TABLE_A1[name])


def test_every_module_constant_is_checked_against_the_page():
    tree = ast.parse(Path(c.__file__).read_text())
    names = {
        node.target.id
        for node in tree.body
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
    }
    assert names - DERIVED == set(TABLE_A1), sorted(names ^ (set(TABLE_A1) | DERIVED))
