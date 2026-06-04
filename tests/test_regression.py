"""Regression tests locking in LIDAR's variant-classification behavior.

These tests assert that the packaged feature-extraction + Random Forest
inference (`lidar.lidar.format_gatk` / `format_dv`) keeps exactly the
variants it is expected to. The golden values were verified to match the
original Jupyter-notebook implementation (Hybrid/Hybrid_Model.py) on the
same synthetic inputs, so any future change that silently alters the
feature column order, parsing, or model wiring will fail here.

Run with:  pytest -q
"""
from pathlib import Path

import pytest

from lidar.lidar import (
    format_dv,
    format_gatk,
    load_model,
    package_resource_path,
)

FIXTURES = Path(__file__).parent / "fixtures"

# Golden truth: POS values the pipeline should KEEP (predicted true variant).
# Verified equal to the original notebook (Hybrid_Model.py) output.
EXPECTED_GATK_KEEP = [100]
EXPECTED_DV_KEEP = [100, 300]


def _kept_positions(df):
    # Column 1 is POS in the headerless VCF dataframe.
    return sorted(int(p) for p in df[1].tolist())


@pytest.fixture(scope="module")
def gatk_clf():
    return load_model(package_resource_path("GATK_RF.pkl", "models"))


@pytest.fixture(scope="module")
def dv_clf():
    return load_model(package_resource_path("DV_RF.pkl", "models"))


def test_models_load(gatk_clf, dv_clf):
    # Guards against the scikit-learn version-skew failure mode.
    assert gatk_clf is not None
    assert dv_clf is not None


def test_gatk_keeps_expected_variants(gatk_clf):
    kept = format_gatk(str(FIXTURES / "test_gatk.vcf"), gatk_clf)
    assert _kept_positions(kept) == EXPECTED_GATK_KEEP


def test_dv_keeps_expected_variants(dv_clf):
    kept = format_dv(str(FIXTURES / "test_dv.vcf"), dv_clf)
    assert _kept_positions(kept) == EXPECTED_DV_KEEP
