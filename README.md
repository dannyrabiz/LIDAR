# LIDAR: Learning-based Integration of DeepVariant and Bayesian Inference for Rare Variant Detection

## Overview

**LIDAR** is a machine-learning–based ensemble framework for **robust germline variant detection** in whole-exome sequencing (WES) and whole-genome sequencing (WGS) data. It integrates outputs from two state-of-the-art germline variant callers—**GATK HaplotypeCaller** and **DeepVariant**—using stacked Random Forest meta-models to improve sensitivity and specificity for **rare pathogenic variants**, particularly in **non-ideal clinical samples**.

---
### Quickstart

```bash
lidar gatk.vcf.gz dv.vcf.gz work output --ref-sdf hg19.sdf
```
---

## Key Features

- Stacked ensemble architecture integrating Bayesian (GATK) and deep-learning (DeepVariant) variant callers  
- Caller-specific Random Forest models trained on rich per-variant annotations  
- Explicit handling of low-coverage and heterogeneous samples common in clinical sequencing  
- Rescue of rare pathogenic variants missed by standard pipelines  
- Compatible with both WES and WGS workflows  
- Modular design for downstream extension or retraining  

---

## Method Summary

1. Variants are independently called using **GATK HaplotypeCaller** and **DeepVariant**.
2. Variants detected by both callers are retained directly.
3. Variants detected by only one caller are evaluated using caller-specific Random Forest classifiers.
4. The final output is a merged VCF containing high-confidence germline variants with improved recall for rare pathogenic alleles.

---

## Installation

Clone the repository and install the package. Either of the following works.

**With pip:**

```bash
git clone https://github.com/dannyrabiz/LIDAR.git
cd LIDAR
pip install .
```

**With conda:**

```bash
git clone https://github.com/dannyrabiz/LIDAR.git
cd LIDAR
conda env create -f environment.yml
conda activate lidar
```

Both methods install the `lidar` command-line tool and the Python package.

---

## Requirements

- Python ≥ 3.8
- [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools) (`rtg`) on your `PATH` — used to reconcile variant representation between callers
- GATK ≥ 4.x and DeepVariant ≥ 1.x — used upstream to *produce* the input VCFs (not called by LIDAR itself)
- GLnexus (for joint genotyping, if applicable)

Python dependencies (`numpy`, `pandas`, `scikit-learn`) are installed automatically.

### Model compatibility

> ⚠️ The pretrained Random Forest models bundled with LIDAR were trained with
> **scikit-learn 1.2.x**. scikit-learn changed its internal tree format in 1.3,
> so the models will not unpickle under scikit-learn ≥ 1.3. The dependency pin
> (`scikit-learn>=1.2,<1.3`) handles this automatically; if you manage your own
> environment, install a compatible version. To use a newer scikit-learn you
> would need to retrain the models (see `Hybrid_Training_Testing.ipynb`).

---

## Reference Genome (RTG SDF)

LIDAR requires an RTG Tools reference SDF corresponding to the same genome
build used for variant calling (e.g., hg19 or GRCh38).

Create the SDF once using:

```bash
rtg format -o hg19.sdf hg19.fa
```
---

## Usage

Run LIDAR on a single sample:

```bash
lidar \
  sample.gatk.vcf.gz \
  sample.dv.vcf.gz \
  workdir \
  output_prefix \
  --ref-sdf /path/to/reference.sdf
```
Alternatively, without installing the CLI:

```
python -m lidar.lidar \
  sample.gatk.vcf.gz \
  sample.dv.vcf.gz \
  workdir \
  output_prefix \
  --ref-sdf /path/to/reference.sdf
```
---


## Input

1. GATK VCF – Variant calls generated using GATK best practices
2. DeepVariant VCF – Variant calls generated using DeepVariant (optionally via GLnexus)
3. SampleID – Unique identifier for the sample
4. Output VCF – Destination path for the LIDAR output

---

## Output

- A VCF file containing:
  - Variants called by both GATK and DeepVariant 
  - High-confidence variants rescued by the LIDAR ensemble
    
---

## Intended Use

LIDAR is intended for:

- Germline variant analysis in **clinical or translational sequencing studies**
- Cohorts with **variable sequencing depth or batch effects**
- Detection of **rare pathogenic or likely pathogenic variants**
- Complementing, not replacing, standard GATK/DeepVariant workflows

⚠️ **LIDAR is a research tool and is not intended for direct clinical diagnosis without appropriate validation.**

---

## Repository Contents

- `lidar/lidar.py` – Main LIDAR inference script and CLI
- `lidar/resources/models/` – Pretrained Random Forest classifiers (`GATK_RF.pkl`, `DV_RF.pkl`)
- `lidar/resources/headers/` – VCF header templates for the output
- `Hybrid_Training_Testing.ipynb` – Notebook used to build truth sets and train the models
- `VariantCallingComparison.ipynb` – Analysis/figure-generation notebook
- `tests/` – Regression tests verifying the package reproduces the original notebook logic
- `pyproject.toml` / `environment.yml` – Packaging and environment definitions

---

## Testing

```bash
pip install -e ".[test]"
pytest -q
```

The tests confirm the bundled models load and that the packaged inference keeps
exactly the variants the original Jupyter-notebook implementation does.

---

## Citation

If you use LIDAR in your work, please cite:

> Rabizadeh et al. Manuscript in preparation.

(A DOI and final citation will be added upon publication.)

---

## Contact

For questions, feedback, or collaboration inquiries:

📧 **drabiza1@jhu.edu**

---

## License

This project is licensed under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.
