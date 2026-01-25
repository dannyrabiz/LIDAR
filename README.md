# LIDAR: Learning-based Integration of DeepVariant and Bayesian Inference for Rare Variant Detection

## Overview

**LIDAR** is a machine-learning–based ensemble framework for **robust germline variant detection** in whole-exome sequencing (WES) and whole-genome sequencing (WGS) data. It integrates outputs from two state-of-the-art germline variant callers—**GATK HaplotypeCaller** and **DeepVariant**—using stacked Random Forest meta-models to improve sensitivity and specificity for **rare pathogenic variants**, particularly in **non-ideal clinical samples**.

While modern germline variant callers achieve high overall accuracy (>99%) on typical WES/WGS data, their performance degrades substantially for **rare variants**, **low-coverage samples**, and **heterogeneous clinical cohorts**. LIDAR addresses this gap by learning caller-specific error modes and selectively rescuing true variants missed by individual pipelines.

The framework was developed and validated using a large cohort of patients with metastatic prostate cancer and is designed to generalize to real-world clinical sequencing data.

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


Clone the repository and install the Python environment:

```bash
git clone https://github.com/<org>/lidar.git
cd lidar
conda env create -f environment.yml
conda activate lidar
```
---

## Requirements

- Python ≥ 3.8  
- GATK ≥ 4.x  
- DeepVariant ≥ 1.x  
- GLnexus (for joint genotyping, if applicable)

Python dependencies such as `scikit-learn`, `numpy`, and `pandas` should be installed in your environment.

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

```python -m lidar.lidar \
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

- `lidar.py` – Main LIDAR inference script  
- `models/` – Pretrained Random Forest classifiers  
- `utils/` – Helper functions for feature extraction and VCF handling  
- `README.md` – This file  

---

## Citation

If you use LIDAR in your work, please cite:

> *[Manuscript title]*  
> Rabizadeh et al., *in preparation*

(A DOI and final citation will be added upon publication.)

---

## Contact

For questions, feedback, or collaboration inquiries:

📧 **drabiza1@jhu.edu**

---

## License

This project is licensed under the **MIT License**.  
See the [LICENSE](LICENSE.md) file for details.
