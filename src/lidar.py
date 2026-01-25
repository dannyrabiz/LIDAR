#!/usr/bin/env python3
import argparse
import pickle
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from importlib import resources as importlib_resources  # py>=3.9
except Exception:
    importlib_resources = None

vcf_header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]


def read_vcf_as_df(vcf_path):
    return pd.read_csv(vcf_path, sep="\t", comment="#", header=None)


def parse_info(info_str):
    if not isinstance(info_str, str) or info_str == ".":
        return {}
    out = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = True
    return out


def safe_float(x, default=0.0):
    try:
        if x is True:
            return 1.0
        if x is False or x is None:
            return default
        return float(x)
    except Exception:
        return default


def normalize_gt(gt):
    if not isinstance(gt, str) or gt == ".":
        return 0.0
    gt = gt.split(":")[0]
    mapping = {
        "0/1": 0.0,
        "1/0": 0.0,
        "0|1": 0.0,
        "1|0": 0.0,
        "1/1": 1.0,
        "1|1": 1.0,
        "./1": 1.0,
        "1/.": 1.0,
        "./.": 0.0,
        "0/0": 0.0,
        "0|0": 0.0,
    }
    return mapping.get(gt, 0.0)


def run_cmd_bash(command):
    subprocess.run(["bash", "-lc", command], check=True)


def package_resource_path(filename, subdir):
    pkg = "lidar"

    if importlib_resources is not None:
        try:
            p = importlib_resources.files(pkg).joinpath("resources", subdir, filename)
            return Path(p)
        except Exception:
            pass

    here = Path(__file__).resolve().parent
    candidates = [
        here / "resources" / subdir / filename,
        here.parent / "resources" / subdir / filename,
        here.parent.parent / "resources" / subdir / filename,
    ]
    for c in candidates:
        if c.exists():
            return c

    raise FileNotFoundError(
        f"Could not locate resource '{filename}'. Provide an explicit path via CLI (e.g., --dv-model/--gatk-model/--header-path)."
    )


def format_dv(vcf_path, dv_clf):
    df = read_vcf_as_df(vcf_path)
    fmt = df[9].astype(str).str.split(":", expand=True)

    gt = fmt[0] if 0 in fmt.columns else ""
    dp = fmt[1] if 1 in fmt.columns else "0"
    ad_field = fmt[2] if 2 in fmt.columns else "0,0"
    gq = fmt[3] if 3 in fmt.columns else "0"

    ad_split = ad_field.astype(str).str.split(",", expand=True)
    rd = ad_split[0] if 0 in ad_split.columns else "0"
    ad = ad_split[1] if 1 in ad_split.columns else "0"

    info_dicts = df[7].astype(str).apply(parse_info)
    aq = info_dicts.apply(lambda d: safe_float(d.get("AQ", 0.0), 0.0))

    train = pd.DataFrame()
    train["SNP_INDEL"] = df[3].astype(str).apply(len) - df[4].astype(str).apply(len)
    train["Qual"] = pd.to_numeric(df[5], errors="coerce").fillna(0.0)
    train["AQ"] = aq

    train["GT"] = gt.apply(normalize_gt)
    train["DP"] = pd.to_numeric(dp, errors="coerce").fillna(0.0)
    train["GQ"] = pd.to_numeric(gq, errors="coerce").fillna(0.0)
    train["RD"] = pd.to_numeric(rd, errors="coerce").fillna(0.0)
    train["AD"] = pd.to_numeric(ad, errors="coerce").fillna(0.0)

    train["VAF"] = np.where(train["DP"] > 0, train["AD"] / train["DP"], 0.0)
    train["Caller"] = 1.0

    train = train.replace([np.inf, -np.inf], 0.0).fillna(0.0).astype(float)

    results = dv_clf.predict(train.to_numpy())
    return df.loc[np.where(results == 0)[0]]


def format_gatk(vcf_path, gatk_clf):
    df = read_vcf_as_df(vcf_path)
    fmt = df[9].astype(str).str.split(":", expand=True)

    gt = fmt[0] if 0 in fmt.columns else ""
    ad_field = fmt[1] if 1 in fmt.columns else "0,0"
    dp = fmt[2] if 2 in fmt.columns else "0"
    gq = fmt[3] if 3 in fmt.columns else "0"

    ad_split = ad_field.astype(str).str.split(",", expand=True)
    rd = ad_split[0] if 0 in ad_split.columns else "0"
    ad = ad_split[1] if 1 in ad_split.columns else "0"

    info_dicts = df[7].astype(str).apply(parse_info)
    info_str = df[7].astype(str)

    train = pd.DataFrame()
    train["SNP_INDEL"] = df[3].astype(str).apply(len) - df[4].astype(str).apply(len)
    train["Filter"] = (df[6].astype(str) != "PASS").astype(float)
    train["Negative"] = info_str.str.contains("NEGATIVE", na=False).replace(True, 0).replace(False, 1).astype(float)
    train["Positive"] = info_str.str.contains("POSITIVE", na=False).replace(True, 0).replace(False, 1).astype(float)
    train["Qual"] = pd.to_numeric(df[5], errors="coerce").fillna(0.0)
    train["SOR"] = info_dicts.apply(lambda d: safe_float(d.get("SOR", 0.0), 0.0))
    train["FS"] = info_dicts.apply(lambda d: safe_float(d.get("FS", 0.0), 0.0))
    train["QD"] = info_dicts.apply(lambda d: safe_float(d.get("QD", 0.0), 0.0))
    train["MQRANK"] = info_dicts.apply(lambda d: safe_float(d.get("MQRankSum", 0.0), 0.0))

    train["GT"] = gt.apply(normalize_gt)
    train["DP"] = pd.to_numeric(dp, errors="coerce").fillna(0.0)
    train["GQ"] = pd.to_numeric(gq, errors="coerce").fillna(0.0)
    train["RD"] = pd.to_numeric(rd, errors="coerce").fillna(0.0)
    train["AD"] = pd.to_numeric(ad, errors="coerce").fillna(0.0)

    train["VAF"] = np.where(train["DP"] > 0, train["AD"] / train["DP"], 0.0)
    train["Caller"] = 1.0

    train = train.replace([np.inf, -np.inf], 0.0).fillna(0.0).astype(float)

    results = gatk_clf.predict(train.to_numpy())
    return df.loc[np.where(results == 0)[0]]


def LIDAR(
    gatk_samp,
    dv_samp,
    workdir,
    output_file,
    ref_sdf,
    header_path=None,
    sample_placeholder="46A",
    sample_id=None,
    dv_model_path=None,
    gatk_model_path=None,
    conda_env_cmd="ml anaconda; conda activate bio;",
    keep_workdir=False,
):
    workdir_path = Path(workdir)

    if header_path is None:
        header_path = package_resource_path("FULL_GATK_Header.txt", "headers")
    else:
        header_path = Path(header_path)

    if dv_model_path is None:
        dv_model_path = package_resource_path("DV_RF.pkl", "models")
    else:
        dv_model_path = Path(dv_model_path)

    if gatk_model_path is None:
        gatk_model_path = package_resource_path("GATK_RF.pkl", "models")
    else:
        gatk_model_path = Path(gatk_model_path)

    with open(dv_model_path, "rb") as f:
        dv_clf = pickle.load(f)
    with open(gatk_model_path, "rb") as f:
        gatk_clf = pickle.load(f)

    ref_sdf = Path(ref_sdf)
    if not ref_sdf.exists():
        raise FileNotFoundError(f"RTG SDF not found: {ref_sdf}")

    cmd = (
        f"{conda_env_cmd}"
        f"rtg vcfeval -b {gatk_samp} -c {dv_samp} "
        f"--all-records --squash-ploidy --no-roc "
        f"-t {ref_sdf} -o {workdir_path}"
    )
    run_cmd_bash(cmd)

    gatk_df = format_gatk(workdir_path / "fn.vcf.gz", gatk_clf)
    dv_df = format_dv(workdir_path / "fp.vcf.gz", dv_clf)
    both_df = read_vcf_as_df(workdir_path / "tp-baseline.vcf.gz")

    all_df = pd.concat([both_df, gatk_df, dv_df], axis=0, ignore_index=True)
    all_df.columns = vcf_header
    all_df = all_df.sort_values(["CHROM", "POS"])

    if sample_id is None:
        sample_id = workdir_path.name

    header_text = header_path.read_text().rstrip("\n")
    header_text = header_text.replace(sample_placeholder, sample_id)

    output_vcf = Path(f"{output_file}.vcf")
    with open(output_vcf, "w") as vcf:
        vcf.write(header_text + "\n")
    all_df.to_csv(output_vcf, sep="\t", mode="a", index=False, header=None)

    if not keep_workdir:
        for p in sorted(workdir_path.glob("**/*"), reverse=True):
            try:
                p.unlink()
            except IsADirectoryError:
                p.rmdir()
        try:
            workdir_path.rmdir()
        except Exception:
            pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gatk_samp")
    parser.add_argument("dv_samp")
    parser.add_argument("workdir")
    parser.add_argument("output_file")
    parser.add_argument("--header-path", default=None)
    parser.add_argument("--sample-placeholder", default="46A")
    parser.add_argument("--sample-id", default=None)
    parser.add_argument("--dv-model", default=None)
    parser.add_argument("--gatk-model", default=None)
    parser.add_argument(  "--ref-sdf",required=True,
    help="Path to RTG reference SDF (create with: rtg format -o <ref>.sdf <ref>.fa)"  )
    parser.add_argument("--conda-env-cmd", default="ml anaconda; conda activate bio;")
    parser.add_argument("--keep-workdir", action="store_true")

    args = parser.parse_args()

    LIDAR(
        args.gatk_samp,
        args.dv_samp,
        args.workdir,
        args.output_file,
        header_path=args.header_path,
        sample_placeholder=args.sample_placeholder,
        sample_id=args.sample_id,
        dv_model_path=args.dv_model,
        gatk_model_path=args.gatk_model,
        ref_sdf=args.ref_sdf,
        conda_env_cmd=args.conda_env_cmd,
        keep_workdir=args.keep_workdir,
    )


if __name__ == "__main__":
    main()
