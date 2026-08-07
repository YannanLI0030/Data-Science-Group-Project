#!/usr/bin/env python3
"""
fix_manifest.py — repair the stray brace and restructure provenance fields.

Usage:
    python fix_manifest.py staging_gz/manifest.json
"""
import json, re, sys, shutil
from pathlib import Path

# Canonical upstream source per file id. Portal-level URLs only: deep links
# rot, and none of these were downloaded directly anyway.
UPSTREAM = {
    "1_4_hpa_rna_celline.tsv": (
        "Human Protein Atlas, RNA cell line data",
        "https://www.proteinatlas.org/about/download",
        "HPA release not stated by provider - confirm with supervisor"),
    "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv": (
        "DepMap Omics expression, TPM log2(x+1), profile level",
        "https://depmap.org/portal/data_page/?tab=allData",
        "DepMap release not stated by provider - confirm with supervisor"),
    "3_GEOexpression.txt": (
        "NCBI GEO expression matrix",
        "https://www.ncbi.nlm.nih.gov/geo/",
        "GSE accession(s) not stated by provider - confirm with supervisor"),
    "4_Harmonized_MS_CCLE_Gygi_subsetted.csv": (
        "CCLE quantitative proteomics (Gygi lab), harmonised and subsetted",
        "https://depmap.org/portal/data_page/?tab=allData",
        "Harmonisation and subsetting performed by provider, not by this project"),
    "5_OmicsFusionFilteredSupplementary.csv": (
        "DepMap Omics fusion calls, filtered",
        "https://depmap.org/portal/data_page/?tab=allData",
        "DepMap release not stated by provider - confirm with supervisor"),
    "6_OmicsSomaticMutationsProfile.csv": (
        "DepMap Omics somatic mutations, profile level",
        "https://depmap.org/portal/data_page/?tab=allData",
        "DepMap release not stated by provider - confirm with supervisor"),
    "7_cellosaurus.csv": (
        "Cellosaurus cell line knowledge resource",
        "https://www.cellosaurus.org/",
        "Cellosaurus version not stated by provider - confirm with supervisor"),
    "8_DepMap_OmicsProfiles.csv": (
        "DepMap Omics profiles, ProfileID to ModelID bridge",
        "https://depmap.org/portal/data_page/?tab=allData",
        "DepMap release not stated by provider - confirm with supervisor"),
    "9_DepMap_sample_info.csv": (
        "DepMap sample / model annotation table",
        "https://depmap.org/portal/data_page/?tab=allData",
        "Filename follows pre-22Q2 naming (sample_info); may have been renamed by provider"),
    "10_GEOInfo.txt": (
        "NCBI GEO sample metadata, GSM to Cellosaurus mapping",
        "https://www.ncbi.nlm.nih.gov/geo/",
        "Mapping may have been assembled by provider - confirm with supervisor"),
    "11_hpa_rna_celline_description.tsv": (
        "Human Protein Atlas cell line descriptions",
        "https://www.proteinatlas.org/about/download",
        "HPA release not stated by provider - confirm with supervisor"),
    "12_CCLE_metabolomics_20190502.csv": (
        "CCLE metabolomics, 2019-05-02 release",
        "https://depmap.org/portal/data_page/?tab=allData",
        "Release date embedded in filename: 2019-05-02"),
    "13_CCLE_miRNA_20181103.gct": (
        "CCLE miRNA expression, 2018-11-03 release",
        "https://depmap.org/portal/data_page/?tab=allData",
        "Release date embedded in filename: 2018-11-03"),
    "14_OmicsGlobalSignatures.csv": (
        "DepMap Omics global genomic signatures (MSI, CIN, ploidy)",
        "https://depmap.org/portal/data_page/?tab=allData",
        "DepMap release not stated by provider - confirm with supervisor"),
}

ACQUIRED_FROM = ("University of Bristol OneDrive, shared by project supervisor "
                 "(AZ_project_data.zip, 2026-DSM-AZ folder)")
ACQUIRED_DATE = "2026-05-28"


def repair_text(text: str) -> tuple[str, int]:
    """Remove a stray '{' that splits one entry into two."""
    pattern = re.compile(r'("source_note":\s*"[^"]*",)\s*\n\s*\{\s*\n(\s*"download_url")')
    fixed, n = pattern.subn(r'\1\n\2', text)
    return fixed, n


def main() -> int:
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "staging_gz/manifest.json")
    text = path.read_text(encoding="utf-8")

    try:
        json.loads(text)
        print("JSON already parses, no brace repair needed.")
    except json.JSONDecodeError as exc:
        print(f"JSON broken at line {exc.lineno}: {exc.msg}")
        text, n = repair_text(text)
        print(f"Removed {n} stray brace(s).")
        json.loads(text)  # raises if still broken
        print("Repaired and parses cleanly.")

    shutil.copy(path, path.with_suffix(".json.bak"))
    print(f"Backup written to {path.with_suffix('.json.bak')}")

    m = json.loads(text)
    m["provenance_note"] = (
        "All 14 files were provided by the project supervisor as a single "
        "archive rather than downloaded individually. upstream_source records "
        "the canonical origin of each dataset; acquired_from records how this "
        "project obtained it. Where a release version is not stated, this is "
        "recorded explicitly rather than guessed."
    )

    unknown = []
    for e in m["files"]:
        src, url, note = UPSTREAM.get(
            e["file"], ("UNKNOWN", "UNKNOWN", "not in lookup table"))
        if src == "UNKNOWN":
            unknown.append(e["file"])
        e.pop("download_url", None)
        e.pop("download_date", None)
        e["upstream_source"] = src
        e["upstream_url"] = url
        e["version_note"] = note
        e["acquired_from"] = ACQUIRED_FROM
        e["acquired_date"] = ACQUIRED_DATE

    path.write_text(json.dumps(m, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nRewrote {path} with {len(m['files'])} entries.")
    if unknown:
        print("Files not found in lookup table:", unknown)

    todo = [e["file"] for e in m["files"] if "confirm with supervisor" in e["version_note"]]
    print(f"\n{len(todo)} files still need a release version from the supervisor:")
    for f in todo:
        print(f"  - {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
