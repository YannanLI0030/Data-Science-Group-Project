#!/usr/bin/env python3
"""
aws_prepare_and_upload.py
=========================
CellLineSelector - Prepare local data and upload to S3.

Local folder names keep their spaces (so data_loader.py needs no changes).
S3 keys use underscores (so the CLI and boto3 need no quoting gymnastics).
The mapping between the two lives in FOLDER_MAP below.

What it does, in order:
    1. Verify all 14 expected files exist locally (report anything missing).
    2. Gzip each file into a staging directory (skips files already staged).
    3. Compute sha256 / size / line count for every file -> manifest.json
    4. Upload staging directory + manifest to s3://<bucket>/raw/<version>/
    5. Verify by listing what actually landed in S3.

Usage
-----
    pip install boto3

    # See what is actually on disk, change nothing:
    python aws_prepare_and_upload.py --data-dir ./data --inspect

    # Verify, gzip, build manifest, but do not upload:
    python aws_prepare_and_upload.py --data-dir ./data --bucket YOUR_BUCKET --dry-run

    # Full run:
    python aws_prepare_and_upload.py --data-dir ./data --bucket YOUR_BUCKET

    # Re-running is safe: staging and upload both skip work already done.

Credentials come from `aws configure` (~/.aws/credentials) or the standard
AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY environment variables.
NEVER hard-code them in this file.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

# --------------------------------------------------------------------------
# Local folder name  ->  S3 folder name.
# Left side must match your disk exactly, spaces included.
# Right side is what appears in the S3 key.
# --------------------------------------------------------------------------

FOLDER_MAP: dict[str, str] = {
    "gene expression":     "gene_expression",
    "gene properties":     "gene_properties",
    "nomenclature":        "nomenclature",
    "non gene expression": "non_gene_expression",
}

# Folders that exist locally but must NOT go into raw/ (derived data).
IGNORE_FOLDERS = {"merged", "parquet", "staging_gz"}

# Expected filenames, keyed by LOCAL folder name.
EXPECTED_FILES: dict[str, list[str]] = {
    "gene expression": [
        "1_4_hpa_rna_celline.tsv",
        "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv",
        "3_GEOexpression.txt",
        "4_Harmonized_MS_CCLE_Gygi_subsetted.csv",
    ],
    "gene properties": [
        "5_OmicsFusionFilteredSupplementary.csv",
        "6_OmicsSomaticMutationsProfile.csv",
    ],
    "nomenclature": [
        "7_cellosaurus.csv",
        "8_DepMap_OmicsProfiles.csv",
        "9_DepMap_sample_info.csv",
        "10_GEOInfo.txt",
        "11_hpa_rna_celline_description.tsv",
    ],
    "non gene expression": [
        "12_CCLE_metabolomics_20190502.csv",
        "13_CCLE_miRNA_20181103.gct",
        "14_OmicsGlobalSignatures.csv",
    ],
}

# Provenance notes that go into the manifest. Fill in the real download URLs
# and dates before the report - this is the data governance record.
SOURCE_NOTES: dict[str, str] = {
    "gene_expression": "HPA RNA cell line, DepMap OmicsExpression TPM log2(x+1), GEO series matrix, CCLE Gygi harmonised MS proteomics",
    "gene_properties": "DepMap OmicsFusion filtered, DepMap OmicsSomaticMutations profile-level",
    "nomenclature": "Cellosaurus registry, DepMap OmicsProfiles bridge, DepMap sample info, GEO GSM metadata, HPA cell line descriptions",
    "non_gene_expression": "CCLE metabolomics 2019-05-02, CCLE miRNA 2018-11-03, DepMap OmicsGlobalSignatures",
}

GZIP_LEVEL = 6           # 6 is the sweet spot; 9 is much slower for ~2% gain
CHUNK = 8 * 1024 * 1024  # 8 MB read chunks


# --------------------------------------------------------------------------
# Inspect mode: report what is actually on disk, change nothing
# --------------------------------------------------------------------------

def inspect(data_dir: Path) -> int:
    if not data_dir.is_dir():
        print(f"Not a directory: {data_dir}")
        return 1

    print(f"Inspecting {data_dir.resolve()}\n")
    expected_lookup = {
        folder: set(names) for folder, names in EXPECTED_FILES.items()
    }

    all_ok = True
    for child in sorted(data_dir.iterdir()):
        if child.name.startswith("."):
            continue

        if child.is_file():
            print(f"[loose file, ignored]  {child.name}  "
                  f"({child.stat().st_size / 1e6:,.1f} MB)")
            continue

        if child.name in IGNORE_FOLDERS:
            n = sum(1 for _ in child.rglob("*") if _.is_file())
            print(f"[derived, not uploaded]  {child.name}/  ({n} files)")
            continue

        known = child.name in FOLDER_MAP
        s3name = FOLDER_MAP.get(child.name, "???")
        tag = f"-> s3: {s3name}/" if known else "NOT IN FOLDER_MAP"
        print(f"\n{child.name}/   {tag}")

        wanted = expected_lookup.get(child.name, set())
        seen = set()
        for f in sorted(child.iterdir()):
            if f.name.startswith(".") or not f.is_file():
                continue
            seen.add(f.name)
            mark = "ok " if f.name in wanted else "EXTRA"
            print(f"   {mark}  {f.stat().st_size / 1e6:9,.1f} MB  {f.name}")

        for missing in sorted(wanted - seen):
            print(f"   MISS   {'':>9}     {missing}")
            all_ok = False

        if not known:
            all_ok = False

    print()
    if all_ok:
        print("Layout matches expectations. Safe to run --dry-run next.")
    else:
        print("Mismatches above. Either rename the files on disk, or edit")
        print("EXPECTED_FILES / FOLDER_MAP at the top of this script to match.")
    return 0


# --------------------------------------------------------------------------
# Step 1: verify
# --------------------------------------------------------------------------

def verify_layout(data_dir: Path) -> tuple[list[tuple[Path, PurePosixPath]], list[str]]:
    """Return ([(source_path, s3_relative_path)], missing_descriptions)."""
    found: list[tuple[Path, PurePosixPath]] = []
    missing: list[str] = []

    for local_folder, filenames in EXPECTED_FILES.items():
        s3_folder = FOLDER_MAP[local_folder]
        folder = data_dir / local_folder
        for name in filenames:
            path = folder / name
            if path.is_file():
                found.append((path, PurePosixPath(s3_folder) / name))
            else:
                missing.append(f"{local_folder}/{name}")

    return found, missing


# --------------------------------------------------------------------------
# Step 2: gzip into staging
# --------------------------------------------------------------------------

def gzip_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".partial")
    with open(src, "rb") as fin, gzip.open(tmp, "wb", compresslevel=GZIP_LEVEL) as fout:
        shutil.copyfileobj(fin, fout, length=CHUNK)
    tmp.replace(dst)  # atomic: a killed run never leaves a valid-looking stub


def stage_all(found, staging: Path) -> list[Path]:
    staged: list[Path] = []
    total = len(found)

    for i, (src, rel) in enumerate(found, start=1):
        dst = staging / (str(rel) + ".gz")

        if dst.is_file() and dst.stat().st_mtime >= src.stat().st_mtime:
            print(f"  [{i}/{total}] skip (already staged)  {rel}")
            staged.append(dst)
            continue

        src_mb = src.stat().st_size / 1e6
        print(f"  [{i}/{total}] gzip {rel}  ({src_mb:,.0f} MB) ...", end="", flush=True)
        gzip_file(src, dst)
        dst_mb = dst.stat().st_size / 1e6
        ratio = src_mb / dst_mb if dst_mb else 0
        print(f" -> {dst_mb:,.0f} MB  ({ratio:.1f}x)")
        staged.append(dst)

    return staged


# --------------------------------------------------------------------------
# Step 3: manifest
# --------------------------------------------------------------------------

def sha256_and_lines(path: Path, count_lines: bool) -> tuple[str, int | None]:
    """One pass over the file for both the checksum and the line count."""
    h = hashlib.sha256()
    lines = 0
    with open(path, "rb") as f:
        while chunk := f.read(CHUNK):
            h.update(chunk)
            if count_lines:
                lines += chunk.count(b"\n")
    return h.hexdigest(), (lines if count_lines else None)


def build_manifest(found, staged, version: str, bucket: str) -> dict:
    entries = []
    total = len(found)

    for i, ((src, rel), gz) in enumerate(zip(found, staged), start=1):
        print(f"  [{i}/{total}] checksum {rel} ...", end="", flush=True)

        raw_sha, raw_lines = sha256_and_lines(src, count_lines=True)
        gz_sha, _ = sha256_and_lines(gz, count_lines=False)
        s3_folder = rel.parent.name

        entries.append({
            "file": rel.name,
            "local_path": str(src),
            "s3_folder": s3_folder,
            "s3_key": f"raw/{version}/{rel}.gz",
            "raw_bytes": src.stat().st_size,
            "raw_sha256": raw_sha,
            "raw_line_count": raw_lines,
            "gz_bytes": gz.stat().st_size,
            "gz_sha256": gz_sha,
            "source_note": SOURCE_NOTES.get(s3_folder, ""),
            "download_url": "TODO: fill in before the report",
            "download_date": "TODO: fill in before the report",
        })
        print(" ok")

    return {
        "dataset": "CellLineSelector raw data",
        "version": version,
        "bucket": bucket,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "file_count": len(entries),
        "total_raw_bytes": sum(e["raw_bytes"] for e in entries),
        "total_gz_bytes": sum(e["gz_bytes"] for e in entries),
        "note": (
            "raw/ is immutable. Never overwrite these objects. Derived parquet "
            "belongs under curated/. Local folder names contain spaces; S3 keys "
            "use underscores. Line counts are physical newline counts and may "
            "exceed logical record counts where fields contain embedded newlines."
        ),
        "files": entries,
    }


# --------------------------------------------------------------------------
# Steps 4 and 5: upload and verify
# --------------------------------------------------------------------------

def upload_all(staged: list[Path], staging: Path, manifest_path: Path,
               bucket: str, version: str) -> None:
    import boto3
    from boto3.s3.transfer import TransferConfig
    from botocore.exceptions import ClientError

    s3 = boto3.client("s3")
    # Multipart above 64 MB with 8 parallel threads. Without this, a 1 GB
    # single-part PUT has to restart from zero if the connection drops.
    cfg = TransferConfig(
        multipart_threshold=64 * 1024 * 1024,
        multipart_chunksize=64 * 1024 * 1024,
        max_concurrency=8,
        use_threads=True,
    )

    jobs = [(gz, f"raw/{version}/{gz.relative_to(staging).as_posix()}") for gz in staged]
    jobs.append((manifest_path, f"raw/{version}/manifest.json"))

    total = len(jobs)
    for i, (path, key) in enumerate(jobs, start=1):
        size = path.stat().st_size

        try:
            head = s3.head_object(Bucket=bucket, Key=key)
            if head["ContentLength"] == size:
                print(f"  [{i}/{total}] skip (already in S3)  {key}")
                continue
        except ClientError as exc:
            if exc.response["Error"]["Code"] not in ("404", "NoSuchKey"):
                raise

        done = [0]

        def progress(n: int, _done=done, _size=size, _key=key) -> None:
            _done[0] += n
            pct = 100 * _done[0] / _size if _size else 100
            print(f"\r  [{i}/{total}] {_key}  {pct:5.1f}%", end="", flush=True)

        s3.upload_file(str(path), bucket, key, Config=cfg, Callback=progress)
        print()


def verify_remote(bucket: str, version: str, expected: int) -> None:
    import boto3

    s3 = boto3.client("s3")
    paginator = s3.get_paginator("list_objects_v2")
    objects = []
    for page in paginator.paginate(Bucket=bucket, Prefix=f"raw/{version}/"):
        objects.extend(page.get("Contents", []))

    total_bytes = sum(o["Size"] for o in objects)
    print(f"\nS3 now holds {len(objects)} objects, {total_bytes / 1e9:.2f} GB")
    for o in sorted(objects, key=lambda x: x["Key"]):
        print(f"  {o['Size'] / 1e6:9,.1f} MB  {o['Key']}")

    if len(objects) != expected:
        print(f"\nWARNING: expected {expected} objects, found {len(objects)}.")


# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=Path("data"),
                    help="Local folder holding the four data subfolders (default: ./data)")
    ap.add_argument("--staging", type=Path, default=Path("staging_gz"),
                    help="Where gzipped copies are written (default: ./staging_gz)")
    ap.add_argument("--bucket", help="Target S3 bucket name (not needed for --inspect)")
    ap.add_argument("--version", default="v1_2026-08",
                    help="Dataset version folder under raw/ (default: v1_2026-08)")
    ap.add_argument("--inspect", action="store_true",
                    help="List what is actually on disk and exit. Changes nothing.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Verify, gzip and build the manifest, but do not upload")
    args = ap.parse_args()

    if args.inspect:
        return inspect(args.data_dir)

    if not args.bucket:
        ap.error("--bucket is required unless you pass --inspect")

    print("=" * 70)
    print("Step 1/5  Verifying local layout")
    print("=" * 70)
    found, missing = verify_layout(args.data_dir)
    print(f"  found {len(found)} of 14 expected files")
    if missing:
        print("\n  MISSING:")
        for m in missing:
            print(f"    - {m}")
        print("\n  Run with --inspect to see what is actually there.")
        return 1

    print("\n" + "=" * 70)
    print("Step 2/5  Gzipping into staging")
    print("=" * 70)
    staged = stage_all(found, args.staging)

    print("\n" + "=" * 70)
    print("Step 3/5  Building manifest (sha256 + line counts)")
    print("=" * 70)
    manifest = build_manifest(found, staged, args.version, args.bucket)
    manifest_path = args.staging / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    raw_gb = manifest["total_raw_bytes"] / 1e9
    gz_gb = manifest["total_gz_bytes"] / 1e9
    print(f"\n  raw {raw_gb:.2f} GB  ->  gz {gz_gb:.2f} GB  ({raw_gb / gz_gb:.1f}x smaller)")
    print(f"  manifest written to {manifest_path}")

    if args.dry_run:
        print("\n--dry-run set: stopping before upload.")
        return 0

    print("\n" + "=" * 70)
    print(f"Step 4/5  Uploading to s3://{args.bucket}/raw/{args.version}/")
    print("=" * 70)
    upload_all(staged, args.staging, manifest_path, args.bucket, args.version)

    print("\n" + "=" * 70)
    print("Step 5/5  Verifying remote")
    print("=" * 70)
    verify_remote(args.bucket, args.version, expected=len(staged) + 1)

    print("\nDone. raw/ is now the immutable reference copy.")
    print("Keep your local files: S3 is the archive, not your working directory.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
