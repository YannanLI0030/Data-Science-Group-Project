# Project data access

The raw biological datasets are not stored in Git because of their size. The
final application and all reported experiments use local copies of the project
data; AWS or other cloud access is not required.

## Local data layout

Place the supplied data under `data_s3/` and preserve the original
subdirectory and file names. The expected files and their SHA-256 checksums are
recorded in `data_manifest/manifest.json`.

The repository contains the application code, reviewed benchmarks, experiment
configurations and frozen result tables. Raw datasets and generated caches
remain outside version control.

A local installation should follow this structure:

```text
data_s3/
├── gene expression/
├── gene properties/
├── nomenclature/
└── non gene expression/
```

## Data verification

To verify a local file on macOS or Linux, run:

```bash
shasum -a 256 "data_s3/nomenclature/9_DepMap_sample_info.csv"
```

On Windows, run:

```powershell
certutil -hashfile "data_s3\nomenclature\9_DepMap_sample_info.csv" SHA256
```

Compare the result with the corresponding `raw_sha256` entry in
`data_manifest/manifest.json`. A mismatch indicates that the file is incomplete
or comes from a different data snapshot.

## Running locally

After placing the data under `data_s3/`, install the dependencies and start the
application from the repository root:

```bash
python -m pip install -r requirements.txt
python start.py
```

The command-line workflow is also available:

```bash
python start.py --cli \
  --target_gene EGFR \
  --exclusion_gene ABCB1 \
  --disease "cervical cancer" \
  --non-interactive
```

See `README.md` for the supported query modes and additional examples.

## Historical note

An earlier development workflow used a shared S3 archive to distribute the raw
files between team members. The final application and published experiments no
longer depend on that service. The former bucket, IAM and credential
instructions are therefore omitted because they are not required to reproduce
the submitted system and may no longer reflect the current cloud
configuration.
