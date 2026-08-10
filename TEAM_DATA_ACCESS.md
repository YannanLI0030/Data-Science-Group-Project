# Accessing the project data from S3

The 14 raw data files live in S3 so that everyone works from an identical,
checksummed copy. This page is for team members who need to run the scoring
code locally.

**Bucket:** `s3://celllineselector-team39/raw/v1_2026-08/`
**Region:** `eu-west-2` (London)
**Access:** read-only. Only Yannan can write to `raw/`, by design - that layer
is an immutable reference copy.

---

## One-time setup

### 1. Get your own AWS credentials

Sign in to the AWS console, then:

    IAM -> Users -> (your user) -> Security credentials -> Create access key
    -> Command Line Interface (CLI)

The secret is shown **once**. Save it to a password manager immediately.

Never accept a key that someone else generated for you, and never paste a
secret key into Slack, WhatsApp, email, or a chat with an AI assistant. If a
secret is ever exposed, delete that key and create a new one.

If you cannot see the `Create access key` button, you do not have permission
yet - ask Yannan to add you to the `celllineselector-readers` group.

### 2. Install and configure the tools

    pip install awscli boto3 s3fs
    aws configure

Fill in:

    AWS Access Key ID:     (yours)
    AWS Secret Access Key: (yours)
    Default region name:   eu-west-2
    Default output format: json

Check it worked:

    aws sts get-caller-identity
    aws s3 ls s3://celllineselector-team39/raw/v1_2026-08/

The first command should print your own user ARN. The second should list four
folders plus `manifest.json`.

### 3. Pull the data down once

Do **not** stream from S3 on every run. One of the files is 908 MB
uncompressed, and `data_loader.py` re-downloads on each call. Sync once
instead:

    cd <repo root>
    aws s3 sync s3://celllineselector-team39/raw/v1_2026-08/ ./data_s3/
    gunzip -r ./data_s3/

That transfers about 925 MB and expands to roughly 3.8 GB, so make sure you
have around 5 GB free. Add `data_s3/` to `.gitignore` if it is not already
covered.

---

## Running the code

`data_loader.py` takes its root path from the `CLS_DATA_ROOT` environment
variable, falling back to `data/`:

    # macOS / Linux
    export CLS_DATA_ROOT=./data_s3

    # Windows PowerShell
    $env:CLS_DATA_ROOT = "./data_s3"

Then, as usual:

    import sys; sys.path.insert(0, "src")
    from data_loader import CellLineDataLoader

    loader = CellLineDataLoader()
    print(loader.sample_info.shape)      # expect (1840, 29)

Note that `sample_info` is a property, so there are no parentheses.

Nothing about the existing API changed. `CellLineDataLoader("data")` still
reads from a local folder exactly as before, so code written against the old
version keeps working.

### Reading straight from S3 (occasional use only)

Setting `CLS_DATA_ROOT` to an `s3://` URI also works and needs no local copy:

    export CLS_DATA_ROOT='s3://celllineselector-team39/raw/v1_2026-08'

The loader then reads the `.gz` objects directly and maps the local folder
names (which contain spaces) onto the S3 keys (which use underscores). This is
convenient for a quick check on a small file, but slow for anything large. Use
the synced local copy for real work.

---

## Verifying you got the right data

`data_manifest/manifest.json` in this repo records, for every file: sha256 of
both the raw and gzipped form, byte size, physical line count, the upstream
source, and how the project obtained it.

To confirm a local file matches the archive:

    shasum -a 256 data_s3/nomenclature/9_DepMap_sample_info.csv     # macOS/Linux
    certutil -hashfile <path> SHA256                                 # Windows

Compare against `raw_sha256` for that entry in the manifest. If they differ,
the file is truncated or from a different version - re-sync rather than
debugging downstream.

---

## Notes on provenance

All 14 files were supplied by the project supervisor as a single archive
rather than downloaded individually, so release versions are not stated for
most of them. The manifest records this explicitly (`version_note`) instead of
guessing. If anyone confirms the DepMap release with Daniel or Luigi, update
the manifest and re-upload it.

---

## Cost

Storage is about 925 MiB, roughly USD 0.02 per month. Egress within the free
tier covers normal use. Syncing the full dataset repeatedly is the only thing
that would show up on a bill, which is another reason to sync once and keep
the local copy.
