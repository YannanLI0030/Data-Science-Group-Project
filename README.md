# CellLineSelector

CellLineSelector is a local multi-omics recommendation system for selecting
human cell lines by target gene, target Protein, optional exclusion gene and
disease or tissue context. It combines a deterministic ranker with an optional
post-ranking explanation Agent.

## Final application

`dynamic_cellline_selector_gene_protein.py` is the canonical ranking backend.
It loads the project datasets on demand, harmonises identifiers, applies the
disease filter and writes the ranked CSV. `api_server.py` exposes the same
backend to `web/index.html`; the UI does not contain a second scorer.

The Agent runs only after the ranking has been saved. It receives structured
evidence cards and cannot change candidate eligibility, scores, recommendation
levels or rank order. Claims without valid evidence IDs, and numerical claims
not supported by their cited evidence, are dropped before display. The
optional full mode also checks candidate-entity consistency. The offline
scripted explanation is the default and needs no API key.

Three query modes are supported:

- `GENE_MULTIOMICS`: target gene only;
- `PROTEIN_ONLY`: target Protein only;
- `COMBINED`: target gene and Protein.

For gene and combined queries, A0 assigns 0.55 to RNA, 0.30 to Protein and
0.15 to Confidence. Missing biological modalities are removed from the
biological denominator. Protein-only queries use 0.85 Protein and 0.15
Confidence. Confidence combines completeness (0.40), source support (0.35)
and RNA--Protein consistency (0.25). The exclusion penalty is capped at 0.30.

## Setup

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

The 3.6 GB raw dataset is not stored in Git. Place the supplied files under
`data_s3/` using their original four subdirectories. The expected filenames
and checksums are recorded in `data_manifest/manifest.json`.

Hosted explanation models are optional. Copy `.env.example` to `.env` only if
you want to configure an OpenAI or Anthropic key. Never commit `.env`. Keys
entered in the UI remain in page memory and are sent only to the local API for
the current request.

## Run the application

Start the local UI:

```bash
python start.py
```

The browser opens at `http://127.0.0.1:8000/`. To start the server without
opening a browser:

```bash
python api_server.py --no-browser
```

Run the command-line workflow:

```bash
python start.py --cli \
  --target_gene EGFR \
  --exclusion_gene ABCB1 \
  --disease "cervical cancer" \
  --non-interactive
```

Use `--target_protein` instead of `--target_gene` for Protein-only mode, or
supply both for a combined query. Add `--no-agent` when only the deterministic
ranking is required.

## Tests

The application tests do not require the full raw dataset:

```bash
python -m unittest discover -s tests -v
```

They cover the Agent boundary, citation and numerical grounding, Protein-only
scoring, and the UI/API query contract. The main ablation checks can be run
separately:

```bash
python scripts/ablation_runner.py --self-test
python scripts/dynamic_ablation_runner.py --self-test
python scripts/exclusion_penalty_ablation_runner.py --self-test
```

## Evaluation record

The evaluation stages answer different questions and should not be collapsed
into one “best model” claim:

| Stage | Purpose | Main boundary |
|---|---|---|
| V4 development | Generate component and weight hypotheses from 50 labels across 10 genes | Model-informed pool; only 2 verified negatives |
| Dynamic ablation | Check 14 interventions across two 100-gene panels | No relevance labels; structural sensitivity only |
| Frozen V5 | Compare 10 configurations on 12 reviewed genes | 99 positives and 90 unknowns; no verified negatives |
| Penalty V1 | Test P00/P15/P30/P45 mechanics on 5 fixed queries | No pair-specific Gold; no optimal-cap claim |

V5 did not reproduce the V4 advantage of removing direct Protein, so A0 was
retained. Removing the full Confidence block or source support caused clearer
NDCG@5 reductions, while other variants showed endpoint- or gene-specific
trade-offs. The Penalty experiment passed all 155 implementation checks and
showed that P30 changes rankings, but it did not establish 0.30 as biologically
optimal.

Detailed records are in:

- `docs/dynamic_ablation_unlabelled.md`;
- `docs/dynamic_ablation_v5_freeze.md`;
- `docs/v5_holdout_evaluation.md`;
- `docs/exclusion_penalty_structural_v1.md`;
- `docs/weight_search_v4_development.md`;
- `docs/experiment_provenance.md`.

## Repository map

```text
api_server.py                              local UI/API
dynamic_cellline_selector_gene_protein.py canonical dynamic ranker
start.py                                  UI and CLI entry point
src/data_loader.py                        local raw-data loader
src/data_merger.py                        identifier harmonisation
src/agentic/                              post-ranking explanation layer
web/index.html                            browser interface
tests/                                    application boundary tests
scripts/                                  ablation and benchmark tooling
benchmarks/                               reviewed pools and Gold files
results/                                  frozen experiment outputs
docs/                                     methods, limits and provenance
```

`src/merged_cellline_selector.py` is retained for the earlier merged-table
experiments. It is not the final UI/API ranking backend.

## Scope and limitations

- Protein and GEO coverage remain sparse and uneven across genes.
- V5 evaluates early ranking within a frozen candidate union, not open-world
  recall across every human gene and cell line.
- V5 contains no verified negatives and cannot support a negative-sink claim.
- The V4 weight search is development-stage model selection, not completed
  optimisation of production weights.
- Penalty V1 verifies implementation and ranking behaviour, not biological
  suitability or the optimal value of the cap.
- Experimental use still requires assay-specific and laboratory validation.

See `RELEASE_NOTES.md` for the final application checks and
`docs/MERGE_NOTES.md` for component ownership and runtime boundaries.
