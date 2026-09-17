# Release notes - merged final package

These checks were run on the complete local release package. The raw datasets,
gene caches and generated run files are excluded from Git; their absence from
the repository does not change the recorded application results below.

## Verified workflow

The merged package was tested on the included raw datasets with:

```powershell
python start.py --cli --target_gene EGFR --exclusion_gene ABCB1 `
  --disease "cervical cancer" --top_n 3 --non-interactive
```

The deterministic ranking returned:

1. HCS2 - final score 0.9071
2. MS751 - final score 0.8793
3. CASKI - final score 0.8198

The Agent then emitted nine grounded explanatory claims and changed no ranking
field. A second run with `--no-agent` produced a byte-identical ranking CSV
(SHA-256 `ED3735D11AFBF823D0707F9F8730DE8D448E957B3E7A24FE6A6092A6DB647E3E`).

## Included safeguards

- the dynamic recommender is the sole ranking authority;
- Agent evidence is created only after scoring;
- L1 rejects unknown evidence IDs;
- L2 rejects numbers absent from cited evidence;
- full mode also rejects cross-candidate entity swaps;
- ranking fields are hashed before and after explanation;
- an Agent/API failure cannot suppress the original report or CSV.

## Expected data notice

The supplementary loader reports four miRNA columns that have no exact
`CCLE_Name` match in File 9. They are skipped and reported explicitly; this is
the behaviour of the supplied original loader, not an Agent error.

## Complete UI release (2026-09-10)

- Added a responsive three-step UI under `web/index.html`.
- Replaced the teammate snapshot ranker API with the canonical dynamic
  recommender in `api_server.py`.
- Added immutable `run_id` records: Agent explanations can only consume an
  existing saved ranking.
- Fixed candidate details so every row opens its own evidence chain instead of
  a static demonstration cell line.
- Added front-end model settings for scripted, OpenAI, Anthropic and custom
  OpenAI-compatible endpoints. UI-entered API keys are memory-only and are not
  persisted.
- Added an explicit protein-only scoring test. A live `EGFR` protein / `lung
  cancer` run returned NCIH1573 first and displayed the verified formula
  `0.85 × 1.0000 + 0.15 × 0.8675 = 0.9801`; RNA was labelled supporting-only.
- Re-ran the original `EGFR + exclude ABCB1 + cervical cancer` workflow through
  the UI. The first three remained HCS2 (0.9071), MS751 (0.8793) and CASKI
  (0.8198).
