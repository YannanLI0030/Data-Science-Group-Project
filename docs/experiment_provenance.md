# Experiment provenance and reproduction boundaries

The committed manifests record the exact inputs, scripts and outputs used for
each frozen run. Result-table hashes were rechecked before publication and
still match their manifests.

Some scripts were later edited to shorten comments or replace personal default
paths with repository-relative discovery. Those maintenance changes alter a
file hash even when the numerical method is unchanged. The manifests were not
rewritten, because doing so would falsely describe a later file as the script
used for the original run.

For publication, user-specific absolute paths in manifests and saved notebook
output were replaced with repository-relative paths. Files from the separate
runtime and review workspaces use the logical prefixes `external/runtime_snapshot/`,
`external/data_snapshot/` and `external/review_workspace/`. This cleanup changes
only path labels. Recorded SHA-256 values, row counts, labels, configurations
and numerical results remain those of the frozen runs. Hashes that identify a
manifest itself therefore refer to the original run-time copy, before this
publication-only path normalisation.

| Experiment | Exact run-time version in Git history | Later change |
|---|---|---|
| Dynamic 100-gene ablation and V5 pool build | `6ebd934` | Comments only |
| V5 label freeze and evaluation | `5cf8aa8` | Repository-relative paths |
| Exclusion-penalty structural run | `4acd4ca` | Comments and external-path discovery |
| V4 development weight search | `0a66a2c` | Comments only |

The final team application uses the later portable production script from the
2026-09-11 release. Relative paths, non-interactive CLI handling and the
downstream Agent were added after the Penalty run; the `score_candidates()`
implementation used by the structural comparison was unchanged. The original
production scorer remains identified by its hash in the Penalty manifest.
The subsequent final archive translated API, UI and data-merger messages into
English. Those text changes alter file hashes but do not change the scoring
logic or test structure.

For example, the exact dynamic runner can be inspected with:

```bash
git show 6ebd934:scripts/dynamic_ablation_runner.py
```

The current scripts are the convenient versions for a new local run. Use a new
output directory where the runner provides one, and do not overwrite frozen
results. A full rerun requires the external integrated data recorded by hash in
the relevant manifest. The Penalty experiment additionally requires the
separate dynamic production scorer and gene caches. AWS is not required.

The large `dynamic_ablation_unlabelled_top10.csv` intermediate is intentionally
not committed. Its 116,545 rows and SHA-256 hash are recorded in the dynamic
manifest. The frozen V5 review file, internal audit, Gold file and evaluation
outputs derived from that stage are committed.

The Node scripts under `scripts/v5/` use `@oai/artifact-tool`, which was part of
the original analysis environment. Their tracked outputs and manifests are the
portable evidence record; they are not required to inspect the reported V5
metrics.
