# Merge notes and component ownership

## Canonical scientific backend

The following components come from the dynamic multi-omics implementation and
are authoritative in the merged system:

- structured gene/protein/exclusion/disease inputs;
- on-demand loading from raw datasets 1-14;
- identifier harmonisation in `src/data_loader.py` and `src/data_merger.py`;
- disease hard filtering, missing-value handling and gene cache;
- `score_candidates()`, recommendation levels and final ordering;
- evidence trace, supplementary context and CSV output.

## Teammate Agent contribution adapted into the final runtime

The output layer keeps the teammate Agent project's central separation of
concerns:

- narrow deterministic or hosted writer interface;
- content-addressed evidence IDs;
- cited claim schema;
- L1 citation and L2 numeric grounding;
- optional structural entity checking;
- rejection of unsupported claims;
- auditable run manifest and ranking-integrity record.

These features were adapted to operate on the canonical dynamic recommender's
output. The earlier Agent percentile ranker and snapshot ranking path are not
called, because the agreed project requirement makes the dynamic recommender the
sole ranking authority.

## Runtime boundary

`src/agentic/output_agent.py` has no import of the scoring module, data loader or
raw matrices. It is invoked only after `score_candidates()` and receives a list
of already ranked dictionaries. Before and after explanation, it hashes rank,
DepMap ID, cell-line name, final score, confidence and recommendation level. A
mismatch raises an exception while preserving the deterministic report and CSV.
