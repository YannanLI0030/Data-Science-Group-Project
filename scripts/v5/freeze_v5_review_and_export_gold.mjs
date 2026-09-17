import crypto from "node:crypto";
import { execFile as execFileCallback } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { promisify } from "node:util";
import { fileURLToPath } from "node:url";
import { Workbook } from "@oai/artifact-tool";

const scriptPath = fileURLToPath(import.meta.url);
const baselineRepo = path.resolve(path.dirname(scriptPath), "../..");
const workingPath = process.env.V5_REVIEW_FILE ?? path.join(
  baselineRepo,
  "benchmarks/candidate_pool_v5_holdout_review_completed.csv",
);
const baselineCommit = "6ebd934";
const baselineRepoPath = "benchmarks/candidate_pool_v5_holdout_review.csv";
const poolManifestPath = path.join(
  baselineRepo,
  "benchmarks/candidate_pool_v5_holdout_manifest.json",
);
const outputDir = process.env.V5_OUTPUT_DIR ?? path.join(baselineRepo, "benchmarks");
const completedPath = path.join(outputDir, "candidate_pool_v5_holdout_review_completed.csv");
const goldPath = path.join(outputDir, "gold_standard_v5_holdout.csv");
const freezeManifestPath = path.join(outputDir, "gold_standard_v5_holdout_manifest.json");
const completedPreviewPath = path.join(os.tmpdir(), "v5_review_completed_rows181_189.png");
const goldPreviewPath = path.join(os.tmpdir(), "gold_standard_v5_holdout_preview.png");

const reviewerHeaders = [
  "review_id", "benchmark_role", "review_priority", "priority_reason", "gene",
  "DepMap_ID", "cell_line", "lineage", "disease", "subtype", "query_mode",
  "query_scope", "judgement", "benchmark_task", "evidence_type", "source_url",
  "evidence_summary", "verified", "review_notes",
];
const protectedColumns = reviewerHeaders.slice(0, 12);
const requiredReviewFields = reviewerHeaders.slice(12);
const goldHeaders = [
  "gene", "expected_cell_line", "expected_depmap_id", "relation",
  "benchmark_task", "evidence_type", "source_url", "evidence_summary",
  "source_hint", "verified", "notes",
];
const execFile = promisify(execFileCallback);
const gitExecutable = process.env.GIT_EXECUTABLE ?? "git";

function sha256Bytes(bytes) {
  return crypto.createHash("sha256").update(bytes).digest("hex");
}

function csvEscape(value) {
  const text = value == null ? "" : String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function csvBytes(rows) {
  const body = rows.map(row => row.map(csvEscape).join(",")).join("\r\n") + "\r\n";
  return Buffer.from(`\uFEFF${body}`, "utf8");
}

async function loadCsv(filePath, sheetName) {
  const bytes = await fs.readFile(filePath);
  return loadCsvBytes(bytes, sheetName, filePath);
}

async function loadCsvBytes(bytes, sheetName, sourceLabel) {
  const workbook = await Workbook.fromCSV(bytes.toString("utf8"), { sheetName });
  const sheet = workbook.worksheets.getItem(sheetName);
  const values = sheet.getUsedRange(true).values.map(row =>
    row.map(value => value == null ? "" : String(value))
  );
  if (values.length === 0) throw new Error(`Empty CSV: ${sourceLabel}`);
  values[0][0] = values[0][0].replace(/^\uFEFF/, "");
  return { bytes, workbook, sheet, values };
}

function indexColumns(headers) {
  return new Map(headers.map((header, index) => [header, index]));
}

function requireExactHeaders(actual, expected, label) {
  if (JSON.stringify(actual) !== JSON.stringify(expected)) {
    throw new Error(`${label} headers differ from the frozen schema`);
  }
}

function relationCounts(rows, col) {
  const counts = { positive: 0, negative: 0, unknown: 0 };
  for (const row of rows) counts[row[col.get("judgement")]] += 1;
  return counts;
}

const working = await loadCsv(workingPath, "completed_review");
const { stdout: baselineStdout } = await execFile(
  gitExecutable,
  ["-C", baselineRepo, "show", `${baselineCommit}:${baselineRepoPath}`],
  { encoding: "buffer", maxBuffer: 10 * 1024 * 1024 },
);
const baseline = await loadCsvBytes(
  Buffer.from(baselineStdout),
  "blank_review_baseline",
  `${baselineCommit}:${baselineRepoPath}`,
);
const poolManifest = JSON.parse(await fs.readFile(poolManifestPath, "utf8"));

requireExactHeaders(working.values[0], reviewerHeaders, "Working review");
requireExactHeaders(baseline.values[0], reviewerHeaders, "Blank reviewer baseline");
const workRows = working.values.slice(1);
const baselineRows = baseline.values.slice(1);
const col = indexColumns(reviewerHeaders);

if (workRows.length !== 189) throw new Error(`Expected 189 review rows, found ${workRows.length}`);
if (baselineRows.length !== 189) throw new Error(`Expected 189 baseline rows, found ${baselineRows.length}`);
if (poolManifest.counts?.candidate_rows !== 189) throw new Error("Candidate-pool manifest no longer declares 189 rows");

const baselineHash = sha256Bytes(baseline.bytes);
const expectedBaselineHash = poolManifest.outputs?.review_pool?.sha256;
if (baselineHash !== expectedBaselineHash) {
  throw new Error(`Blank reviewer baseline SHA-256 mismatch: ${baselineHash} != ${expectedBaselineHash}`);
}

for (let i = 0; i < workRows.length; i += 1) {
  const workKey = workRows[i][col.get("review_id")];
  const baselineKey = baselineRows[i][col.get("review_id")];
  if (workKey !== baselineKey) throw new Error(`Review order/key changed at data row ${i + 1}`);
  for (const name of protectedColumns) {
    if (workRows[i][col.get(name)] !== baselineRows[i][col.get(name)]) {
      throw new Error(`Protected field ${name} changed for ${workKey}`);
    }
  }
}

const seenReviewIds = new Set();
const seenGeneDepMap = new Set();
for (const row of workRows) {
  const reviewId = row[col.get("review_id")];
  const gene = row[col.get("gene")];
  const depmapId = row[col.get("DepMap_ID")];
  const judgement = row[col.get("judgement")].trim().toLowerCase();
  const verified = row[col.get("verified")].trim().toLowerCase();
  const benchmarkTask = row[col.get("benchmark_task")].trim().toLowerCase();
  if (seenReviewIds.has(reviewId)) throw new Error(`Duplicate review_id: ${reviewId}`);
  seenReviewIds.add(reviewId);
  const pair = `${gene}|${depmapId}`;
  if (seenGeneDepMap.has(pair)) throw new Error(`Duplicate gene/DepMap ID: ${pair}`);
  seenGeneDepMap.add(pair);
  if (!/^ACH-\d{6}$/.test(depmapId)) throw new Error(`Invalid DepMap ID for ${reviewId}: ${depmapId}`);
  if (!new Set(["positive", "negative", "unknown"]).has(judgement)) {
    throw new Error(`Invalid judgement for ${reviewId}: ${judgement}`);
  }
  if (!new Set(["expression_suitability", "both"]).has(benchmarkTask)) {
    throw new Error(`Invalid benchmark task for ${reviewId}: ${benchmarkTask}`);
  }
  if ((judgement === "unknown" && verified !== "no") ||
      (judgement !== "unknown" && verified !== "yes")) {
    throw new Error(`Judgement/verified mismatch for ${reviewId}`);
  }
  for (const name of requiredReviewFields) {
    if (!row[col.get(name)].trim()) throw new Error(`Blank ${name} for ${reviewId}`);
  }
}

const counts = relationCounts(workRows, col);
if (counts.positive !== 99 || counts.negative !== 0 || counts.unknown !== 90) {
  throw new Error(`Unexpected label counts: ${JSON.stringify(counts)}`);
}

const observedGeneCounts = Object.fromEntries(
  [...new Set(workRows.map(row => row[col.get("gene")]))].sort().map(gene => [
    gene,
    workRows.filter(row => row[col.get("gene")] === gene).length,
  ])
);
const expectedGeneCounts = Object.fromEntries(
  poolManifest.counts.per_gene.map(entry => [entry.gene, entry.candidate_rows])
);
if (JSON.stringify(observedGeneCounts) !== JSON.stringify(
  Object.fromEntries(Object.entries(expectedGeneCounts).sort(([a], [b]) => a.localeCompare(b)))
)) {
  throw new Error("Per-gene candidate counts differ from the frozen candidate-pool manifest");
}

working.sheet.getRangeByIndexes(0, 0, working.values.length, reviewerHeaders.length).values = working.values;
working.workbook.recalculate();
const completedBytes = csvBytes(working.values);

const verifiedRows = workRows.filter(row =>
  row[col.get("verified")].trim().toLowerCase() === "yes" &&
  ["positive", "negative"].includes(row[col.get("judgement")].trim().toLowerCase())
);
const goldRows = verifiedRows.map(row => {
  const evidenceSummary = row[col.get("evidence_summary")].trim();
  const sourceUrl = row[col.get("source_url")].trim();
  return [
    row[col.get("gene")].trim().toUpperCase(),
    row[col.get("cell_line")].trim(),
    row[col.get("DepMap_ID")].trim().toUpperCase(),
    row[col.get("judgement")].trim().toLowerCase(),
    row[col.get("benchmark_task")].trim().toLowerCase(),
    row[col.get("evidence_type")].trim(),
    sourceUrl,
    evidenceSummary,
    `${evidenceSummary}; Sources: ${sourceUrl}`,
    "yes",
    row[col.get("review_notes")].trim(),
  ];
}).sort((a, b) =>
  a[0].localeCompare(b[0]) || a[3].localeCompare(b[3]) || a[2].localeCompare(b[2])
);

if (goldRows.length !== 99) throw new Error(`Expected 99 gold rows, found ${goldRows.length}`);
if (new Set(goldRows.map(row => `${row[0]}|${row[2]}`)).size !== goldRows.length) {
  throw new Error("Gold export contains duplicate gene/DepMap IDs");
}
const goldWorkbook = Workbook.create();
const goldSheet = goldWorkbook.worksheets.add("gold_standard_v5_holdout");
goldSheet.getRangeByIndexes(0, 0, goldRows.length + 1, goldHeaders.length).values = [goldHeaders, ...goldRows];
goldWorkbook.recalculate();
const goldValues = goldSheet.getUsedRange(true).values.map(row =>
  row.map(value => value == null ? "" : String(value))
);
const goldBytes = csvBytes(goldValues);

await fs.mkdir(outputDir, { recursive: true });
for (const outputPath of [completedPath, goldPath, freezeManifestPath]) {
  try {
    await fs.access(outputPath);
    throw new Error(`Refusing to overwrite frozen output: ${outputPath}`);
  } catch (error) {
    if (error?.code !== "ENOENT") throw error;
  }
}
await fs.writeFile(completedPath, completedBytes);
await fs.writeFile(goldPath, goldBytes);

const positiveGenes = [...new Set(goldRows.filter(row => row[3] === "positive").map(row => row[0]))].sort();
const goldPerGene = Object.fromEntries(positiveGenes.map(gene => [
  gene,
  {
    positive: goldRows.filter(row => row[0] === gene && row[3] === "positive").length,
    negative: goldRows.filter(row => row[0] === gene && row[3] === "negative").length,
  },
]));
const freezeManifest = {
  schema_version: 1,
  created_utc: new Date().toISOString(),
  status: "labels_frozen_before_internal_audit_unblinding",
  experiment_role: "targeted independent-label component holdout",
  source_working_review: {
    path: workingPath,
    rows: workRows.length,
    sha256: sha256Bytes(working.bytes),
  },
  frozen_blank_baseline: {
    repository: baselineRepo,
    commit: baselineCommit,
    repository_path: baselineRepoPath,
    sha256: baselineHash,
    expected_sha256: expectedBaselineHash,
    protected_columns_verified_unchanged: protectedColumns,
  },
  checks: {
    review_rows: workRows.length,
    review_columns: reviewerHeaders.length,
    unique_review_ids: seenReviewIds.size,
    unique_gene_depmap_pairs: seenGeneDepMap.size,
    genes: Object.keys(observedGeneCounts).length,
    relation_counts: counts,
    missing_required_review_fields: 0,
    invalid_depmap_ids: 0,
    judgement_verified_mismatches: 0,
    internal_audit_opened_before_label_freeze: false,
  },
  outputs: {
    completed_review: {
      path: completedPath,
      rows: workRows.length,
      sha256: sha256Bytes(completedBytes),
    },
    gold_standard: {
      path: goldPath,
      rows: goldRows.length,
      genes: positiveGenes.length,
      relation_counts: { positive: goldRows.length, negative: 0 },
      per_gene: goldPerGene,
      sha256: sha256Bytes(goldBytes),
    },
  },
  interpretation_limits: [
    "Unknown rows are excluded from gold and are not treated as negatives.",
    "No independently verified negative labels were found; negative-sink metrics are unavailable for V5.",
    "Frozen reviewer labels must not be edited after internal rankings or configuration provenance are unblinded.",
  ],
};
await fs.writeFile(freezeManifestPath, `${JSON.stringify(freezeManifest, null, 2)}\n`, "utf8");

const completedCheck = await working.workbook.inspect({
  kind: "table",
  range: "completed_review!A182:S190",
  include: "values,formulas",
  tableMaxRows: 10,
  tableMaxCols: 19,
  tableMaxCellChars: 300,
  maxChars: 30000,
});
const goldCheck = await goldWorkbook.inspect({
  kind: "table",
  range: `gold_standard_v5_holdout!A1:K${goldRows.length + 1}`,
  include: "values,formulas",
  tableMaxRows: 12,
  tableMaxCols: 11,
  tableMaxCellChars: 200,
  maxChars: 30000,
});
const completedErrors = await working.workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 300 },
  summary: "completed-review formula error scan",
});
const goldErrors = await goldWorkbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 300 },
  summary: "gold formula error scan",
});
console.log(completedCheck.ndjson);
console.log(goldCheck.ndjson);
console.log(completedErrors.ndjson);
console.log(goldErrors.ndjson);

const completedPreview = await working.workbook.render({
  sheetName: "completed_review",
  range: "A182:S190",
  scale: 1,
  format: "png",
});
await fs.writeFile(completedPreviewPath, new Uint8Array(await completedPreview.arrayBuffer()));
const goldPreview = await goldWorkbook.render({
  sheetName: "gold_standard_v5_holdout",
  range: "A1:K14",
  scale: 1,
  format: "png",
});
await fs.writeFile(goldPreviewPath, new Uint8Array(await goldPreview.arrayBuffer()));

console.log(JSON.stringify({
  completedPath,
  goldPath,
  freezeManifestPath,
  completedPreviewPath,
  goldPreviewPath,
  counts,
  goldRows: goldRows.length,
  goldGenes: positiveGenes.length,
  goldPerGene,
  completedSha256: sha256Bytes(completedBytes),
  goldSha256: sha256Bytes(goldBytes),
}, null, 2));
