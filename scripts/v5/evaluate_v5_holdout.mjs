import crypto from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { Workbook } from "@oai/artifact-tool";

const root = "/Users/liyannan/Documents/Codex/2026-08-17/zhe";
const repo = "/Users/liyannan/Desktop/Data-Science-Group-Project";
const goldPath = path.join(root, "outputs/v5_holdout_review/gold_standard_v5_holdout.csv");
const completedReviewPath = path.join(root, "outputs/v5_holdout_review/candidate_pool_v5_holdout_review_completed.csv");
const labelFreezeManifestPath = path.join(root, "outputs/v5_holdout_review/gold_standard_v5_holdout_manifest.json");
const designPath = path.join(repo, "config/ablation_v5_frozen_design.json");
const poolManifestPath = path.join(repo, "benchmarks/candidate_pool_v5_holdout_manifest.json");
const auditPath = path.join(repo, "results/v5_holdout_internal/candidate_pool_v5_holdout_audit.csv");
const outputDir = path.join(root, "outputs/v5_holdout_evaluation");
const summaryPath = path.join(outputDir, "v5_holdout_summary.csv");
const detailPath = path.join(outputDir, "v5_holdout_by_gene_config.csv");
const contrastsPath = path.join(outputDir, "v5_holdout_controlled_contrasts.csv");
const tieAuditPath = path.join(outputDir, "v5_holdout_tie_audit.csv");
const evaluationManifestPath = path.join(outputDir, "v5_holdout_evaluation_manifest.json");
const previewDir = path.join(root, "work/v5_holdout_evaluation_previews");
const scriptPath = fileURLToPath(import.meta.url);

const ks = [3, 5, 10];
const baselineConfig = "A0_team_baseline";
const primaryMetric = "NDCG@5_verified_pool";
const bootstrapIterations = 5000;
const bootstrapSeed = 20260909;
const tolerance = 1e-12;

function sha256(bytes) {
  return crypto.createHash("sha256").update(bytes).digest("hex");
}

function stripBom(value) {
  return String(value ?? "").replace(/^\uFEFF/, "");
}

function csvEscape(value) {
  if (value == null || (typeof value === "number" && Number.isNaN(value))) return "";
  const text = String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function csvBytes(rows) {
  return Buffer.from(`\uFEFF${rows.map(row => row.map(csvEscape).join(",")).join("\r\n")}\r\n`, "utf8");
}

async function loadCsv(filePath, sheetName) {
  const bytes = await fs.readFile(filePath);
  const workbook = await Workbook.fromCSV(bytes.toString("utf8"), { sheetName });
  const sheet = workbook.worksheets.getItem(sheetName);
  const values = sheet.getUsedRange(true).values.map(row =>
    row.map(value => value == null ? "" : String(value))
  );
  if (!values.length) throw new Error(`Empty CSV: ${filePath}`);
  values[0][0] = stripBom(values[0][0]);
  const headers = values[0];
  const col = new Map(headers.map((header, index) => [header, index]));
  return { filePath, bytes, workbook, sheet, values, headers, col, rows: values.slice(1) };
}

function requireColumns(table, required) {
  const missing = required.filter(name => !table.col.has(name));
  if (missing.length) throw new Error(`${table.filePath} is missing columns: ${missing.join(", ")}`);
}

function asNumber(value, label) {
  const parsed = Number(value);
  if (!Number.isFinite(parsed)) throw new Error(`Expected finite number for ${label}, found ${value}`);
  return parsed;
}

function mean(values) {
  if (!values.length) return Number.NaN;
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}

function quantile(sortedValues, q) {
  if (!sortedValues.length) return Number.NaN;
  const index = (sortedValues.length - 1) * q;
  const lower = Math.floor(index);
  const upper = Math.ceil(index);
  if (lower === upper) return sortedValues[lower];
  const weight = index - lower;
  return sortedValues[lower] * (1 - weight) + sortedValues[upper] * weight;
}

function median(values) {
  return quantile([...values].sort((a, b) => a - b), 0.5);
}

function seedFor(label) {
  const prefix = crypto.createHash("sha256").update(label).digest("hex").slice(0, 8);
  return (bootstrapSeed ^ Number.parseInt(prefix, 16)) >>> 0;
}

function mulberry32(seed) {
  let state = seed >>> 0;
  return () => {
    state |= 0;
    state = (state + 0x6D2B79F5) | 0;
    let value = Math.imul(state ^ (state >>> 15), 1 | state);
    value = (value + Math.imul(value ^ (value >>> 7), 61 | value)) ^ value;
    return ((value ^ (value >>> 14)) >>> 0) / 4294967296;
  };
}

function pairedBootstrap(deltas, label) {
  if (!deltas.length) return { meanDelta: Number.NaN, ciLow: Number.NaN, ciHigh: Number.NaN };
  const observed = mean(deltas);
  if (deltas.every(value => Math.abs(value) <= tolerance)) {
    return { meanDelta: observed, ciLow: 0, ciHigh: 0 };
  }
  const random = mulberry32(seedFor(label));
  const samples = [];
  for (let iteration = 0; iteration < bootstrapIterations; iteration += 1) {
    let total = 0;
    for (let i = 0; i < deltas.length; i += 1) {
      total += deltas[Math.floor(random() * deltas.length)];
    }
    samples.push(total / deltas.length);
  }
  samples.sort((a, b) => a - b);
  return {
    meanDelta: observed,
    ciLow: quantile(samples, 0.025),
    ciHigh: quantile(samples, 0.975),
  };
}

function groupByExactScore(items) {
  const groups = [];
  for (const item of items) {
    const last = groups.at(-1);
    if (!last || last[0].score !== item.score) groups.push([item]);
    else last.push(item);
  }
  return groups;
}

function tieAwareAtK(groups, positiveIds, k) {
  let expectedHits = 0;
  let dcg = 0;
  let position = 1;
  for (const group of groups) {
    if (position > k) break;
    const groupSize = group.length;
    const groupPositive = group.filter(item => positiveIds.has(item.depmapId)).length;
    const slots = Math.min(groupSize, k - position + 1);
    const positiveFraction = groupPositive / groupSize;
    expectedHits += slots * positiveFraction;
    for (let rank = position; rank < position + slots; rank += 1) {
      dcg += positiveFraction / Math.log2(rank + 1);
    }
    position += groupSize;
  }
  let ideal = 0;
  for (let rank = 1; rank <= Math.min(k, positiveIds.size); rank += 1) {
    ideal += 1 / Math.log2(rank + 1);
  }
  return {
    expectedHits,
    precisionLowerBound: expectedHits / k,
    recallWithinUnion: expectedHits / positiveIds.size,
    ndcgObservedLabels: ideal ? dcg / ideal : Number.NaN,
  };
}

function firstPositiveTieMidrank(groups, positiveIds) {
  let position = 1;
  for (const group of groups) {
    const end = position + group.length - 1;
    if (group.some(item => positiveIds.has(item.depmapId))) return (position + end) / 2;
    position = end + 1;
  }
  return null;
}

function combination(n, k) {
  if (k < 0 || k > n) return 0;
  const reducedK = Math.min(k, n - k);
  let result = 1;
  for (let i = 1; i <= reducedK; i += 1) result = result * (n - reducedK + i) / i;
  return result;
}

function expectedReciprocalRankAtK(groups, positiveIds, k) {
  let before = 0;
  for (const group of groups) {
    if (before >= k) return 0;
    const groupSize = group.length;
    const positiveCount = group.filter(item => positiveIds.has(item.depmapId)).length;
    if (positiveCount === 0) {
      before += groupSize;
      continue;
    }
    const denominator = combination(groupSize, positiveCount);
    const maxJ = Math.min(groupSize - positiveCount + 1, k - before);
    let expected = 0;
    for (let j = 1; j <= maxJ; j += 1) {
      const probability = combination(groupSize - j, positiveCount - 1) / denominator;
      expected += probability / (before + j);
    }
    return expected;
  }
  return 0;
}

function largestTieSize(groups) {
  return Math.max(...groups.map(group => group.length));
}

function round(value, digits = 8) {
  if (!Number.isFinite(value)) return "";
  const factor = 10 ** digits;
  return Math.round((value + Number.EPSILON) * factor) / factor;
}

function countDirections(deltas) {
  let wins = 0;
  let ties = 0;
  let losses = 0;
  for (const delta of deltas) {
    if (delta > tolerance) wins += 1;
    else if (delta < -tolerance) losses += 1;
    else ties += 1;
  }
  return { wins, ties, losses };
}

function assertClose(actual, expected, label, epsilon = 1e-10) {
  if (Math.abs(actual - expected) > epsilon) {
    throw new Error(`${label}: expected ${expected}, observed ${actual}`);
  }
}

function selfTestMetrics() {
  const perfect = groupByExactScore([
    { depmapId: "P1", score: 1 },
    { depmapId: "N1", score: 0.5 },
    { depmapId: "N2", score: 0.4 },
  ]);
  const perfectMetrics = tieAwareAtK(perfect, new Set(["P1"]), 3);
  assertClose(perfectMetrics.precisionLowerBound, 1 / 3, "perfect P@3");
  assertClose(perfectMetrics.recallWithinUnion, 1, "perfect Recall@3");
  assertClose(perfectMetrics.ndcgObservedLabels, 1, "perfect NDCG@3");
  assertClose(1 / firstPositiveTieMidrank(perfect, new Set(["P1"])), 1, "perfect MRR");
  assertClose(expectedReciprocalRankAtK(perfect, new Set(["P1"]), 3), 1, "perfect expected MRR");

  const tied = groupByExactScore([
    { depmapId: "N0", score: 1 },
    { depmapId: "P1", score: 0.5 },
    { depmapId: "N1", score: 0.5 },
    { depmapId: "N2", score: 0.5 },
  ]);
  const tiedMetrics = tieAwareAtK(tied, new Set(["P1"]), 2);
  assertClose(tiedMetrics.expectedHits, 1 / 3, "partial tie expected hits");
  assertClose(tiedMetrics.precisionLowerBound, 1 / 6, "partial tie P@2");
  assertClose(tiedMetrics.recallWithinUnion, 1 / 3, "partial tie Recall@2");
  assertClose(firstPositiveTieMidrank(tied, new Set(["P1"])), 3, "tie midrank");
  assertClose(expectedReciprocalRankAtK(tied, new Set(["P1"]), 4), (1 / 2 + 1 / 3 + 1 / 4) / 3, "exact expected tied MRR");

  const noHit = groupByExactScore([
    { depmapId: "N1", score: 1 },
    { depmapId: "N2", score: 0.5 },
  ]);
  const noHitMetrics = tieAwareAtK(noHit, new Set(["P1"]), 2);
  assertClose(noHitMetrics.expectedHits, 0, "no-hit expected hits");
  if (firstPositiveTieMidrank(noHit, new Set(["P1"])) !== null) throw new Error("no-hit MRR should be censored to zero");
  assertClose(expectedReciprocalRankAtK(noHit, new Set(["P1"]), 2), 0, "no-hit expected MRR");
}

function rowsToWorkbook(headers, records, sheetName) {
  const workbook = Workbook.create();
  const sheet = workbook.worksheets.add(sheetName);
  const rows = [headers, ...records.map(record => headers.map(header => record[header] ?? ""))];
  sheet.getRangeByIndexes(0, 0, rows.length, headers.length).values = rows;
  workbook.recalculate();
  return { workbook, sheet, rows };
}

async function verifyWorkbook(workbook, sheetName, address, label) {
  const table = await workbook.inspect({
    kind: "table",
    range: `${sheetName}!${address}`,
    include: "values,formulas",
    tableMaxRows: 15,
    tableMaxCols: 30,
    tableMaxCellChars: 120,
    maxChars: 30000,
  });
  const errors = await workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
    options: { useRegex: true, maxResults: 300 },
    summary: `${label} formula error scan`,
  });
  console.log(table.ndjson);
  console.log(errors.ndjson);
}

selfTestMetrics();

const [gold, audit] = await Promise.all([
  loadCsv(goldPath, "gold_input"),
  loadCsv(auditPath, "audit_input"),
]);
const [labelFreezeManifest, poolManifest, design] = await Promise.all([
  fs.readFile(labelFreezeManifestPath, "utf8").then(JSON.parse),
  fs.readFile(poolManifestPath, "utf8").then(JSON.parse),
  fs.readFile(designPath, "utf8").then(JSON.parse),
]);

const inputHashes = {
  gold: sha256(gold.bytes),
  completedReview: sha256(await fs.readFile(completedReviewPath)),
  internalAudit: sha256(audit.bytes),
  frozenDesign: sha256(await fs.readFile(designPath)),
};
const expectedHashes = {
  gold: labelFreezeManifest.outputs.gold_standard.sha256,
  completedReview: labelFreezeManifest.outputs.completed_review.sha256,
  internalAudit: poolManifest.outputs.internal_audit.sha256,
  frozenDesign: poolManifest.design_sha256,
};
for (const name of Object.keys(expectedHashes)) {
  if (inputHashes[name] !== expectedHashes[name]) {
    throw new Error(`${name} SHA-256 mismatch: ${inputHashes[name]} != ${expectedHashes[name]}`);
  }
}
if (labelFreezeManifest.status !== "labels_frozen_before_internal_audit_unblinding") {
  throw new Error(`Unexpected label-freeze status: ${labelFreezeManifest.status}`);
}
if (design.status !== "frozen_before_v5_candidate_labelling") {
  throw new Error(`Unexpected frozen-design status: ${design.status}`);
}

const configs = design.primary_pool_configs;
if (!Array.isArray(configs) || configs.length !== 10 || !configs.includes(baselineConfig)) {
  throw new Error("Frozen design must contain the ten primary configs including A0");
}
const configDefinitions = new Map(design.config_definitions.map(item => [item.name, item]));
if (configs.some(config => !configDefinitions.has(config))) throw new Error("Missing primary config definition");

requireColumns(gold, ["gene", "expected_depmap_id", "relation", "verified"]);
requireColumns(audit, ["review_id", "gene", "DepMap_ID", "selection_stratum", ...configs.flatMap(config => [`rank_${config}`, `score_${config}`])]);

const goldPositiveByGene = new Map();
const goldPairs = new Set();
for (const row of gold.rows) {
  const gene = row[gold.col.get("gene")].trim().toUpperCase();
  const depmapId = row[gold.col.get("expected_depmap_id")].trim().toUpperCase();
  const relation = row[gold.col.get("relation")].trim().toLowerCase();
  const verified = row[gold.col.get("verified")].trim().toLowerCase();
  if (relation !== "positive" || verified !== "yes") throw new Error(`V5 gold contains non-positive or unverified row: ${gene}/${depmapId}`);
  if (!/^ACH-\d{6}$/.test(depmapId)) throw new Error(`Invalid gold DepMap ID: ${depmapId}`);
  const pair = `${gene}|${depmapId}`;
  if (goldPairs.has(pair)) throw new Error(`Duplicate gold pair: ${pair}`);
  goldPairs.add(pair);
  if (!goldPositiveByGene.has(gene)) goldPositiveByGene.set(gene, new Set());
  goldPositiveByGene.get(gene).add(depmapId);
}
if (gold.rows.length !== 99 || goldPositiveByGene.size !== 12) {
  throw new Error(`Expected 99 positives across 12 genes, found ${gold.rows.length} across ${goldPositiveByGene.size}`);
}

const auditPairs = new Set();
const auditRowsByGene = new Map();
for (const row of audit.rows) {
  const gene = row[audit.col.get("gene")].trim().toUpperCase();
  const depmapId = row[audit.col.get("DepMap_ID")].trim().toUpperCase();
  const pair = `${gene}|${depmapId}`;
  if (auditPairs.has(pair)) throw new Error(`Duplicate audit pair: ${pair}`);
  auditPairs.add(pair);
  if (!auditRowsByGene.has(gene)) auditRowsByGene.set(gene, []);
  auditRowsByGene.get(gene).push(row);
  for (const config of configs) {
    const rankText = row[audit.col.get(`rank_${config}`)].trim();
    const scoreText = row[audit.col.get(`score_${config}`)].trim();
    if ((rankText === "") !== (scoreText === "")) throw new Error(`Rank/score null mismatch for ${pair}/${config}`);
  }
}
if (audit.rows.length !== 189 || auditRowsByGene.size !== 12) {
  throw new Error(`Expected 189 audit rows across 12 genes, found ${audit.rows.length} across ${auditRowsByGene.size}`);
}
const missingGoldPairs = [...goldPairs].filter(pair => !auditPairs.has(pair));
if (missingGoldPairs.length) throw new Error(`Gold pairs absent from frozen audit: ${missingGoldPairs.slice(0, 10).join(", ")}`);
if ([...goldPositiveByGene.keys()].some(gene => !auditRowsByGene.has(gene))) throw new Error("Gold gene absent from audit");

const detailRecords = [];
const retrievalAudit = {};
for (const gene of [...goldPositiveByGene.keys()].sort()) {
  const positiveIds = goldPositiveByGene.get(gene);
  const geneAuditRows = auditRowsByGene.get(gene);
  const strata = new Set(geneAuditRows.map(row => row[audit.col.get("selection_stratum")].trim()));
  if (strata.size !== 1) throw new Error(`${gene} has inconsistent selection_stratum values`);
  const selectionStratum = [...strata][0];
  retrievalAudit[gene] = {};
  for (const config of configs) {
    const items = geneAuditRows
      .filter(row => row[audit.col.get(`score_${config}`)].trim() !== "")
      .map(row => ({
        depmapId: row[audit.col.get("DepMap_ID")].trim().toUpperCase(),
        score: asNumber(row[audit.col.get(`score_${config}`)], `${gene}/${config}/score`),
        storedRank: asNumber(row[audit.col.get(`rank_${config}`)], `${gene}/${config}/rank`),
      }))
      .sort((a, b) => a.storedRank - b.storedRank || a.depmapId.localeCompare(b.depmapId));
    if (items.length < 10) throw new Error(`${gene}/${config} has only ${items.length} Top-10 candidates`);
    for (let i = 0; i < items.length; i += 1) {
      if (items[i].storedRank !== i + 1) throw new Error(`${gene}/${config} stored ranks are not contiguous at ${items[i].depmapId}`);
      if (i > 0 && items[i].score > items[i - 1].score) throw new Error(`${gene}/${config} scores increase with rank`);
    }
    const groups = groupByExactScore(items);
    const cutoffScore = items.at(-1).score;
    const cutoffTieSize = groups.at(-1).length;
    const firstMidrank = firstPositiveTieMidrank(groups, positiveIds);
    const expectedMrr = expectedReciprocalRankAtK(groups, positiveIds, 10);
    const retrievedPositive = items.filter(item => positiveIds.has(item.depmapId)).length;
    const metrics = {};
    for (const k of ks) metrics[k] = tieAwareAtK(groups, positiveIds, k);
    const record = {
      gene,
      selection_stratum: selectionStratum,
      config,
      description: configDefinitions.get(config).description,
      n_candidates_top10_with_ties: items.length,
      n_gold_positive_within_union: positiveIds.size,
      n_positive_retrieved_top10_with_ties: retrievedPositive,
      n_unknown_retrieved_top10_with_ties: items.length - retrievedPositive,
      positive_coverage_top10_with_all_ties: round(retrievedPositive / positiveIds.size),
      labelled_positive_fraction_retrieved: round(retrievedPositive / items.length),
      first_positive_tie_midrank: firstMidrank ?? "",
      "MRR@10_expected_tie": round(expectedMrr),
      "MRR@10_midrank_legacy": round(firstMidrank == null ? 0 : 1 / firstMidrank),
      top_score: items[0].score,
      cutoff_score: cutoffScore,
      largest_exact_score_tie: largestTieSize(groups),
      cutoff_exact_score_tie: cutoffTieSize,
    };
    for (const k of ks) {
      record[`P@${k}_verified_lower_bound`] = round(metrics[k].precisionLowerBound);
      record[`Recall@${k}_verified_pool`] = round(metrics[k].recallWithinUnion);
      record[`NDCG@${k}_verified_pool`] = round(metrics[k].ndcgObservedLabels);
      record[`expected_positive_hits@${k}_tie_aware`] = round(metrics[k].expectedHits);
    }
    detailRecords.push(record);
    retrievalAudit[gene][config] = {
      n_candidates_top10_with_ties: items.length,
      cutoff_score: cutoffScore,
      cutoff_exact_score_tie: cutoffTieSize,
    };
  }
}

if (detailRecords.length !== 120) throw new Error(`Expected 120 gene/config records, found ${detailRecords.length}`);

const detailByConfigGene = new Map(detailRecords.map(record => [`${record.config}|${record.gene}`, record]));
const summaryRecords = [];
const contrastRecords = [];
const tieAuditRecords = detailRecords.map(record => ({
  gene: record.gene,
  selection_stratum: record.selection_stratum,
  config: record.config,
  n_candidates_top10_with_ties: record.n_candidates_top10_with_ties,
  top_score: record.top_score,
  cutoff_score: record.cutoff_score,
  largest_exact_score_tie: record.largest_exact_score_tie,
  cutoff_exact_score_tie: record.cutoff_exact_score_tie,
  has_any_exact_score_tie: record.largest_exact_score_tie > 1 ? "yes" : "no",
  has_cutoff_exact_score_tie: record.cutoff_exact_score_tie > 1 ? "yes" : "no",
}));
const tieChecks = {
  queries_with_10_candidates: tieAuditRecords.filter(record => record.n_candidates_top10_with_ties === 10).length,
  queries_with_11_candidates: tieAuditRecords.filter(record => record.n_candidates_top10_with_ties === 11).length,
  queries_with_any_exact_score_tie: tieAuditRecords.filter(record => record.has_any_exact_score_tie === "yes").length,
  queries_with_cutoff_exact_score_tie: tieAuditRecords.filter(record => record.has_cutoff_exact_score_tie === "yes").length,
};
if (
  tieChecks.queries_with_10_candidates !== 118 ||
  tieChecks.queries_with_11_candidates !== 2 ||
  tieChecks.queries_with_any_exact_score_tie !== 4 ||
  tieChecks.queries_with_cutoff_exact_score_tie !== 2
) {
  throw new Error(`Unexpected frozen tie structure: ${JSON.stringify(tieChecks)}`);
}
const metricNames = [
  ...ks.flatMap(k => [`P@${k}_verified_lower_bound`, `Recall@${k}_verified_pool`, `NDCG@${k}_verified_pool`]),
  "MRR@10_expected_tie",
  "MRR@10_midrank_legacy",
];
for (const config of configs) {
  const records = detailRecords.filter(record => record.config === config);
  const summary = {
    config,
    description: configDefinitions.get(config).description,
    aggregation: "macro mean across 12 frozen V5 genes",
    n_genes: records.length,
    n_gold_positive_within_union: records.reduce((sum, record) => sum + record.n_gold_positive_within_union, 0),
    n_positive_retrieved_top10_with_ties: records.reduce((sum, record) => sum + record.n_positive_retrieved_top10_with_ties, 0),
    micro_positive_coverage_top10_with_all_ties: round(
      records.reduce((sum, record) => sum + record.n_positive_retrieved_top10_with_ties, 0) /
      records.reduce((sum, record) => sum + record.n_gold_positive_within_union, 0)
    ),
    macro_positive_coverage_top10_with_all_ties: round(mean(records.map(record => record.positive_coverage_top10_with_all_ties))),
    "micro_Recall@10_verified_pool": round(
      records.reduce((sum, record) => sum + record["expected_positive_hits@10_tie_aware"], 0) /
      records.reduce((sum, record) => sum + record.n_gold_positive_within_union, 0)
    ),
    mean_candidates_top10_with_ties: round(mean(records.map(record => record.n_candidates_top10_with_ties))),
  };
  for (const metric of metricNames) summary[metric] = round(mean(records.map(record => record[metric])));

  for (const metric of [primaryMetric, "Recall@10_verified_pool", "MRR@10_expected_tie"]) {
    const values = records.map(record => record[metric]);
    const short = metric === primaryMetric ? "NDCG@5" : metric === "Recall@10_verified_pool" ? "Recall@10" : "MRR@10";
    const sorted = [...values].sort((a, b) => a - b);
    summary[`${short}_median`] = round(median(values));
    summary[`${short}_q1`] = round(quantile(sorted, 0.25));
    summary[`${short}_q3`] = round(quantile(sorted, 0.75));
  }

  for (const metric of [primaryMetric, "Recall@10_verified_pool", "MRR@10_expected_tie"]) {
    const deltas = records.map(record =>
      record[metric] - detailByConfigGene.get(`${baselineConfig}|${record.gene}`)[metric]
    );
    const bootstrap = pairedBootstrap(deltas, `${config}|${metric}`);
    const directions = countDirections(deltas);
    const metricLabel = metric === primaryMetric ? "NDCG@5" : metric === "Recall@10_verified_pool" ? "Recall@10" : "MRR@10";
    summary[`delta_${metricLabel}_vs_A0`] = round(bootstrap.meanDelta);
    summary[`delta_${metricLabel}_ci_low`] = round(bootstrap.ciLow);
    summary[`delta_${metricLabel}_ci_high`] = round(bootstrap.ciHigh);
    summary[`${metricLabel}_wins_vs_A0`] = directions.wins;
    summary[`${metricLabel}_ties_vs_A0`] = directions.ties;
    summary[`${metricLabel}_losses_vs_A0`] = directions.losses;
  }
  summaryRecords.push(summary);
}

const controlledContrasts = [
  { id: "C00_A0_vs_RNA_mean", question: "Multi-omic A0 versus RNA-only mean reference", reference: "B0_rna_mean", treatment: "A0_team_baseline" },
  { id: "C01_remove_direct_protein", question: "Remove direct Protein while retaining Protein-derived Confidence", reference: "A0_team_baseline", treatment: "A1a_no_direct_protein" },
  { id: "C02_remove_all_protein_evidence", question: "Remove direct Protein and all Protein-derived Confidence", reference: "A0_team_baseline", treatment: "A1b_no_protein_evidence" },
  { id: "C03_remove_protein_confidence_only", question: "Retain direct Protein but remove Protein-derived Confidence", reference: "A0_team_baseline", treatment: "A1c_no_protein_confidence" },
  { id: "C04_remove_protein_conf_after_direct", question: "Remove Protein-derived Confidence after direct Protein is already removed", reference: "A1a_no_direct_protein", treatment: "A1b_no_protein_evidence" },
  { id: "C05_remove_direct_without_protein_conf", question: "Remove direct Protein when Protein-derived Confidence is already absent", reference: "A1c_no_protein_confidence", treatment: "A1b_no_protein_evidence" },
  { id: "C06_remove_full_confidence", question: "Remove the complete Confidence block", reference: "A0_team_baseline", treatment: "A2_no_confidence" },
  { id: "C07_remove_completeness", question: "Remove Confidence completeness component", reference: "A0_team_baseline", treatment: "A2a_no_conf_completeness" },
  { id: "C08_remove_source_support", question: "Remove Confidence source-support component", reference: "A0_team_baseline", treatment: "A2b_no_conf_source_support" },
  { id: "C09_remove_consistency", question: "Remove Confidence RNA-Protein consistency component", reference: "A0_team_baseline", treatment: "A2c_no_conf_consistency" },
  { id: "C10_enable_adaptive_trust", question: "Enable correlation-based adaptive RNA trust", reference: "A0_team_baseline", treatment: "A4_adaptive_trust" },
];
for (const contrast of controlledContrasts) {
  for (const metric of [primaryMetric, "Recall@10_verified_pool", "MRR@10_expected_tie"]) {
    const referenceValues = [];
    const treatmentValues = [];
    const deltas = [];
    for (const gene of [...goldPositiveByGene.keys()].sort()) {
      const referenceValue = detailByConfigGene.get(`${contrast.reference}|${gene}`)[metric];
      const treatmentValue = detailByConfigGene.get(`${contrast.treatment}|${gene}`)[metric];
      referenceValues.push(referenceValue);
      treatmentValues.push(treatmentValue);
      deltas.push(treatmentValue - referenceValue);
    }
    const bootstrap = pairedBootstrap(deltas, `${contrast.id}|${metric}`);
    const directions = countDirections(deltas);
    contrastRecords.push({
      contrast_id: contrast.id,
      controlled_question: contrast.question,
      reference_config: contrast.reference,
      treatment_config: contrast.treatment,
      metric,
      delta_direction: "treatment_minus_reference",
      n_genes: deltas.length,
      reference_macro_mean: round(mean(referenceValues)),
      treatment_macro_mean: round(mean(treatmentValues)),
      mean_delta: round(bootstrap.meanDelta),
      median_delta: round(median(deltas)),
      paired_gene_bootstrap_ci_low: round(bootstrap.ciLow),
      paired_gene_bootstrap_ci_high: round(bootstrap.ciHigh),
      wins_for_treatment: directions.wins,
      ties: directions.ties,
      losses_for_treatment: directions.losses,
      bootstrap_iterations: bootstrapIterations,
      bootstrap_seed: seedFor(`${contrast.id}|${metric}`),
    });
  }
}

const detailHeaders = [
  "gene", "selection_stratum", "config", "description",
  "n_candidates_top10_with_ties", "n_gold_positive_within_union",
  "n_positive_retrieved_top10_with_ties", "n_unknown_retrieved_top10_with_ties",
  "positive_coverage_top10_with_all_ties", "labelled_positive_fraction_retrieved",
  "P@3_verified_lower_bound", "Recall@3_verified_pool", "NDCG@3_verified_pool", "expected_positive_hits@3_tie_aware",
  "P@5_verified_lower_bound", "Recall@5_verified_pool", "NDCG@5_verified_pool", "expected_positive_hits@5_tie_aware",
  "P@10_verified_lower_bound", "Recall@10_verified_pool", "NDCG@10_verified_pool", "expected_positive_hits@10_tie_aware",
  "first_positive_tie_midrank", "MRR@10_expected_tie", "MRR@10_midrank_legacy", "top_score", "cutoff_score",
  "largest_exact_score_tie", "cutoff_exact_score_tie",
];
const summaryHeaders = [
  "config", "description", "aggregation", "n_genes", "n_gold_positive_within_union",
  "n_positive_retrieved_top10_with_ties", "micro_positive_coverage_top10_with_all_ties",
  "macro_positive_coverage_top10_with_all_ties", "micro_Recall@10_verified_pool", "mean_candidates_top10_with_ties",
  "P@3_verified_lower_bound", "Recall@3_verified_pool", "NDCG@3_verified_pool",
  "P@5_verified_lower_bound", "Recall@5_verified_pool", "NDCG@5_verified_pool",
  "P@10_verified_lower_bound", "Recall@10_verified_pool", "NDCG@10_verified_pool",
  "MRR@10_expected_tie", "MRR@10_midrank_legacy",
  "NDCG@5_median", "NDCG@5_q1", "NDCG@5_q3",
  "Recall@10_median", "Recall@10_q1", "Recall@10_q3",
  "MRR@10_median", "MRR@10_q1", "MRR@10_q3",
  "delta_NDCG@5_vs_A0", "delta_NDCG@5_ci_low", "delta_NDCG@5_ci_high",
  "NDCG@5_wins_vs_A0", "NDCG@5_ties_vs_A0", "NDCG@5_losses_vs_A0",
  "delta_Recall@10_vs_A0", "delta_Recall@10_ci_low", "delta_Recall@10_ci_high",
  "Recall@10_wins_vs_A0", "Recall@10_ties_vs_A0", "Recall@10_losses_vs_A0",
  "delta_MRR@10_vs_A0", "delta_MRR@10_ci_low", "delta_MRR@10_ci_high",
  "MRR@10_wins_vs_A0", "MRR@10_ties_vs_A0", "MRR@10_losses_vs_A0",
];
const contrastHeaders = [
  "contrast_id", "controlled_question", "reference_config", "treatment_config", "metric",
  "delta_direction", "n_genes", "reference_macro_mean", "treatment_macro_mean", "mean_delta",
  "median_delta", "paired_gene_bootstrap_ci_low", "paired_gene_bootstrap_ci_high",
  "wins_for_treatment", "ties", "losses_for_treatment", "bootstrap_iterations", "bootstrap_seed",
];
const tieAuditHeaders = [
  "gene", "selection_stratum", "config", "n_candidates_top10_with_ties", "top_score", "cutoff_score",
  "largest_exact_score_tie", "cutoff_exact_score_tie", "has_any_exact_score_tie", "has_cutoff_exact_score_tie",
];

const summaryBook = rowsToWorkbook(summaryHeaders, summaryRecords, "v5_summary");
const detailBook = rowsToWorkbook(detailHeaders, detailRecords, "v5_by_gene_config");
const contrastBook = rowsToWorkbook(contrastHeaders, contrastRecords, "v5_controlled_contrasts");
const tieAuditBook = rowsToWorkbook(tieAuditHeaders, tieAuditRecords, "v5_tie_audit");

function columnName(count) {
  let number = count;
  let output = "";
  while (number > 0) {
    number -= 1;
    output = String.fromCharCode(65 + (number % 26)) + output;
    number = Math.floor(number / 26);
  }
  return output;
}
function previewRange(headers, records, maxDataRows = 14) {
  return `A1:${columnName(headers.length)}${Math.min(records.length, maxDataRows) + 1}`;
}

await verifyWorkbook(summaryBook.workbook, "v5_summary", previewRange(summaryHeaders, summaryRecords, 10), "summary");
await verifyWorkbook(detailBook.workbook, "v5_by_gene_config", previewRange(detailHeaders, detailRecords), "detail");
await verifyWorkbook(contrastBook.workbook, "v5_controlled_contrasts", previewRange(contrastHeaders, contrastRecords), "controlled contrasts");
await verifyWorkbook(tieAuditBook.workbook, "v5_tie_audit", previewRange(tieAuditHeaders, tieAuditRecords), "tie audit");

const outputSpecs = [
  { path: summaryPath, bytes: csvBytes(summaryBook.rows), workbook: summaryBook.workbook, sheet: "v5_summary", range: previewRange(summaryHeaders, summaryRecords, 10), preview: "summary.png", rows: summaryRecords.length, columns: summaryHeaders.length },
  { path: detailPath, bytes: csvBytes(detailBook.rows), workbook: detailBook.workbook, sheet: "v5_by_gene_config", range: previewRange(detailHeaders, detailRecords), preview: "detail.png", rows: detailRecords.length, columns: detailHeaders.length },
  { path: contrastsPath, bytes: csvBytes(contrastBook.rows), workbook: contrastBook.workbook, sheet: "v5_controlled_contrasts", range: previewRange(contrastHeaders, contrastRecords), preview: "controlled_contrasts.png", rows: contrastRecords.length, columns: contrastHeaders.length },
  { path: tieAuditPath, bytes: csvBytes(tieAuditBook.rows), workbook: tieAuditBook.workbook, sheet: "v5_tie_audit", range: previewRange(tieAuditHeaders, tieAuditRecords), preview: "tie_audit.png", rows: tieAuditRecords.length, columns: tieAuditHeaders.length },
];
for (const spec of outputSpecs) {
  try {
    await fs.access(spec.path);
    throw new Error(`Refusing to overwrite one-time V5 evaluation output: ${spec.path}`);
  } catch (error) {
    if (error?.code !== "ENOENT") throw error;
  }
}
try {
  await fs.access(evaluationManifestPath);
  throw new Error(`Refusing to overwrite one-time V5 evaluation manifest: ${evaluationManifestPath}`);
} catch (error) {
  if (error?.code !== "ENOENT") throw error;
}

await fs.mkdir(outputDir, { recursive: true });
await fs.mkdir(previewDir, { recursive: true });
for (const spec of outputSpecs) {
  await fs.writeFile(spec.path, spec.bytes);
  const preview = await spec.workbook.render({ sheetName: spec.sheet, range: spec.range, scale: 1, format: "png" });
  await fs.writeFile(path.join(previewDir, spec.preview), new Uint8Array(await preview.arrayBuffer()));
}

const evaluationManifest = {
  schema_version: 1,
  created_utc: new Date().toISOString(),
  status: "one_time_v5_holdout_evaluation_completed",
  experiment_role: "targeted independent-label component holdout",
  not_a_weight_search: true,
  frozen_primary_configs: configs,
  secondary_configs_excluded: design.secondary_unlabelled_sensitivity_configs ?? ["A3_equal_bio", "A3b_reversed_bio", "A5_v3_structure", "A6_v3_full"],
  baseline_config: baselineConfig,
  primary_metric: "macro mean NDCG@5_verified_pool across the 12 frozen genes",
  secondary_metrics: metricNames.filter(metric => metric !== primaryMetric),
  metric_scope: {
    candidate_scope: "Each configuration's frozen Top-10 with exact-score ties, evaluated against verified positives in the frozen ten-configuration candidate union.",
    query_scope: "GENE_MULTIOMICS queries only, with no disease filter. The 12 genes were prospectively selected by coverage stratum and structural discriminability.",
    unknown_policy: "Unknown candidates are unlabelled, not negatives. They occupy observed rank positions and receive zero observed gain; therefore P@k is explicitly a lower bound.",
    precision: "Expected verified-positive hits in the first k slots under random ordering within exact-score ties, divided by k; reported as P@k_verified_lower_bound.",
    recall: "Expected verified-positive hits in the first k slots divided by every verified positive for that gene in the frozen 189-row ten-configuration candidate union; reported as Recall@k_verified_pool.",
    ndcg: "Binary-gain NDCG using expected DCG under random ordering within exact-score ties and an IDCG containing all verified positives for the gene in the frozen pool; reported as NDCG@k_verified_pool.",
    mrr: "Exact expected reciprocal rank of the first verified positive under random ordering within its score-tie group, truncated at rank 10. MRR@10_midrank_legacy is retained only as an audit diagnostic.",
    negative_sink: "Unavailable because V5 contains no independently verified negative labels.",
  },
  aggregation: {
    primary_summary: "Unweighted macro mean, plus median and interquartile range, across the 12 frozen V5 genes.",
    paired_comparison: `Predeclared controlled gene-level deltas; summary also reports each configuration versus ${baselineConfig}.`,
    bootstrap: {
      iterations: bootstrapIterations,
      base_seed: bootstrapSeed,
      interval: "percentile 95% CI over genes",
      warning: "Intervals are descriptive and are not multiplicity-adjusted significance tests.",
    },
  },
  input_integrity: {
    all_sha256_checks_passed: true,
    gold: { path: goldPath, rows: gold.rows.length, genes: goldPositiveByGene.size, sha256: inputHashes.gold },
    completed_review: { path: completedReviewPath, sha256: inputHashes.completedReview },
    internal_audit: { path: auditPath, rows: audit.rows.length, sha256: inputHashes.internalAudit },
    frozen_design: { path: designPath, sha256: inputHashes.frozenDesign },
    label_freeze_manifest: { path: labelFreezeManifestPath, sha256: sha256(await fs.readFile(labelFreezeManifestPath)) },
    candidate_pool_manifest: { path: poolManifestPath, sha256: sha256(await fs.readFile(poolManifestPath)) },
    evaluation_script: { path: scriptPath, sha256: sha256(await fs.readFile(scriptPath)) },
  },
  checks: {
    metric_self_tests_passed: true,
    gold_pairs_all_present_in_audit: true,
    duplicate_gold_pairs: 0,
    duplicate_audit_pairs: 0,
    rank_score_null_mismatches: 0,
    config_gene_rows: detailRecords.length,
    each_config_gene_has_at_least_10_candidates: true,
    stored_top10_ranks_contiguous_and_scores_nonincreasing: true,
    candidate_pool_rows: audit.rows.length,
    verified_positive_rows: gold.rows.length,
    unknown_rows: audit.rows.length - gold.rows.length,
    global_verified_label_coverage: round(gold.rows.length / audit.rows.length),
    ...tieChecks,
    negatives: 0,
  },
  controlled_contrasts: controlledContrasts,
  outputs: Object.fromEntries(outputSpecs.map(spec => [path.basename(spec.path), {
    path: spec.path,
    rows: spec.rows,
    columns: spec.columns,
    sha256: sha256(spec.bytes),
  }])),
  interpretation_limits: [
    "The 12 genes were targeted for coverage strata and structural discriminability; they are not a random sample of all genes.",
    "The candidate universe is the frozen union of Top-10-with-ties from the ten primary configurations, so results measure comparative early-ranking performance within that union rather than open-world recall over all cell lines.",
    "V5 contains 99 positives and no verified negatives. It cannot evaluate exclusion behavior or negative sinking.",
    "Unknown labels make precision a lower bound and may change observed NDCG if later evidence is added.",
    "Most positive evidence is expression detection evidence, often independent proteomics; a positive label does not prove experimental superiority.",
    "Any configuration change motivated by these results converts V5 from holdout evidence into development evidence and requires a new untouched holdout for confirmation.",
  ],
  retrieval_audit: retrievalAudit,
};
await fs.writeFile(evaluationManifestPath, `${JSON.stringify(evaluationManifest, null, 2)}\n`, "utf8");

for (const spec of outputSpecs) {
  const savedBytes = await fs.readFile(spec.path);
  if (sha256(savedBytes) !== sha256(spec.bytes)) throw new Error(`Saved hash mismatch: ${spec.path}`);
  const saved = await loadCsv(spec.path, `saved_${spec.sheet}`);
  if (saved.rows.length !== spec.rows || saved.headers.length !== spec.columns) {
    throw new Error(`Saved dimensions mismatch for ${spec.path}`);
  }
  const savedErrors = await saved.workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
    options: { useRegex: true, maxResults: 300 },
    summary: `${path.basename(spec.path)} saved-file formula error scan`,
  });
  console.log(savedErrors.ndjson);
}

console.log(JSON.stringify({
  status: evaluationManifest.status,
  inputHashes,
  outputs: evaluationManifest.outputs,
  summary: summaryRecords,
}, null, 2));
