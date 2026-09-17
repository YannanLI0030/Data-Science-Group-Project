from __future__ import annotations

import copy
import unittest

from src.agentic.grounding import verify
from src.agentic.output_agent import OutputAgent, ranking_fingerprint
from src.agentic.schema import Claim, EvidenceCard


class AgentBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.rows = [
            {
                "rank": 1,
                "DepMap_ID": "ACH-000001",
                "cellLine": "HCS2",
                "finalScore": 0.91,
                "confidenceScore": 0.41,
                "recommendationLevel": "Strongly Recommended",
                "rnaScore": 1.0,
                "proteinScore": None,
                "biologicalScore": 1.0,
                "exclusionPenalty": 0.0048,
                "completenessScore": 0.5,
                "sourceSupportScore": 0.25,
                "hasDepMapRNA": True,
                "hasHpaRNA": False,
                "hasGeoRNA": False,
                "hasProteomics": False,
            },
            {
                "rank": 2,
                "DepMap_ID": "ACH-000002",
                "cellLine": "MS751",
                "finalScore": 0.88,
                "confidenceScore": 0.50,
                "recommendationLevel": "Strongly Recommended",
                "rnaScore": 0.95,
                "proteinScore": None,
                "biologicalScore": 0.95,
                "exclusionPenalty": 0.0,
                "completenessScore": 0.5,
                "sourceSupportScore": 0.25,
                "hasDepMapRNA": True,
                "hasHpaRNA": False,
                "hasGeoRNA": False,
                "hasProteomics": False,
            },
        ]
        self.request = {
            "target_gene": "EGFR", "target_protein": None,
            "exclusion_gene": "ABCB1", "disease": "cervical cancer",
            "query_mode": "GENE_MULTIOMICS", "top_n": 2,
        }

    def test_output_agent_cannot_change_ranking(self):
        before_rows = copy.deepcopy(self.rows)
        before_hash = ranking_fingerprint(self.rows)
        result = OutputAgent().run(
            self.request, self.rows,
            [{"dataset": "DepMap_RNAseq", "type": "PRIMARY_rna_expression",
              "value": 6.868, "unit": "log2(TPM+1)"}],
            [{"alternativeCellLine": "MS751", "similarityScore": 0.99}],
            {"global_signatures": {}, "top_mirna": [], "top_metabolites": []},
        )
        self.assertEqual(self.rows, before_rows)
        self.assertEqual(ranking_fingerprint(self.rows), before_hash)
        self.assertEqual(result.ranking_fingerprint, before_hash)

    def test_unknown_citation_is_dropped(self):
        card = EvidenceCard("known", "RANKING", "Score is 0.91.", "test", values=(0.91,))
        verdict = verify([Claim("Score is 0.91.", "SUPPORT", ("invented",))], [card])[0]
        self.assertFalse(verdict.accepted)
        self.assertEqual(verdict.failed_layer, "L1_citation")

    def test_unsupported_number_is_dropped(self):
        card = EvidenceCard("known", "RANKING", "Score is 0.91.", "test", values=(0.91,))
        verdict = verify([Claim("Score is 7.30.", "SUPPORT", ("known",))], [card])[0]
        self.assertFalse(verdict.accepted)
        self.assertEqual(verdict.failed_layer, "L2_numeric")


if __name__ == "__main__":
    unittest.main()
