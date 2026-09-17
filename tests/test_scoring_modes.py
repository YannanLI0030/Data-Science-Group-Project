from __future__ import annotations

import unittest

from dynamic_cellline_selector_gene_protein import score_candidates


class ScoringModeTests(unittest.TestCase):
    @staticmethod
    def rows():
        return [
            {
                "DepMap_ID": "ACH-A",
                "cellLine": "A",
                "rnaExpr": 100.0,
                "protExpr": 10.0,
                "exclusionExpr": 0.0,
                "nRna": 1,
                "nProt": 1,
                "nExclusion": 1,
                "hasDepMapRNA": True,
                "hasHpaRNA": False,
                "hasGeoRNA": False,
                "hasProteomics": True,
            },
            {
                "DepMap_ID": "ACH-B",
                "cellLine": "B",
                "rnaExpr": 0.0,
                "protExpr": 1.0,
                "exclusionExpr": 0.0,
                "nRna": 1,
                "nProt": 1,
                "nExclusion": 1,
                "hasDepMapRNA": True,
                "hasHpaRNA": False,
                "hasGeoRNA": False,
                "hasProteomics": True,
            },
        ]

    def test_protein_only_uses_085_protein_and_015_confidence(self):
        scored = score_candidates(self.rows(), "PROTEIN_ONLY")
        for row in scored:
            expected = max(
                0.0,
                min(
                    1.0,
                    0.85 * row["proteinScore"]
                    + 0.15 * row["confidenceScore"]
                    - row["exclusionPenalty"],
                ),
            )
            self.assertAlmostEqual(row["finalScore"], expected, places=4)
            self.assertEqual(row["biologicalScore"], row["proteinScore"])

    def test_protein_only_rejects_missing_protein(self):
        rows = self.rows()
        rows[0]["protExpr"] = None
        rows[0]["nProt"] = 0
        rows[0]["hasProteomics"] = False
        scored = score_candidates(rows, "PROTEIN_ONLY")
        self.assertEqual([x["DepMap_ID"] for x in scored], ["ACH-B"])


if __name__ == "__main__":
    unittest.main()
