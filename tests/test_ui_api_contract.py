from __future__ import annotations

import unittest

from api_server import _clean_query, _model_options, _weights


class UiApiContractTests(unittest.TestCase):
    def test_query_accepts_gene_protein_or_both_and_requires_disease(self):
        self.assertEqual(
            _clean_query({"target_gene": "egfr", "disease": "lung"})["target_gene"],
            "EGFR",
        )
        self.assertEqual(
            _clean_query({"target_protein": "HER2", "disease": "breast"})[
                "target_protein"
            ],
            "ERBB2",
        )
        with self.assertRaises(ValueError):
            _clean_query({"target_gene": "EGFR", "disease": ""})

    def test_protein_only_contract_exposes_085_and_zero_direct_rna(self):
        weights = _weights("PROTEIN_ONLY")
        self.assertEqual(weights["protein"], 0.85)
        self.assertEqual(weights["confidence"], 0.15)
        self.assertEqual(weights["rna"], 0.0)

    def test_custom_model_endpoint_is_validated(self):
        opts = _model_options(
            {
                "provider": "openai_compatible",
                "model": "local-model",
                "endpoint": "http://127.0.0.1:11434/v1/chat/completions",
            }
        )
        self.assertEqual(opts["provider"], "openai_compatible")
        with self.assertRaises(ValueError):
            _model_options(
                {"provider": "openai_compatible", "endpoint": "not-a-url"}
            )


if __name__ == "__main__":
    unittest.main()
