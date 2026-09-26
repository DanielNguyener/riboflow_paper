"""`code/clustering/validate_cluster_pseudogene_counts.py`: consensus-record typing and counting."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code" / "clustering"))
import validate_cluster_pseudogene_counts as vp  # noqa: E402

# v34 pseudogene genes on chr1 +: a processed one at [100, 200), a transcribed_unprocessed one at
# [180, 400) (overlaps the first), a unitary one at [1000, 1100); nothing on chr1 -.
V34 = {("chr1", "+"): (np.array([100, 180, 1000]), np.array([200, 400, 1100]),
                       ["processed", "unprocessed", "other"])}
CONSENSUS = pd.DataFrame([
    ("chr1", "+", 110, 160, "p1", "ENSG1"),     # inside the processed gene only
    ("chr1", "+", 150, 390, "p2", "ENSG1"),     # overlaps both; the unprocessed one more
    ("chr1", "+", 1010, 1090, "p3", "ENSG2"),   # unitary -> other
    ("chr1", "-", 110, 160, "p4", "ENSG1"),     # right place, wrong strand -> unclassified
    ("chr1", "+", 600, 700, "p5", "ENSG3"),     # no overlap -> unclassified
    ("chr2", "+", 0, 10, "p6", "ENSG1"),        # chromosome without pseudogenes -> unclassified
], columns=["chrom", "strand", "start", "end", "pseudogene_id", "parent_id"])
CLUSTERS = pd.DataFrame({"gene": ["A", "B", "C"], "gene_id": ["ENSG1.7", "ENSG2.1", "ENSG9.2"],
                         "transcript_id": ["tA", "tB", "tC"], "cluster": [2, 1, 1]})


def test_category_tests_unprocessed_before_processed():
    assert vp.category("transcribed_unprocessed_pseudogene") == "unprocessed"
    assert vp.category("translated_processed_pseudogene") == "processed"
    assert vp.category("unitary_pseudogene") == "other"


def test_classify_by_best_same_strand_overlap():
    assert vp.classify(CONSENSUS, V34) == ["processed", "unprocessed", "other",
                                           "unclassified", "unclassified", "unclassified"]


def test_annotate_counts_per_parent_and_zero_for_absent():
    consensus = CONSENSUS.assign(category=vp.classify(CONSENSUS, V34))
    t = vp.annotate(CLUSTERS, consensus).set_index("gene")
    assert t.loc["A", ["n_processed", "n_unprocessed", "n_other", "n_unclassified", "n_total"]].tolist() \
        == [1, 1, 0, 2, 4]
    assert t.loc["B", ["n_other", "n_total"]].tolist() == [1, 1]
    assert not t.loc["B", "has_processed"] and t.loc["B", "has_any"]
    assert t.loc["C", "n_total"] == 0 and not t.loc["C", "has_any"]
