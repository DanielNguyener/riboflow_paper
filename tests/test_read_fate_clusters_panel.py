"""`code/panels/plot_read_fate_clusters.py`: the R merge <-> scipy linkage conversion and
the cut check that guards Figure 6."""
import sys
from pathlib import Path

import numpy as np
from scipy.cluster import hierarchy

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code" / "panels"))
import plot_read_fate_clusters as rf  # noqa: E402

# Two tight pairs far apart plus two loners: hclust would merge (1,2), (3,4), then the pairs,
# then the loners, in R's convention (negative = leaf, positive = earlier merge).
MERGE = np.array([[-1, -2], [-3, -4], [1, 2], [-5, -6], [3, 4]])
HEIGHT = np.array([0.1, 0.2, 1.0, 1.5, 5.0])


def test_round_trip_and_sizes():
    Z = rf.r_merge_to_linkage(MERGE, HEIGHT)
    assert Z.shape == (5, 4)
    assert Z[:, 3].tolist() == [2, 2, 4, 2, 6]           # cluster sizes
    assert Z[2, :2].tolist() == [6, 7]                   # merge 3 joins merges 1 and 2 (ids n + 0, n + 1)
    merge, height = rf.linkage_to_r_merge(Z)
    # scipy sorts the two ids, so the columns may swap; the unordered pair is what matters
    assert [sorted(r) for r in merge.tolist()] == [sorted(r) for r in MERGE.tolist()]
    assert np.allclose(height, HEIGHT)


def test_cut_matches_partition():
    Z = rf.r_merge_to_linkage(MERGE, HEIGHT)
    cut = hierarchy.fcluster(Z, 2, "maxclust")
    assert rf.same_partition(cut, [7, 7, 7, 7, 9, 9])
    assert rf.same_partition(cut, [2, 2, 2, 2, 1, 1])    # labels are irrelevant
    assert not rf.same_partition(cut, [1, 1, 2, 2, 3, 3])
    # the pairs joined (1.0) before the loners did (1.5), so k = 3 splits the loners
    assert rf.same_partition(hierarchy.fcluster(Z, 3, "maxclust"), [1, 1, 1, 1, 2, 3])
    assert rf.same_partition(hierarchy.fcluster(Z, 4, "maxclust"), [1, 1, 2, 2, 3, 4])


def test_cluster_colours_follow_the_dominant_fate():
    import pandas as pd
    centroids = pd.DataFrame({"cluster": [2, 1], "n": [1, 1],
                              **{"pct_" + c: [0.0, 0.0] for c in rf.COMPONENTS}})
    centroids.loc[0, "pct_genome_only_multimapped"] = 80.0   # cluster 2
    centroids.loc[1, "pct_shared_genome_unique"] = 90.0      # cluster 1
    assert rf.cluster_colours(centroids) == [rf.COLOURS[0], rf.COLOURS[3]]
