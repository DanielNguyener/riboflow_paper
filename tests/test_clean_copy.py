"""Clean-copy reproduction: the panels build from a distribution-only copy of the tree.

A pass in the working directory proves nothing about a fresh checkout. `results/` is a
gigabyte of regenerated output and is gitignored; if a panel silently depends on a file
that happens to be lying around there, it works here and fails for everyone else.

So: copy the tree minus everything `.gitignore` excludes, run the panels there, and see
what actually happens.

The two coverage panels need a coverage HDF5, which is not distributed (see the README). They are expected to fail -- with an actionable message, which is
asserted, not with a traceback. Supply one and all 20 panels and all five figures are
exercised:

    RIBOFLOW_PAPER_COVERAGE_H5=results/coverage/HeLa.shared_coverage.h5 \\
        python -m pytest tests/test_clean_copy.py -q

The last tests check that, from the clean copy, every published TIFF under
figures/published/ is reproduced byte for byte.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

#: Panels that need only files shipped in the repository.
SELF_CONTAINED = ["fig02A", "fig02B", "fig03C", "fig03D",
                  "fig04", "fig04A", "fig04B", "fig04C",
                  "fig05A", "fig05B", "fig05C", "fig05D",
                  "fig05A_plos", "fig05B_plos", "fig05C_plos", "fig05D_plos",
                  "fig06A", "fig06B"]
#: Figures reproducible from shipped tables alone; 3 needs the coverage HDF5.
SELF_CONTAINED_FIGURES = [2, 4, 5, 6]
#: Panels that need a coverage HDF5, which is a documented pipeline product.
NEEDS_COVERAGE = ["fig03A", "fig03B"]

COVERAGE_H5 = os.environ.get("RIBOFLOW_PAPER_COVERAGE_H5")


# ── ignore classification, by real git ───────────────────────────────────────
# This used to be a hand-written matcher with a "last matching rule wins" loop. That is
# not git's rule: git never descends into an excluded DIRECTORY, so a negation cannot
# re-admit a file beneath one, and the matcher happily reported the opposite. It said the
# checksum manifest shipped under a `results/` + negation .gitignore that in reality
# ignored it.
#
# So ask git. `git check-ignore` is a pure path matcher -- the paths need not exist -- and
# it is the same implementation that decides what a real clone contains.

def _ignored_paths(paths, gitignore=None):
    """The subset of `paths` that git would ignore, as a set of relative strings.

    Runs in a throwaway repository seeded with this repo's `.gitignore`. The real
    repository is never initialised, staged or otherwise touched.
    """
    import subprocess
    import tempfile

    paths = [str(p) for p in paths]
    if not paths:
        return set()
    text = (gitignore or (REPO / ".gitignore")).read_text()
    if gitignore is None:
        # `.git/info/exclude` is part of what a real clone excludes; it holds the
        # local-only entries that are not in the published `.gitignore`.
        local = REPO / ".git" / "info" / "exclude"
        if local.exists():
            text += "\n" + local.read_text()
    with tempfile.TemporaryDirectory() as scratch:
        root = Path(scratch)
        subprocess.run(["git", "init", "-q", str(root)], check=True,
                       capture_output=True)
        (root / ".gitignore").write_text(text)
        result = subprocess.run(
            ["git", "-C", str(root), "check-ignore", "--stdin"],
            input="\n".join(paths), text=True, capture_output=True)
        # exit 0: some ignored; 1: none ignored; anything else is a real failure
        if result.returncode not in (0, 1):
            raise AssertionError("git check-ignore failed: %s" % result.stderr)
        return {line.strip() for line in result.stdout.splitlines() if line.strip()}


def ships(relative, gitignore=None) -> bool:
    """True when git would include `relative` in a clone."""
    return str(relative) not in _ignored_paths([relative], gitignore)


SHIPS = ("code/make_panels.py",
         "code/make_figures.py",
         "code/te_route/normalization.R",
         "data/read_taxonomy/taxonomy/taxonomy_all.tsv",
         "data/ribo_seq_qc/genome/tables/readlen_window_qc.csv",
         "data/coverage/concordance/region_concordance_per_sample.tsv",
         "data/alignment_fate/gene_partition_route7.tsv",
         "data/alignment_fate/locus_LRRFIP1.npz",
         "data/annotation/orf_catalog.tsv",
         "data/ribo_rna/counts/ribo_counts_genome.csv",
         "data/te_route/tables/per_gene_delta.tsv",
         "data/te_route/housekeeping/Housekeeping_GenesHuman.csv",
         "supporting_information/S1_Table/samples.csv",
         "supporting_information/S1_Table/build_s1_table.py",
         "config/panel_manifest.yaml",
         "config/inputs.example.yaml",
         "benchmark/summarize_benchmarks.py",
         "figures/panel_references/fig02A_readlen_psite_selection.pdf",
         # The published figures ship: they are what the clean-copy test reproduces.
         "figures/published/Fig2.tif", "figures/published/Fig6.tif",
         "figures/published/Fig4_plos.pdf",
         "results/coverage/coverage_checksums.tsv",
         "tests/conftest.py", "tests/test_clean_copy.py",
         "pytest.ini", "requirements-dev.txt")

DOES_NOT_SHIP = ("figures/published/Figure5_assembled.pdf",
                 "figures/published/Fig2.png",
                 "results/coverage/HeLa.shared_coverage.h5",
                 "results/coverage/concordance/region_concordance_per_sample.tsv",
                 "results/panels/fig02A.pdf",
                 "results/te_route/tables/per_gene_delta.tsv",
                 "results/alignment_fate/HeLa.gene_read_partition_reads.tsv",
                 "results/ribo_seq_qc/genome/tables/readlen_window_qc.csv",
                 "results/.cache/annotation/coverage_annotation.pkl",
                 "_build/regen/x.tsv",
                 "code/__pycache__/x.pyc",
                 ".pytest_cache/v/cache/lastfailed",
                 ".DS_Store",
                 # NOTE: CLAUDE.md and .claude/ are NOT here. They are excluded through
                 # .git/info/exclude, which is local to a clone and never pushed, so the
                 # published .gitignore stays free of them. The test below checks that
                 # mechanism against the real repository instead.
                 "config/local.yaml",
                 # a stray HDF5 under data/ must NOT slip past *.h5
                 "data/example.h5",
                 "sample.bam", "sample.bam.bai", "sample.bam.csi", "sample.cram",
                 "x.pkl", "x.tar.gz")


def test_the_release_ships_exactly_what_it_should():
    """Every path that must be in a clone, and every path that must not."""
    ignored = _ignored_paths(SHIPS + DOES_NOT_SHIP)
    wrongly_ignored = [p for p in SHIPS if p in ignored]
    wrongly_shipped = [p for p in DOES_NOT_SHIP if p not in ignored]
    assert not wrongly_ignored, "must ship but is ignored: %s" % (wrongly_ignored,)
    assert not wrongly_shipped, "must NOT ship but would: %s" % (wrongly_shipped,)


def test_the_checksum_manifest_is_the_only_shipped_file_under_results():
    """The one versioned generated file.

    `results/` cannot simply be ignored with a negation for this file -- git would never
    look inside the excluded directory. The rule set opens each level by its contents.
    """
    manifest = "results/coverage/coverage_checksums.tsv"
    assert ships(manifest)

    generated = [str(p.relative_to(REPO)) for p in (REPO / "results").rglob("*")
                 if p.is_file()]
    if len(generated) <= 20:
        pytest.skip("no regenerated results/ tree to check against (fresh clone)")
    ignored = _ignored_paths(generated)
    leaked = sorted(set(generated) - ignored - {manifest})
    assert not leaked, "these generated files would ship: %s" % leaked[:10]


def test_a_negation_cannot_reopen_an_excluded_parent(tmp_path):
    """Calibration: the exact mistake the old matcher made must still be detectable.

    A `results/` + `!results/.../file` .gitignore looks correct under "last rule wins" and
    is wrong in git. If this ever passes, the test above has stopped meaning anything.
    """
    broken = tmp_path / "broken_gitignore"
    broken.write_text("results/\n!results/coverage/coverage_checksums.tsv\n")
    assert not ships("results/coverage/coverage_checksums.tsv", broken), \
        "git re-admitted a file under an excluded directory -- it must not"


# ── the clean copy ───────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def clean_copy(tmp_path_factory):
    destination = tmp_path_factory.mktemp("clean")
    candidates = [p for p in sorted(REPO.rglob("*")) if p.is_file()
                  and not any(part in {".git", "__pycache__", ".pytest_cache"}
                              for part in p.relative_to(REPO).parts)]
    relatives = [str(p.relative_to(REPO)) for p in candidates]
    ignored = _ignored_paths(relatives)          # one git call for the whole tree
    copied = 0
    for source, relative in zip(candidates, relatives):
        if relative in ignored:
            continue
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        copied += 1
    assert copied > 100, "the clean copy is implausibly small: %d files" % copied
    return destination


def test_the_clean_copy_excludes_generated_output(clean_copy):
    assert not (clean_copy / "results" / "panels").exists()
    assert not list((clean_copy / "results").glob("**/*.h5")) if \
        (clean_copy / "results").exists() else True
    assert not (clean_copy / "_build").exists()
    assert not (clean_copy / "CLAUDE.md").exists()


def test_the_clean_copy_contains_the_pipeline_and_its_data(clean_copy):
    for required in ("code/make_panels.py", "code/make_tables.py",
                     "code/coverage/build_shared_coverage.py",
                     "config/panel_manifest.yaml", "config/cohort_manifest.tsv",
                     "data/alignment_fate/gene_partition_route7.tsv",
                     "data/te_route/tables/per_gene_delta.tsv",
                     "figures/published/Fig5.tif",
                     "supporting_information/S1_Table/samples.csv",
                     "supporting_information/S1_Table/build_s1_table.py"):
        assert (clean_copy / required).exists(), required


def run_panels(clean_copy, panels):
    return subprocess.run(
        [sys.executable, "code/make_panels.py", *panels, "--force"],
        capture_output=True, text=True, cwd=str(clean_copy))


@pytest.fixture(scope="module")
def self_contained_run(clean_copy):
    return run_panels(clean_copy, SELF_CONTAINED)


def test_every_self_contained_panel_builds_in_a_clean_copy(self_contained_run):
    output = self_contained_run.stdout + self_contained_run.stderr
    assert self_contained_run.returncode == 0, output[-4000:]
    for panel in SELF_CONTAINED:
        assert re.search(r"^%s\s+OK" % panel, output, re.MULTILINE), \
            "%s did not build:\n%s" % (panel, output[-4000:])


def test_the_self_contained_panels_wrote_their_files(clean_copy, self_contained_run):
    assert self_contained_run.returncode == 0
    written = list((clean_copy / "results" / "panels").glob("*.pdf"))
    assert len(written) >= len(SELF_CONTAINED)


@pytest.mark.parametrize("panel", NEEDS_COVERAGE)
def test_a_coverage_panel_fails_with_an_actionable_message(clean_copy, panel):
    """Not with a traceback. The coverage HDF5 is a documented product, and a reader who
    does not have one should be told how to make one."""
    if COVERAGE_H5:
        pytest.skip("a coverage file was supplied; the success path is tested below")
    result = run_panels(clean_copy, [panel])
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "does not exist" in output
    assert "build_shared_coverage.py" in output, \
        "the failure must say how to produce the missing input:\n%s" % output[-3000:]
    assert "Traceback" not in output, "it should refuse, not crash:\n%s" % output[-3000:]


@pytest.mark.skipif(not COVERAGE_H5,
                    reason="set RIBOFLOW_PAPER_COVERAGE_H5 to exercise the coverage panels")
def test_all_twenty_panels_build_when_a_coverage_file_is_supplied(clean_copy):
    """The complete claim: given the one documented pipeline product that is too large to
    distribute, a clean copy reproduces every panel."""
    source = Path(COVERAGE_H5).resolve()
    assert source.exists(), source
    target = clean_copy / "results" / "coverage" / source.name
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        shutil.copy2(source, target)

    result = run_panels(clean_copy, SELF_CONTAINED + NEEDS_COVERAGE)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-4000:]
    assert "20/20 panels produced" in output, output[-2000:]


# ── the published figures, byte for byte ─────────────────────────────────────

def _tiff_pixels(path):
    import numpy as np
    from PIL import Image
    with Image.open(path) as image:
        return image.size, image.mode, tuple(round(v) for v in image.info.get("dpi", (0, 0))), \
            np.asarray(image.convert("RGB")).copy()


@pytest.fixture(scope="module")
def assembled(clean_copy, self_contained_run):
    """Every self-contained figure, assembled in the clean copy into a scratch directory."""
    assert self_contained_run.returncode == 0
    out = clean_copy / "_assembled"
    figures = list(SELF_CONTAINED_FIGURES)
    if COVERAGE_H5:
        source = Path(COVERAGE_H5).resolve()
        target = clean_copy / "results" / "coverage" / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            shutil.copy2(source, target)
        figures.append(3)
    argv = [sys.executable, "code/assemble_figures.py", "--check", "--output-dir", str(out)]
    for number in figures:
        argv += ["--figure", str(number)]
    result = subprocess.run(argv, capture_output=True, text=True, cwd=str(clean_copy))
    assert result.returncode == 0, (result.stdout + result.stderr)[-4000:]
    return out, figures


def test_every_self_contained_figure_passes_the_plos_check(assembled):
    out, figures = assembled
    for number in figures:
        assert (out / ("Fig%d.tif" % number)).exists()
        assert (out / ("Fig%d_plos.pdf" % number)).exists()


@pytest.mark.parametrize("number", SELF_CONTAINED_FIGURES + [3])
def test_the_published_tiff_is_reproduced_byte_for_byte(assembled, number):
    """The five TIFFs under figures/published/ are regenerated from shipped inputs and
    compared as files. A pixel compare follows only to make a failure legible."""
    out, figures = assembled
    if number not in figures:
        pytest.skip("Figure %d needs RIBOFLOW_PAPER_COVERAGE_H5" % number)
    published = REPO / "figures" / "published" / ("Fig%d.tif" % number)
    rebuilt = out / ("Fig%d.tif" % number)
    if rebuilt.read_bytes() == published.read_bytes():
        return
    import numpy as np
    (size_a, mode_a, dpi_a, a), (size_b, mode_b, dpi_b, b) = \
        _tiff_pixels(rebuilt), _tiff_pixels(published)
    assert (size_a, mode_a, dpi_a) == (size_b, mode_b, dpi_b), \
        "Fig%d: %s %s %s vs published %s %s %s" % (number, size_a, mode_a, dpi_a,
                                                  size_b, mode_b, dpi_b)
    differing = np.argwhere((a != b).any(axis=2))
    pytest.fail("Fig%d.tif: same geometry, %d pixel(s) differ (first at row %d, col %d)"
                % (number, len(differing), differing[0][0], differing[0][1]))


# ── no panel reaches outside the repository ──────────────────────────────────

def test_no_panel_input_is_an_absolute_path():
    import yaml
    document = yaml.safe_load((REPO / "config" / "panel_manifest.yaml").read_text())
    for panel in document["panels"]:
        for key, value in (panel.get("inputs") or {}).items():
            assert not str(value).startswith("/"), "%s: %s" % (panel["id"], key)


def test_the_working_brief_is_excluded_from_the_real_repository():
    """`CLAUDE.md` and `.claude/` must never be committed.

    They are excluded via `.git/info/exclude` rather than `.gitignore`: that file is local
    to a clone and never pushed, so the published `.gitignore` carries no trace of them.
    Because the mechanism is local, this asserts against the REAL repository rather than a
    temporary one seeded from `.gitignore`.
    """
    import subprocess

    if not (REPO / ".git").exists():
        pytest.skip("not a git repository yet")
    for path in ("CLAUDE.md", ".claude/settings.local.json"):
        if not (REPO / path).exists():
            continue
        result = subprocess.run(["git", "-C", str(REPO), "check-ignore", "-q", path])
        assert result.returncode == 0, "%s would be committed" % path
