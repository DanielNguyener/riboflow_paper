"""The standalone gate: this repository depends on nothing outside itself.

These tests are the reason the claim can be made rather than merely asserted. They fail if
someone reintroduces a path into the original project directory, a default that points at
another local checkout, or a call into a program that no longer exists.

Scope is the WHOLE of `code/`, not a favoured subset.
"""
from __future__ import annotations

import ast
import functools
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "code"

#: The original analysis directory and its part-numbered output tree. Nothing may name
#: either: no path, no filename, no prose mention.
ORIGINAL_PROJECT = re.compile(r"balanced_?25|fig1_outputs")

#: Author- or machine-specific absolute paths.
ABSOLUTE_PATH = re.compile(r"(?:^|[\s\"'=(])(/Users/|/home/|/corral|/scratch/|/work/)")

#: Other local repositories this code must not reach into.
OTHER_CHECKOUTS = ("RiboBase_Repositories", "riboflow_paper_audit_private")

#: The processing pipeline is RiboFlow_v2; its old repository name is not a citation.
DEPRECATED_PIPELINE_NAME = "riboflow_genome"

TEXT_SUFFIXES = {".py", ".md", ".yaml", ".yml", ".tsv", ".csv", ".json", ".cfg",
                 ".ini", ".txt", ".sh", ".bed", ".toml"}

#: Files that legitimately quote a path or the original project name: the S1-table
#: generator and its README document external inputs this repository does not distribute,
#: `numeric_claims.tsv` names the programs each published value came from, and this file
#: has to spell the forbidden patterns out in order to forbid them.
ALLOWED = {
    "supporting_information/S1_Table/build_s1_table.py",
    "supporting_information/S1_Table/samples.csv",
    "supporting_information/S1_Table/README.md",
    "docs/numeric_claims.tsv",
    "tests/test_standalone.py",
}


def tracked_files():
    """Every distributable text file: source, docs, configuration and data."""
    out = []
    for path in REPO.rglob("*"):
        if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
            continue
        relative = path.relative_to(REPO)
        parts = set(relative.parts)
        if parts & {"__pycache__", ".pytest_cache", "results", "_build", ".git",
                    "archive", ".claude", "logs"}:
            continue
        if str(relative) in ALLOWED or relative.name == "CLAUDE.md":
            continue
        if str(relative) == "config/local.yaml":
            # machine-local input paths; gitignored (test_clean_copy checks that), so
            # never distributed -- it is *supposed* to hold an absolute path.
            continue
        out.append(relative)
    return sorted(out)


def source_files():
    return sorted(p.relative_to(REPO) for p in CODE.rglob("*.py")
                  if "__pycache__" not in p.parts)


# ── nothing points at the original project ───────────────────────────────────

def test_there_are_files_to_check():
    """A vacuous pass would be worse than a failure."""
    assert len(tracked_files()) > 50
    assert len(source_files()) > 40


@functools.lru_cache(maxsize=None)
def _lines(relative):
    return (REPO / relative).read_text(errors="replace").splitlines()


#: (rule, files scanned, forbidden pattern). Source-only rules skip ALLOWED themselves.
FORBIDDEN = [
    ("the original project", tracked_files, ORIGINAL_PROJECT),
    ("an author-specific absolute path", source_files, ABSOLUTE_PATH),
    ("the deprecated pipeline repository (say RiboFlow_v2)", tracked_files,
     re.compile(re.escape(DEPRECATED_PIPELINE_NAME))),
    ("another local checkout", tracked_files,
     re.compile("|".join(re.escape(m) for m in OTHER_CHECKOUTS))),
    # Output directories are named for what they hold, not for a figure part.
    ("a part-numbered output directory", source_files,
     re.compile(r'"(0[0-9]|1[0-9])_[a-z_]+"')),
]


@pytest.mark.parametrize("rule,files,pattern", FORBIDDEN, ids=[f[0] for f in FORBIDDEN])
def test_no_file_names(rule, files, pattern):
    offenders = []
    for relative in files():
        if str(relative) in ALLOWED:
            continue
        for number, line in enumerate(_lines(relative), 1):
            if pattern.search(line):
                offenders.append("%s:%d: %s" % (relative, number, line.strip()[:100]))
    assert not offenders, "%s is named by:\n%s" % (rule, "\n".join(offenders))


def test_no_directory_is_named_after_a_figure():
    for directory in CODE.iterdir():
        if directory.is_dir():
            assert not re.match(r"fig\d", directory.name), directory.name


def test_the_suite_guards_the_repository_output_tree():
    """`conftest.py` must carry the session guard that fails if a test writes to
    `results/`. Asserted as source text because the guard's own effect is only visible at
    session teardown, which a test inside that session cannot observe."""
    text = (REPO / "tests" / "conftest.py").read_text()
    assert "repository_results_are_not_written_to" in text, \
        "the results/ isolation fixture is missing from conftest.py"
    assert 'scope="session", autouse=True' in text, \
        "the guard must be session-scoped and automatic, or a test can opt out of it"


def test_entry_point_checks_run_without_the_pipeline_environment():
    """`tests/test_cli.py` asks whether each program refuses with no arguments. That is a
    question about the code, so the subprocess must not inherit `RIBOFLOW_PAPER_*` from
    whoever is running pytest -- with the real GTF and APPRIS exported,
    `build_orf_catalog.py` had everything it needed and wrote into `results/`."""
    text = (REPO / "tests" / "test_cli.py").read_text()
    assert "bare_environment" in text and "RIBOFLOW_PAPER_" in text, \
        "test_cli.run() must strip RIBOFLOW_PAPER_* from the subprocess environment"


# ── the deleted implementations stay deleted ─────────────────────────────────

DELETED_PATHS = [
    ("code/make_examples.py", "(the worked-example export was removed)"),
    ("code/fig03_region_concordance", "code/coverage/"),
    ("code/fig04_ribo_rna_route", "code/ribo_rna/"),
    ("code/fig05_multimap_taxonomy", "code/read_taxonomy/"),
    ("code/fig02_readlen_psite", "code/panels/"),
    ("code/coverage/summarize_psite_placement.py",
     "(went with the reference-offset P-site rule)"),
    ("code/coverage/verify_against_published.py",
     "code/coverage/compute_coverage_concordance.py --compare"),
    ("code/gene_fate_atlas", "(no longer part of this repository)"),
    ("code/gene_discordance", "code/clustering/"),
    # The MAPQ >= 42 transcriptome-uniqueness split left Figures 4-6: transcriptome status
    # is presence in the post-dedup BAM (RiboFlow_v2 MAPQ >= 10), so the read-ID
    # concordance stage and its categories had nothing left to feed.
    ("code/read_taxonomy/compute_concordance.py", "(the taxonomy is present/absent now)"),
    ("code/read_taxonomy/concordance_lib.py", "code/read_taxonomy/reference_lib.py"),
    ("archive", "(no longer part of this repository)"),
    ("docs/audit", "docs/data_availability.md"),
    ("code/ribo_rna/compute_ribo_rna_route.py",
     "code/ribo_rna/count_transcript_reads.py (it re-derived the same universe and "
     "re-counted the same four BAMs)"),
    # The flattening: the two figure sub-projects and the per-figure assemblers they
    # replaced, plus the Figure 3/4E panels the manuscript no longer carries.
    ("TE_Estimation", "code/te_route/ + code/ribo_rna/build_count_matrices.py"),
    ("Alternative_Isoforms", "code/alignment_fate/ + code/panels/ + assemble_figures.py"),
    ("code/assemble_plos.py", "code/assemble_figures.py"),
    ("figures/assemble_figure3.py", "code/assemble_figures.py"),
    ("code/ribo_rna/run_ribo_rna_route.py", "code/ribo_rna/build_count_matrices.py"),
    ("code/alignment_fate/build_transcript_fates.py",
     "code/alignment_fate/build_gene_read_partition.py"),
    ("code/panels/plot_transcript_alignment_fates.py", "code/panels/plot_gene_partition.py"),
    ("code/panels/plot_ribo_rna_scatter.py", "code/te_route/plot_te_route_panels.py"),
    ("code/panels/plot_ribo_rna_route_summary.py", "code/te_route/plot_te_route_panels.py"),
    ("data/ribo_rna/ribo_rna_route_cds.tsv", "data/te_route/tables/route_correlation.tsv"),
    ("data/alignment_fate/HeLa.transcript_alignment_fates.tsv",
     "data/alignment_fate/gene_partition_route7.tsv"),
]


@pytest.mark.parametrize("path,replacement", DELETED_PATHS, ids=lambda v: str(v)[:48])
def test_deleted_paths_stay_gone(path, replacement):
    assert not (REPO / path).exists(), "%s should be gone; use %s" % (path, replacement)


def test_no_program_invokes_a_deleted_one():
    """A stage that shells out to something that no longer exists fails at runtime, not
    at import, so nothing else would catch it."""
    gone = ["cluster_gene_fate.py", "run_fate_cluster.py", "run_fate_category.py",
            "plot_fate_cluster.py", "run_gene_fate_atlas.py", "compute_gene_fate.py",
            "build_gene_covariates.py",
            "build_gene_discordance.py", "cluster_hellinger.R",
            "gene_discordance_report.py", "plot_classification.py", "null_f_alt.py",
            "dedup_dribo.py", "dedup_artifact_test.py", "plot_pseudo_vs_dribo.py",
            "make_examples.py", "analyze_transcript_pseudogene_tie.py",
            "run_region_concordance.py", "run_region_coverage.py",
            "plot_ribo_rna_counts_scatter.py", "plot_pooled_with_example.py",
            "plot_union_combined.py", "verify_audit_baseline.py",
            "summarize_psite_placement.py", "verify_against_published.py",
            "export_example_vectors.py", "prototype_offset_periodicity.py",
           "compute_ribo_rna_route.py",
           "assemble_plos.py", "assemble_figure3.py", "run_ribo_rna_route.py",
           "build_transcript_fates.py", "plot_transcript_alignment_fates.py",
           "plot_ribo_rna_scatter.py", "plot_ribo_rna_route_summary.py",
           "plot_fig5_panels.py",
           # The QC pipeline's own diagnostic plotters, deleted with the rest of the
           # non-panel figure output. `run_pipeline.py` went on shelling out to all four
           # and printing a warning when each failed, so a full run emitted five
           # "failed (exit 2)" lines and still exited 0.
           "plot_frame_periodicity.py", "plot_cds_frame_breakdown.py",
           "plot_region_breakdown.py", "plot_readlen_distribution.py",
           # The clustering exploration, reduced to the Figure 6 path: one R script
           # (ward_cluster.R) and the panel generator in code/panels/.
           "run_clustering.py", "plot_clustering.py", "plot_hclust.py",
           "hellinger_transform.R", "hclust_sweep.R", "choose_k.R", "kmeans_sweep.R",
           "ilr_transform.R", "clustering_lib.R",
           "compute_concordance.py", "concordance_lib.py"]
    offenders = []
    for relative in source_files():
        text = (REPO / relative).read_text()
        for name in gone:
            if name in text:
                offenders.append("%s mentions %s" % (relative, name))
    assert not offenders, "\n".join(offenders)


def test_the_default_workflow_calls_only_existing_programs():
    """Walk every `sys.executable` invocation in the launchers and assert the target is a
    real file."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("make_tables", CODE / "make_tables.py")
    make_tables = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(make_tables)

    referenced = set()
    for node in ast.walk(ast.parse((CODE / "make_tables.py").read_text())):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) \
                and node.value.endswith(".py"):
            referenced.add(node.value)
    for name in referenced:
        matches = list(CODE.rglob(name))
        assert matches, "make_tables.py names %s, which does not exist" % name


# ── cross-directory imports resolve ──────────────────────────────────────────
# Several modules reach a sibling package by putting its directory on sys.path rather than
# by importing it, because the drivers spawn subprocesses that import the same modules by
# bare name. A path built from a directory that no longer exists is then silent: the insert
# succeeds, and the ImportError only surfaces at run time, with a BAM already open.

def test_every_directory_put_on_sys_path_exists():
    pattern = re.compile(r'(?:REPO|_CODE|_HERE\.parent|parents\[\d\])'
                         r'((?:\s*/\s*"[A-Za-z0-9_]+")+)')
    offenders = []
    for relative in source_files():
        text = (REPO / relative).read_text()
        if "sys.path" not in text:
            continue
        for lineno, line in enumerate(text.splitlines(), 1):
            match = pattern.search(line)
            if not match:
                continue
            parts = re.findall(r'"([A-Za-z0-9_]+)"', match.group(1))
            if parts[0] not in ("code", "common", "read_taxonomy", "coverage",
                                "panels", "ribo_seq_qc", "alignment_fate", "ribo_rna",
                                "te_route", "clustering"):
                continue
            base = REPO if parts[0] == "code" else CODE
            candidate = base.joinpath(*parts)
            if not candidate.exists():
                offenders.append("%s:%d -> %s" % (relative, lineno, candidate))
    assert not offenders, "sys.path built from missing directories:\n" + "\n".join(offenders)


def test_the_alignment_fate_loader_resolves_its_siblings():
    """`load_libraries()` is the only cross-package import in the gene-partition stage, and
    no other test imports it -- a stale directory there breaks Figure 5A and nothing else
    notices."""
    sys.path.insert(0, str(CODE / "alignment_fate"))
    try:
        import gene_read_partition_lib
        modules = gene_read_partition_lib.load_libraries()
    finally:
        sys.path.remove(str(CODE / "alignment_fate"))
    for module in modules:
        assert Path(module.__file__).parent == CODE / "read_taxonomy", module.__file__


# ── configuration is passed, not inherited ───────────────────────────────────

def code_only(path):
    """Source with comments and string literals removed.

    Docstrings discuss environment variables -- that is how a reader learns
    what is configurable. Matching raw text would forbid explaining the thing.
    """
    import tokenize
    pieces = []
    with open(path, "rb") as handle:
        for token in tokenize.tokenize(handle.readline):
            if token.type in (tokenize.COMMENT, tokenize.STRING, tokenize.NL,
                              tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT):
                continue
            pieces.append((token.start[0], token.string))
    return pieces


@pytest.mark.parametrize("relative", source_files(), ids=lambda p: p.name)
def test_no_module_reads_the_environment_at_import_time(relative):
    """An `os.environ` read at module scope is what makes code impossible to configure
    in-process and forces launchers to smuggle values through PYTHONPATH. Reads INSIDE a
    function are fine: that is an explicit fallback, not import-time state."""
    tree = ast.parse((REPO / relative).read_text())
    offenders = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        for child in ast.walk(node):
            if (isinstance(child, ast.Attribute) and isinstance(child.value, ast.Name)
                    and child.value.id == "os" and child.attr == "environ"):
                offenders.append(child.lineno)
    assert not offenders, "%s reads os.environ at MODULE scope (line %s)" % (
        relative, offenders)


@pytest.mark.parametrize("relative", source_files(), ids=lambda p: p.name)
def test_no_module_preregisters_itself_in_sys_modules(relative):
    """Pre-registering a module makes which implementation runs depend on import order."""
    for node in ast.walk(ast.parse((REPO / relative).read_text())):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if (isinstance(target, ast.Subscript)
                    and isinstance(target.value, ast.Attribute)
                    and isinstance(target.value.value, ast.Name)
                    and target.value.value.id == "sys"
                    and target.value.attr == "modules"):
                pytest.fail("%s pre-registers a module in sys.modules at line %d"
                            % (relative, node.lineno))


def test_every_environment_variable_is_namespaced():
    """One prefix, so `env | grep RIBOFLOW_PAPER` shows the whole configuration."""
    pattern = re.compile(r'os\.environ(?:\.get|\.setdefault)?\(\s*"([A-Z_][A-Z0-9_]*)"')
    allowed = {"PYTHONPATH", "MPLBACKEND"}
    found = set()
    for relative in source_files():
        found.update(pattern.findall((REPO / relative).read_text()))
    stray = {name for name in found - allowed
             if not name.startswith("RIBOFLOW_PAPER_")}
    assert not stray, "un-namespaced environment variables: %s" % sorted(stray)


def test_no_environment_variable_defaults_to_a_path_outside_the_repository():
    """The failure this whole gate exists for: a default that silently resolves to a
    directory that only ever existed on one machine."""
    pattern = re.compile(
        r'os\.environ\.get\(\s*"[A-Z_]+"\s*,\s*(.+?)\)\s*$', re.MULTILINE)
    offenders = []
    for relative in source_files():
        for number, line in enumerate((REPO / relative).read_text().splitlines(), 1):
            match = pattern.search(line)
            if match and ORIGINAL_PROJECT.search(match.group(1)):
                offenders.append("%s:%d" % (relative, number))
    assert not offenders, "\n".join(offenders)


# ── the panel manifest ───────────────────────────────────────────────────────

def load_manifest():
    """The manifest as the pipeline sees it: figure keys normalised to strings."""
    sys.path.insert(0, str(REPO / "code"))
    import make_panels
    document, _panels = make_panels.load_manifest(REPO / "config" / "panel_manifest.yaml")
    return document


def figure_stem(number):
    sys.path.insert(0, str(REPO / "code"))
    import make_panels
    return make_panels.figure_stem(number)


def test_the_panel_manifest_declares_exact_output_paths():
    document = load_manifest()
    outputs = [p["output"] for p in document["panels"] if p.get("generator")]
    # S1A-B, 2A-D, 3 (+ its three single panels), 4A-D page-size renders, 5A-B, 6.
    assert len(outputs) == 17, "expected 17 panel outputs, found %d" % len(outputs)
    assert len(set(outputs)) == len(outputs), "duplicate output paths"
    for output in outputs:
        assert not ORIGINAL_PROJECT.search(output), output
        assert "*" not in output and "?" not in output, "%s looks like a glob" % output
        assert output.startswith("results/"), output


def test_every_panel_generator_exists():
    for panel in load_manifest()["panels"]:
        generator = panel.get("generator")
        if generator:
            assert (REPO / generator).exists(), "%s: %s" % (panel["id"], generator)


def test_every_declared_panel_input_is_shipped_or_a_coverage_product():
    """A panel input is either a file in the repository or the coverage HDF5.

    Those are the only two kinds. No panel reads a file another panel wrote: values several
    panels share -- Figure 3's axis maximum, Figure 4's cohort ordering -- are pure
    functions of a table they all declare, so each derives its own. Anything else would be
    a hidden dependency on whatever happens to be lying around in `results/`.
    """
    unexplained = []
    for panel in load_manifest()["panels"]:
        for key, value in (panel.get("inputs") or {}).items():
            if (REPO / value).exists():
                continue
            if value.startswith("results/coverage/"):
                continue          # the documented coverage product, built from BAMs
            unexplained.append("%s: %s -> %s" % (panel["id"], key, value))
    assert not unexplained, "\n".join(unexplained)


def test_no_panel_writes_a_file_another_panel_reads():
    """`results/panels/` holds the final PDFs and nothing else."""
    document = load_manifest()
    assert "derived_inputs" not in document, \
        "a derived-input stage makes panel generation order-dependent"
    for panel in document["panels"]:
        assert "emits" not in panel, "%s emits a side product" % panel["id"]
        for key, value in (panel.get("inputs") or {}).items():
            assert not str(value).startswith("results/panels/"), \
                "%s reads %s, which another panel would have to write first" % (
                    panel["id"], value)


#: Main-text figures by number, supporting figures as `S<n>`; keys are strings.
EXPECTED_FIGURES = {"S1", "2", "3", "4", "5", "6"}


def test_the_figures_block_names_every_published_figure_once():
    """Six figures (five main-text, one supporting), each with a known composer and
    raster rule, output under figures/published/, built from panels the manifest declares."""
    document = load_manifest()
    figures = document["figures"]
    assert set(figures) == EXPECTED_FIGURES, sorted(figures)
    ids = {p["id"] for p in document["panels"] if p.get("generator")}
    outputs = []
    for number, spec in figures.items():
        assert spec["composer"] in ("grid_pair", "fitted_rows", "single_panel", "rows_1to1"), number
        assert spec["raster"]["kind"] in ("matplotlib_pdf", "generator_tiff", "fitz"), number
        assert 300 <= int(spec["raster"]["dpi"]) <= 600, number
        assert spec["output"] == "figures/published/" + figure_stem(number), spec["output"]
        outputs.append(spec["output"])
        rows = spec.get("rows") or [spec["panels"]]
        used = [panel for row in rows for panel in row]
        assert used, "figure %s has no panels" % number
        unknown = [panel for panel in used if panel not in ids]
        assert not unknown, "figure %s uses undeclared panel(s) %s" % (number, unknown)
        if spec["composer"] == "rows_1to1":
            assert len(spec["letters"]) == len(used), number
            assert {"max_width_pt", "margin_pt", "gutter_pt", "row_gap_pt",
                    "letter_pt"} <= set(spec["page"]), number
    assert len(set(outputs)) == len(outputs), "two figures share an output stem"


def test_every_figure_panel_a_rows_figure_places_is_in_results_panels():
    """`rows_1to1` reads `results/panels/<output>.pdf`, so those panels' outputs must be
    manifest outputs -- the composer never renders on its own."""
    document = load_manifest()
    by_id = {p["id"]: p for p in document["panels"] if p.get("generator")}
    for number, spec in document["figures"].items():
        if spec["composer"] != "rows_1to1":
            continue
        for row in spec["rows"]:
            for panel in row:
                assert by_id[panel]["output"].startswith("results/panels/"), panel


def test_the_published_figures_ship():
    for number in sorted(EXPECTED_FIGURES):
        for suffix in (".tif", "_plos.pdf"):
            path = REPO / "figures" / "published" / (figure_stem(number) + suffix)
            assert path.exists(), path


def test_the_launcher_does_not_glob_for_its_outputs():
    """`sorted(BUILD.rglob(name))[0]` let a stale artifact shadow a fresh one."""
    tokens = [text for _line, text in code_only(CODE / "make_panels.py")]
    assert "rglob" not in tokens
    assert "glob" not in tokens


# ── every module still imports ───────────────────────────────────────────────

@pytest.mark.parametrize("relative", source_files(), ids=lambda p: p.name)
def test_every_module_compiles(relative):
    """A rename that missed an importer shows up here rather than three stages in."""
    result = subprocess.run(
        [sys.executable, "-m", "py_compile", str(REPO / relative)],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-1500:]
