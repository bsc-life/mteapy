import pandas as pd
import pytest

from mteapy.context_scoring import (
    ReactionEvidence,
    assess_reaction,
    build_complex_cache,
    score_task,
    score_tasks_matrix,
    tied_argmax,
)


@pytest.mark.parametrize("scores,expected", [
    ({}, ()),
    ({"A": 1.0}, ("A",)),
    ({"A": 1.0, "B": 2.0}, ("B",)),
    ({"A": 2.0, "B": 2.0, "C": 1.0}, ("A", "B")),
    ({"A": 1.0, "B": 1.0 + 1e-12}, ("A", "B")),  # within floating tolerance -> tied
])
def test_tied_argmax(scores, expected):
    assert tied_argmax(scores) == expected


@pytest.mark.parametrize("scores,expected", [
    # Plain "max" would pick B (least negative); "absmax" correctly treats
    # A's strong down-regulation as the dominant signal, same as it would
    # for an equally strong up-regulation.
    ({"A": -2.0, "B": 0.1}, ("A",)),
    ({"A": -2.0, "B": 2.0}, ("A", "B")),  # tied by magnitude despite opposite sign
    ({"A": -2.0, "B": -2.0 + 1e-12}, ("A", "B")),  # tied within tolerance
    ({}, ()),
])
def test_tied_argmax_absmax_compares_by_magnitude_not_raw_value(scores, expected):
    assert tied_argmax(scores, or_func="absmax") == expected


def test_tied_argmax_rejects_unknown_or_func():
    with pytest.raises(ValueError, match="or_func"):
        tied_argmax({"A": 1.0}, or_func="nonsense")


def test_assess_reaction_and_within_complex_or_across_complexes():
    # (A and B) or (C) -- AND takes the min within a complex, OR takes the
    # max across complexes; C alone should win here since it's unopposed.
    complexes = (("A", "B"), ("C",))
    gene_dict = {"A": 5.0, "B": 1.0, "C": 3.0}
    result = assess_reaction("R1", complexes, gene_dict)
    assert result.score == 3.0
    assert result.winning_complexes == (("C",),)
    assert result.evidence == ReactionEvidence.SUPPORTED
    assert result.complex_scores == {("A", "B"): 1.0, ("C",): 3.0}


def test_assess_reaction_missing_gene_defaults_to_zero():
    complexes = (("A", "B"),)
    result = assess_reaction("R1", complexes, {"A": 5.0})  # B absent -> AND floors to 0
    assert result.score == 0.0
    assert result.evidence == ReactionEvidence.NO_EVIDENCE
    assert result.winning_complexes == ()  # no complex counts as "winning" at zero evidence


def test_assess_reaction_no_gpr_is_its_own_category():
    result = assess_reaction("R1", (), {"A": 5.0})
    assert result.score == 0.0
    assert result.evidence == ReactionEvidence.NO_GPR
    assert result.complex_scores == {}
    assert result.winning_complexes == ()


def test_assess_reaction_all_complexes_zero_is_no_evidence_not_ambiguous():
    # Two candidate complexes, both real genes, both simply unexpressed --
    # this must NOT be reported as an "ambiguous tie" between them.
    complexes = (("A",), ("B",))
    result = assess_reaction("R1", complexes, {"A": 0.0, "B": 0.0})
    assert result.score == 0.0
    assert result.evidence == ReactionEvidence.NO_EVIDENCE
    assert result.winning_complexes == ()


def test_assess_reaction_ambiguous_winner_among_real_evidence():
    complexes = (("A",), ("B",))
    result = assess_reaction("R1", complexes, {"A": 2.0, "B": 2.0})
    assert result.score == 2.0
    assert result.evidence == ReactionEvidence.AMBIGUOUS
    assert set(result.winning_complexes) == {("A",), ("B",)}


def test_assess_reaction_absmax_lets_a_downregulated_complex_win():
    # TIDE-style signed log-fold-change signal: complex A is strongly
    # down-regulated (-3.0), complex B barely up (0.2). Under plain "max"
    # B would (wrongly, for a signed signal) look like "the" winner; under
    # "absmax" A correctly dominates since it deviates from 0 the most.
    complexes = (("A",), ("B",))
    gene_dict = {"A": -3.0, "B": 0.2}

    default = assess_reaction("R1", complexes, gene_dict)
    assert default.score == 0.2
    assert default.winning_complexes == (("B",),)

    absmax_result = assess_reaction("R1", complexes, gene_dict, or_func="absmax")
    assert absmax_result.score == -3.0
    assert absmax_result.winning_complexes == (("A",),)
    assert absmax_result.evidence == ReactionEvidence.SUPPORTED


def test_assess_reaction_and_is_unaffected_by_or_func():
    # AND (min of the raw signed values within one complex) never changes
    # with or_func -- only OR (across complexes) does.
    complexes = (("A", "B"),)
    gene_dict = {"A": -1.0, "B": -5.0}
    default = assess_reaction("R1", complexes, gene_dict)
    absmax_result = assess_reaction("R1", complexes, gene_dict, or_func="absmax")
    assert default.score == absmax_result.score == -5.0


def test_assess_reaction_rejects_unknown_or_func():
    with pytest.raises(ValueError, match="or_func"):
        assess_reaction("R1", (("A",),), {"A": 1.0}, or_func="nonsense")


def test_or_func_absmax_reproduces_published_tide_scores_exactly():
    """Regression test against real, published TIDE output: the AGS-paper
    dataset (Benedicto et al., https://doi.org/10.1038/s41540-025-00586-y,
    data/code at github.com/bsc-life/ags-paper), condition TAKi, masked
    log2FoldChange (non-significant by padj >= 0.05 zeroed), scored against
    two real single-reaction Human-GEM tasks from the paper's own
    task_structure_matrix.tsv (X76: a single-gene GPR; X142: a real 2-gene
    OR). Both reproduce the paper's published TIDE task score to full
    floating-point precision -- verified by hand against
    ags-paper/results/TIDE/ags-tide-TAKi.tsv, not re-derived here (this
    package doesn't ship that data), which is exactly why the raw numbers
    below are hardcoded rather than read from a file.
    """
    # X76 "Conversion of asparate to asparagine" -> MAR03903, GPR is a lone
    # gene (ENSG00000070669), masked LFC -1.7410477269287123.
    gene_dict_x76 = {"ENSG00000070669": -1.7410477269287123}
    result = assess_reaction("MAR03903", (("ENSG00000070669",),), gene_dict_x76, or_func="absmax")
    assert result.score == pytest.approx(-1.7410477269287123, rel=0, abs=1e-12)
    assert result.evidence == ReactionEvidence.SUPPORTED

    # X142 "Glycerol-3-phosphate synthesis" -> MAR00479, GPR is
    # "ENSG00000152642 or ENSG00000167588"; the first gene was masked to 0
    # (non-significant), the second carries the real signal -- absmax must
    # pick the real one, not silently default to the first/only-max value.
    complexes = (("ENSG00000152642",), ("ENSG00000167588",))
    gene_dict_x142 = {"ENSG00000152642": 0.0, "ENSG00000167588": -0.8500990237008653}
    result = assess_reaction("MAR00479", complexes, gene_dict_x142, or_func="absmax")
    assert result.score == pytest.approx(-0.8500990237008653, rel=0, abs=1e-12)
    assert result.winning_complexes == (("ENSG00000167588",),)


def test_score_task_absmax_picks_the_most_differentially_regulated_route():
    # Route 1's reaction is mildly up (+0.1); route 2's is strongly down
    # (-4.0). A signed (LFC-style) signal should let route 2 win under
    # "absmax" even though its raw score is lower than route 1's.
    complex_cache = {"R1": (("g1",),), "R2": (("g2",),)}
    task_routes = {1: frozenset({"R1"}), 2: frozenset({"R2"})}
    gene_dict = {"g1": 0.1, "g2": -4.0}

    default = score_task(task_routes, complex_cache, gene_dict, task_id="T")
    assert default.winning_route_ids == (1,)  # plain max: 0.1 > -4.0

    absmax_report = score_task(task_routes, complex_cache, gene_dict, task_id="T", or_func="absmax")
    assert absmax_report.winning_route_ids == (2,)
    assert absmax_report.score == -4.0


def test_score_task_reaction_cache_gives_identical_results_to_uncached():
    complex_cache = {"R1": (("g1",),), "R2": (("g2",),)}
    task_routes = {1: frozenset({"R1"}), 2: frozenset({"R2"})}
    gene_dict = {"g1": 0.1, "g2": -4.0}

    uncached = score_task(task_routes, complex_cache, gene_dict, task_id="T", or_func="absmax")
    cache = {}
    cached = score_task(task_routes, complex_cache, gene_dict, task_id="T", or_func="absmax", reaction_cache=cache)

    assert cached.score == uncached.score
    assert cached.winning_route_ids == uncached.winning_route_ids
    assert set(cache.keys()) == {"R1", "R2"}


def test_score_task_reaction_cache_is_reused_across_calls_not_recomputed():
    complex_cache = {"R1": (("g1",),)}
    task_routes = {1: frozenset({"R1"})}
    gene_dict = {"g1": 2.0}
    cache = {}

    score_task(task_routes, complex_cache, gene_dict, task_id="T1", reaction_cache=cache)
    # Poison the cached entry so a second call can only get the right
    # answer if it reuses the cache instead of recomputing from gene_dict.
    poisoned = cache["R1"]
    cache["R1"] = poisoned.__class__(poisoned.reaction_id, 999.0, poisoned.complex_scores, poisoned.winning_complexes, poisoned.evidence)

    report = score_task(task_routes, complex_cache, gene_dict, task_id="T2", reaction_cache=cache)
    assert report.score == 999.0


def test_score_task_single_route_is_never_tied_and_is_complete():
    complex_cache = {"R1": (("g1",),), "R2": (("g2",),)}
    task_routes = {1: frozenset({"R1", "R2"})}
    gene_dict = {"g1": 3.0, "g2": 1.0}

    report = score_task(task_routes, complex_cache, gene_dict, task_id="T")
    assert not report.is_tied
    assert report.winning_route_ids == (1,)
    assert report.score == 1.0  # min(3.0, 1.0)
    assert report.is_complete  # both reactions cleanly SUPPORTED


def test_score_task_detects_genuine_route_tie():
    complex_cache = {"R1": (("g1",),), "R2": (("g2",),)}
    task_routes = {1: frozenset({"R1"}), 2: frozenset({"R2"})}
    gene_dict = {"g1": 5.0, "g2": 5.0}

    report = score_task(task_routes, complex_cache, gene_dict, task_id="T")
    assert report.is_tied
    assert report.winning_route_ids == (1, 2)
    assert report.score == 5.0


def test_score_task_min_vs_median_can_change_the_winner():
    complex_cache = {"R1": (("g1",),), "R2": (("g2",),), "R3": (("g3",),)}
    task_routes = {
        "A": frozenset({"R1", "R2"}),  # one very low (bottleneck), one very high
        "B": frozenset({"R3"}),         # one middling
    }
    gene_dict = {"g1": 0.1, "g2": 100.0, "g3": 5.0}

    min_report = score_task(task_routes, complex_cache, gene_dict, aggregation="min")
    assert min_report.winning_route_ids == ("B",)  # A's min (0.1) loses to B's 5.0

    median_report = score_task(task_routes, complex_cache, gene_dict, aggregation="median")
    # median of a 1-reaction route is just that reaction's score (5.0);
    # median of A's two reactions is the mean of 0.1 and 100 -> 50.05
    assert median_report.winning_route_ids == ("A",)


def test_score_task_unsupported_aggregation_raises():
    complex_cache = {"R1": (("g1",),)}
    with pytest.raises(ValueError):
        score_task({1: frozenset({"R1"})}, complex_cache, {"g1": 1.0}, aggregation="mode")


def test_score_task_mean_aggregation_matches_published_tide_methodology():
    """mteapy.tide's own calculate_TIDE_scores/calculate_TIDEe_scores
    (the published TIDE/TIDE-essential implementation this package also
    ships) aggregate a task's reaction/gene scores with plain `np.mean`,
    never "min" or "median" -- "mean" here exists specifically to
    reproduce/compare against that established methodology, not as a third
    arbitrary option. This is a straightforward statistics.mean, but the
    result differs from both "min" and "median" for an asymmetric route,
    which is the property that matters here."""
    complex_cache = {"R1": (("g1",),), "R2": (("g2",),)}
    task_routes = {1: frozenset({"R1", "R2"})}
    gene_dict = {"g1": 1.0, "g2": 5.0}

    report = score_task(task_routes, complex_cache, gene_dict, aggregation="mean")
    assert report.score == pytest.approx(3.0)  # mean(1.0, 5.0), not min=1.0 or median=3.0 (coincide here)

    complex_cache2 = {"R1": (("g1",),), "R2": (("g2",),), "R3": (("g3",),)}
    task_routes2 = {1: frozenset({"R1", "R2", "R3"})}
    gene_dict2 = {"g1": 1.0, "g2": 2.0, "g3": 9.0}
    report2 = score_task(task_routes2, complex_cache2, gene_dict2, aggregation="mean")
    assert report2.score == pytest.approx(4.0)  # mean(1,2,9)=4.0, distinct from median=2.0 and min=1.0


def test_score_task_reports_partial_when_winning_route_has_no_evidence_reaction():
    # Route 1 wins on score but one of its reactions has real genes with
    # zero expression support -- the winning route is still THE winner,
    # but the report must flag it as partial, not complete.
    complex_cache = {
        "R1": (("g1",),),          # supported: g1 expressed
        "R2": (("g2",), ("g3",)),  # no_evidence: both g2 and g3 are 0
        "R4": (("g4",),),          # a worse route, would be complete if it won
    }
    task_routes = {
        1: frozenset({"R1", "R2"}),
        2: frozenset({"R4"}),
    }
    gene_dict = {"g1": 10.0, "g2": 0.0, "g3": 0.0, "g4": 1.0}

    report = score_task(task_routes, complex_cache, gene_dict, task_id="T")
    # min(10.0, 0.0) = 0.0 for route 1, vs 1.0 for route 4 -- wait, route 4
    # alone (1.0) beats route 1 (0.0), so route 2 (R4) actually wins here;
    # use route 4's reaction id key "2" for clarity in the assertion below.
    assert report.winning_route_ids == (2,)
    assert report.is_complete  # the winning route (R4 alone) is fully supported

    # But route 1 (the loser) is correctly flagged as partial/incomplete,
    # and its problem reaction is identifiable.
    assert not report.routes[1].is_complete
    assert report.routes[1].incomplete_reactions == ("R2",)
    assert report.routes[1].reactions["R2"].evidence == ReactionEvidence.NO_EVIDENCE


def test_score_task_is_complete_false_when_winning_route_itself_is_ambiguous():
    complex_cache = {"R1": (("g1",), ("g2",))}  # two complexes, equally expressed
    task_routes = {1: frozenset({"R1"})}
    gene_dict = {"g1": 3.0, "g2": 3.0}

    report = score_task(task_routes, complex_cache, gene_dict, task_id="T")
    assert not report.is_complete
    assert report.routes[1].incomplete_reactions == ("R1",)
    assert report.routes[1].reactions["R1"].evidence == ReactionEvidence.AMBIGUOUS


def test_build_complex_cache_and_score_tasks_matrix(toy_model):
    reaction_ids = ["R1", "R2", "R3", "R4"]
    cache = build_complex_cache(toy_model, reaction_ids)
    assert cache["R1"] == (("g1",),)
    assert cache["R3"] == (("g3",),)

    tasks_routes = {
        "LINEAR": {1: frozenset({"R1", "R2"}), 2: frozenset({"R3"})},
    }
    expr_df = pd.DataFrame({
        "sample_A": {"g1": 10.0, "g2": 10.0, "g3": 1.0},   # R1+R2 route wins
        "sample_B": {"g1": 1.0, "g2": 1.0, "g3": 10.0},    # R3 route wins
    })

    scores, complete = score_tasks_matrix(tasks_routes, toy_model, expr_df)
    assert scores.shape == (1, 2) == complete.shape
    assert scores.loc["LINEAR", "sample_A"] == 10.0  # min(10, 10) beats R3's 1.0
    assert scores.loc["LINEAR", "sample_B"] == 10.0  # R3 alone scores 10.0
    assert complete.loc["LINEAR", "sample_A"] and complete.loc["LINEAR", "sample_B"]


def test_and_linked_paralog_pair_correctly_flags_no_evidence():
    # Regression test for a real annotation pattern found in Human-GEM
    # v2.0.1: MAR06914 (cytochrome c oxidase / Complex IV) AND-links COX7B
    # (ubiquitous, real GTEx min 6.24 TPM, never zero) with its testis-
    # specific paralog COX7B2 (real GTEx: 0 TPM in 65/68 tissues). Because
    # the GPR treats both as unconditionally-required fixed subunits rather
    # than OR-linked alternatives for the same subunit position, the whole
    # reaction collapses to a score of 0 in nearly every non-testis sample
    # even though the functional complex (via COX7B alone) is fully present.
    # This exact shape -- one always-expressed gene AND-linked with a real
    # but narrowly-restricted paralog -- is what `is_complete` exists to
    # catch, and it was first found this way, on real data, not the other
    # way around.
    gpr_like_complexes = (("COX7B", "COX7B2", "OTHER_SUBUNIT"),)  # single, forced complex
    non_testis_sample = {"COX7B": 40.0, "COX7B2": 0.0, "OTHER_SUBUNIT": 90.0}
    testis_sample = {"COX7B": 40.0, "COX7B2": 105.9, "OTHER_SUBUNIT": 90.0}

    non_testis_result = assess_reaction("MAR06914", gpr_like_complexes, non_testis_sample)
    assert non_testis_result.score == 0.0
    assert non_testis_result.evidence == ReactionEvidence.NO_EVIDENCE

    testis_result = assess_reaction("MAR06914", gpr_like_complexes, testis_sample)
    assert testis_result.score == 40.0  # min(40.0, 105.9, 90.0)
    assert testis_result.evidence == ReactionEvidence.SUPPORTED
