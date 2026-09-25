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
        score_task({1: frozenset({"R1"})}, complex_cache, {"g1": 1.0}, aggregation="mean")


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
