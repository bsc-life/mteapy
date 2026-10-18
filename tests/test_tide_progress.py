import pytest

from mteapy.context_scoring import RunCancelled, build_complex_cache
from mteapy.tide import calculate_random_TIDE_scores_context_aware


@pytest.fixture
def inputs(toy_model):
    routes = {"1": {0: frozenset({"R1", "R2"}), 1: frozenset({"R3"})}}
    return routes, build_complex_cache(toy_model, {"R1", "R2", "R3"})


GENES = {"g1": 2.0, "g2": -1.5, "g3": 3.0, "g4": -0.5}


def test_progress_is_reported_after_every_permutation(inputs):
    routes, cache = inputs
    seen = []
    with_cb = calculate_random_TIDE_scores_context_aware(
        GENES, routes, cache, "absmax", n_permutations=7, random_seed=1, progress=lambda d, t: seen.append((d, t)))
    assert seen == [(i, 7) for i in range(1, 8)]
    without = calculate_random_TIDE_scores_context_aware(GENES, routes, cache, "absmax", n_permutations=7, random_seed=1)
    assert with_cb.equals(without)


def test_progress_can_cancel_a_serial_run(inputs):
    routes, cache = inputs

    def stop(done, total):
        if done == 3:
            raise RunCancelled()
    with pytest.raises(RunCancelled):
        calculate_random_TIDE_scores_context_aware(GENES, routes, cache, "absmax", n_permutations=50,
                                                   random_seed=1, progress=stop)


def test_parallel_run_gives_the_same_null_and_reports_progress(inputs):
    routes, cache = inputs
    seen = []
    parallel = calculate_random_TIDE_scores_context_aware(
        GENES, routes, cache, "absmax", n_permutations=8, n_cpus=2, random_seed=5,
        progress=lambda d, t: seen.append(d))
    serial = calculate_random_TIDE_scores_context_aware(GENES, routes, cache, "absmax", n_permutations=8,
                                                        n_cpus=2, random_seed=5)
    assert seen == list(range(1, 9)) and parallel.equals(serial)
