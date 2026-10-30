"""Tests for the step 10c state-selection rules.

The selection logic is where this milestone's judgement lives -- which m, on
what statistic, and when to refuse to name one -- so it is tested on hand-built
tables whose right answer is obvious by inspection.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import bout_phase_baseline as phase
import fit_hmm
import select_hmm_states as select


def sweep(claim_a_by_fold, likelihood, grid=(6, 8, 10), shuffled=0.0, n_restarts=4):
    """A fits table: ``claim_a_by_fold[m][fold]`` is that cell's constant value."""
    rows = []
    for m in grid:
        for fold in range(len(claim_a_by_fold[m])):
            for seed in range(n_restarts):
                rows.append({"session": "S", "kind": "real", "n_states": m,
                             "fold": fold, "seed": seed,
                             "claim_a_r": claim_a_by_fold[m][fold],
                             "held_out_log_likelihood": likelihood[m],
                             "median_error_cm": 10.0, "chance_error_cm": 50.0,
                             "median_dwell_s": 0.5, "claim_b_position": 0.3,
                             "claim_b_spectral": 0.2, "converged": True,
                             "crossing_time_s": 0.3})
                rows.append({**rows[-1], "kind": "shuffled",
                             "claim_a_r": shuffled, "median_error_cm": 40.0})
    return pd.DataFrame(rows)


def test_claim_a_picks_its_peak_and_the_likelihood_pick_is_only_recorded():
    frame = sweep({6: [0.3, 0.3, 0.3], 8: [0.8, 0.8, 0.8], 10: [0.5, 0.5, 0.5]},
                  {6: -300.0, 8: -200.0, 10: -100.0})
    table = select.summarise(frame, [6, 8, 10])
    decision = select.choose(table, frame)
    assert decision["m_by_claim_a"] == 8
    assert decision["m_by_likelihood"] == 10
    assert decision["criteria_agree"] is False


def test_a_monotone_likelihood_is_reported_as_still_rising_at_the_cap():
    frame = sweep({6: [0.8, 0.8, 0.8], 8: [0.5, 0.5, 0.5], 10: [0.3, 0.3, 0.3]},
                  {6: -300.0, 8: -200.0, 10: -100.0})
    decision = select.choose(select.summarise(frame, [6, 8, 10]), frame)
    assert decision["likelihood_monotone"] is True
    assert decision["likelihood_rising_at_cap"] is True
    assert decision["likelihood_monotone_run"] == 3


def test_agreement_at_the_top_of_the_grid_is_still_flagged_as_truncation():
    """Both criteria pick the largest m -- which is not independent evidence."""
    frame = sweep({6: [0.3, 0.3, 0.3], 8: [0.5, 0.5, 0.5], 10: [0.8, 0.8, 0.8]},
                  {6: -300.0, 8: -200.0, 10: -100.0})
    decision = select.choose(select.summarise(frame, [6, 8, 10]), frame)
    assert decision["criteria_agree"] is True
    assert decision["likelihood_rising_at_cap"] is True


def test_folds_that_peak_in_different_places_are_called_unstable():
    frame = sweep({6: [0.9, 0.1, 0.1], 8: [0.1, 0.2, 0.2], 10: [0.2, 0.9, 0.9]},
                  {6: -300.0, 8: -200.0, 10: -100.0})
    decision = select.choose(select.summarise(frame, [6, 8, 10]), frame)
    assert decision["folds_unstable"] is True


def test_folds_that_peak_one_step_apart_are_not_called_unstable():
    frame = sweep({6: [0.1, 0.1, 0.1], 8: [0.9, 0.5, 0.9], 10: [0.5, 0.9, 0.5]},
                  {6: -300.0, 8: -200.0, 10: -100.0})
    decision = select.choose(select.summarise(frame, [6, 8, 10]), frame)
    assert decision["folds_unstable"] is False


def test_a_flat_curve_is_marked_flatter_than_the_restart_noise():
    """The argmax of a flat curve moves on nothing, and the table has to say so."""
    rng = np.random.default_rng(0)
    rows = []
    for m in (6, 8, 10):
        for fold in range(3):
            for seed in range(8):
                rows.append({"session": "S", "kind": "real", "n_states": m,
                             "fold": fold, "seed": seed,
                             "claim_a_r": 0.5 + rng.normal(0, 0.2),
                             "held_out_log_likelihood": -100.0 * (4 - m / 4),
                             "median_error_cm": 10.0, "chance_error_cm": 50.0,
                             "median_dwell_s": 0.5, "claim_b_position": 0.3,
                             "claim_b_spectral": 0.2, "converged": True,
                             "crossing_time_s": 0.3})
                rows.append({**rows[-1], "kind": "shuffled", "claim_a_r": 0.0})
    frame = pd.DataFrame(rows)
    decision = select.choose(select.summarise(frame, [6, 8, 10]), frame)
    assert decision["curve_flatter_than_noise"] is True


def test_the_selection_statistic_is_the_median_not_the_best_restart():
    """One lucky restart must not carry an m that is otherwise mediocre."""
    rows = []
    for m in (6, 8):
        for fold in range(3):
            for seed in range(5):
                lucky = (m == 8 and seed == 0)
                rows.append({"session": "S", "kind": "real", "n_states": m,
                             "fold": fold, "seed": seed,
                             "claim_a_r": 0.99 if lucky else (0.2 if m == 8 else 0.5),
                             "held_out_log_likelihood": -100.0,
                             "median_error_cm": 10.0, "chance_error_cm": 50.0,
                             "median_dwell_s": 0.5, "claim_b_position": 0.3,
                             "claim_b_spectral": 0.2, "converged": True,
                             "crossing_time_s": 0.3})
                rows.append({**rows[-1], "kind": "shuffled", "claim_a_r": 0.0})
    frame = pd.DataFrame(rows)
    table = select.summarise(frame, [6, 8])
    assert float(table[table.n_states == 8].claim_a_max.iloc[0]) == pytest.approx(0.99)
    assert select.choose(table, frame)["m_by_claim_a"] == 6


def test_the_shuffled_control_is_summarised_against_the_real_fits():
    frame = sweep({6: [0.4, 0.4, 0.4], 8: [0.6, 0.6, 0.6], 10: [0.5, 0.5, 0.5]},
                  {6: -300.0, 8: -200.0, 10: -100.0}, shuffled=0.45)
    decision = select.choose(select.summarise(frame, [6, 8, 10]), frame)
    assert decision["margin_at_choice"] == pytest.approx(0.15)
    assert decision["shuffled_exceeds_real_frac"] == 0.0


def test_a_control_that_matches_the_data_is_not_hidden_by_the_margin():
    """Every shuffled fit beating every real one has to show up as 100 %."""
    frame = sweep({6: [0.4, 0.4, 0.4], 8: [0.5, 0.5, 0.5], 10: [0.4, 0.4, 0.4]},
                  {6: -300.0, 8: -200.0, 10: -100.0}, shuffled=0.9)
    decision = select.choose(select.summarise(frame, [6, 8, 10]), frame)
    assert decision["shuffled_exceeds_real_frac"] == 1.0
    assert decision["margin_at_choice"] < 0


def test_the_state_grid_is_truncated_at_each_session_s_cap():
    for name, cap in fit_hmm.SESSION_CAPS.items():
        grid = fit_hmm.grid_for(name)
        assert grid, name
        assert max(grid) <= cap
        assert set(grid) <= set(fit_hmm.STATE_GRID)


def test_a_session_with_no_recorded_cap_is_refused_rather_than_guessed():
    with pytest.raises(SystemExit):
        fit_hmm.grid_for("Not_a_session")


def test_the_bout_phase_map_is_the_mean_training_position_at_each_index():
    positions = [np.array([0.0, 1.0, 2.0]), np.array([2.0, 3.0, 4.0])]
    mapped = phase.phase_map(positions, 3, circular=False, extent=10.0)
    assert list(mapped) == [1.0, 2.0, 3.0]


def test_the_bout_phase_map_pools_indices_past_the_cap():
    positions = [np.array([0.0, 10.0, 20.0, 30.0])]
    mapped = phase.phase_map(positions, 2, circular=False, extent=100.0)
    assert mapped[0] == pytest.approx(0.0)
    assert mapped[1] == pytest.approx(20.0)  # bins 1, 2 and 3 pooled
