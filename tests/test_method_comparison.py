"""Tests for step 10e: the cross-phase comparison.

The step's arithmetic is thin -- margins and a join -- so what is tested is the
part that could go wrong silently: that the two methods' effects are computed by
the *same* function on the *same* footing, and that a margin means what the
figure's axis label says.
"""

from __future__ import annotations

import numpy as np
import pytest

import method_comparison as methods
import replay_comparison_hmm as comparison


def test_phase1_and_phase4_effects_use_the_same_statistic():
    """Both sides of the comparison come from one `compare`, not two definitions."""
    assert methods.compare is comparison.compare


def test_the_comparison_is_run_on_the_column_cycle_shuffle_for_both_phases():
    assert methods.SHUFFLE == "column_cycle"


def test_a_margin_is_the_effect_minus_its_own_control():
    effect = comparison.compare(np.zeros(40), np.ones(40))["rank_biserial"]
    control = comparison.compare(np.zeros(40), np.zeros(40))["rank_biserial"]
    assert effect - control == pytest.approx(1.0)


def test_a_control_that_runs_the_other_way_inflates_the_margin():
    """Cicero_09012014's identity shuffle is negative, so its margin exceeds 1.

    Worth pinning because a margin above the statistic's own range is a signal
    that the control, not the effect, is doing the work -- which is exactly what
    D8.8 found on that session and what the 10e write-up has to say about it.
    """
    effect = comparison.compare(np.zeros(40), np.ones(40))["rank_biserial"]
    control = comparison.compare(np.ones(40), np.zeros(40))["rank_biserial"]
    assert effect - control == pytest.approx(2.0)
    assert control < 0
