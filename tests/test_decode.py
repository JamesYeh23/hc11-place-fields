"""Bayesian decoder: likelihood, cross-validation integrity, and recovery."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import requires_data
from hc11 import config
from hc11.decode import (
    chance_baseline,
    decode,
    decode_log_posterior,
    error_vs_ensemble_size,
    posterior_mode,
    prepare_decoder,
    tuning_from_laps,
)
from hc11.fields import select_running_samples, spike_slots
from hc11.laps import segment_laps, segment_laps_from_config
from hc11.track import Track
from test_fields import LINEAR, _place_cell_session

LAP_KW = dict(
    position_smoothing_sigma_s=0.25, direction_threshold_cm_s=5.0,
    min_coverage_frac=0.5, min_lap_duration_s=0.5, max_lap_duration_s=60.0,
)


# ---------------------------------------------------------------------------
# The likelihood itself
# ---------------------------------------------------------------------------


class TestLogPosterior:
    def test_columns_are_normalised(self):
        rng = np.random.default_rng(0)
        tuning = rng.uniform(0.1, 20, size=(6, 12))
        counts = rng.poisson(1.0, size=(6, 5)).astype(float)
        log_post = decode_log_posterior(counts, tuning, 0.25)
        np.testing.assert_allclose(np.exp(log_post).sum(axis=0), 1.0, rtol=1e-12)

    def test_spikes_matching_a_cells_field_point_at_it(self):
        tuning = np.full((3, 10), 0.1)
        tuning[1, 7] = 12.0  # cell 1 fires only at bin 7
        counts = np.zeros((3, 1))
        counts[1, 0] = 3.0  # 12 Hz * 0.25 s = 3 expected spikes there
        assert int(np.argmax(decode_log_posterior(counts, tuning, 0.25)[:, 0])) == 7

    def test_too_few_spikes_is_evidence_against_a_high_rate_position(self):
        """The -dt * sum_i f_i(x) term, which a log-likelihood-only decoder loses.

        One spike where the field predicts 12.5 argues against that position, not
        for it: the Poisson likelihood of 1 given 12.5 is tiny.
        """
        tuning = np.full((1, 10), 0.1)
        tuning[0, 7] = 50.0
        counts = np.zeros((1, 1))
        counts[0, 0] = 1.0
        post = np.exp(decode_log_posterior(counts, tuning, 0.25)[:, 0])
        assert int(np.argmax(post)) != 7
        assert post[7] < post[0]

    def test_two_cells_agree_on_the_overlap(self):
        tuning = np.full((2, 10), 0.1)
        tuning[0, 3:7] = 20.0
        tuning[1, 5:9] = 20.0
        counts = np.ones((2, 1))
        post = np.exp(decode_log_posterior(counts, tuning, 0.25)[:, 0])
        assert 5 <= int(np.argmax(post)) <= 6, "the peak belongs in the overlap"

    def test_silence_favours_low_rate_positions(self):
        """With no spikes the only evidence is the exp(-dt*sum f) term."""
        tuning = np.array([[1.0, 10.0, 1.0]])
        post = np.exp(decode_log_posterior(np.zeros((1, 1)), tuning, 1.0)[:, 0])
        assert post[1] < post[0] and post[1] < post[2]
        assert post[0] == pytest.approx(post[2])

    def test_more_spikes_sharpen_the_posterior(self):
        tuning = np.full((4, 12), 0.5)
        tuning[:, 6] = 20.0
        one = np.exp(decode_log_posterior(np.ones((4, 1)), tuning, 0.25)[:, 0]).max()
        many = np.exp(decode_log_posterior(np.full((4, 1), 5.0), tuning, 0.25)[:, 0]).max()
        assert many > one

    def test_zero_rate_bins_do_not_produce_minus_inf(self):
        tuning = np.zeros((2, 5))
        tuning[0, 2] = 10.0
        log_post = decode_log_posterior(np.ones((2, 1)), tuning, 0.25, rate_floor_hz=0.01)
        assert np.all(np.isfinite(log_post))

    def test_nan_tuning_is_floored_not_propagated(self):
        tuning = np.full((2, 5), 1.0)
        tuning[0, 3] = np.nan
        assert np.all(np.isfinite(decode_log_posterior(np.ones((2, 1)), tuning, 0.25)))

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="cells"):
            decode_log_posterior(np.ones((3, 2)), np.ones((4, 5)), 0.25)

    def test_posterior_mode_returns_bin_centres(self):
        log_post = np.log(np.array([[0.1, 0.8], [0.9, 0.2]]))
        centres = np.array([0.5, 1.5])
        np.testing.assert_array_equal(posterior_mode(log_post, centres), [1.5, 0.5])


# ---------------------------------------------------------------------------
# Cross-validation integrity -- the decoder must never see the lap it decodes
# ---------------------------------------------------------------------------


class TestNoLeakage:
    def _inputs(self, tmp_path, arrays, **kw):
        session = _place_cell_session(tmp_path, arrays, 0.8, **kw)
        laps = segment_laps(session, LINEAR, **LAP_KW)
        mask = np.ones(len(session.position_t), dtype=bool)
        samples = select_running_samples(session, LINEAR, laps, mask, "pos")
        params = config.load_params()
        inputs = prepare_decoder(
            session, LINEAR, samples, laps, np.array([102]), params, time_bin_s=0.25
        )
        return session, laps, samples, params, inputs

    def test_template_for_a_lap_excludes_that_lap(self, tmp_path, synthetic_arrays):
        session, laps, samples, params, inputs = self._inputs(tmp_path, synthetic_arrays)
        slots = {102: spike_slots(session, 102, samples)}
        for lap in inputs.laps:
            held_out = samples.lap_of_sample == lap.index
            expected = tuning_from_laps(
                session, LINEAR, samples, np.array([102]), ~held_out,
                sigma_bins=params["place_fields"]["smoothing_sigma_bins"],
                min_occupancy_s=params["place_fields"]["min_occupancy_s"],
                slots_by_cell=slots,
            )
            np.testing.assert_array_equal(inputs.tuning_by_lap[lap.index], expected)

    def test_a_laps_own_spikes_cannot_reach_its_template(self, tmp_path, synthetic_arrays):
        """Inject a burst into one lap; its own template must not move."""
        session, laps, samples, params, inputs = self._inputs(tmp_path, synthetic_arrays)
        target = inputs.laps[0]
        held_out = samples.lap_of_sample == target.index

        # Fabricate spikes at a position this cell is otherwise silent at, but
        # only during the held-out lap.
        in_lap = samples.indices[held_out]
        fake_times = session.position_t[in_lap]
        merged = np.sort(np.concatenate([session.spike_times, fake_times]))
        order = np.argsort(np.concatenate([session.spike_times, fake_times]), kind="stable")
        merged_ids = np.concatenate(
            [session.spike_ids, np.full(len(fake_times), 102)]
        )[order]
        import dataclasses

        spiked = dataclasses.replace(session, spike_times=merged, spike_ids=merged_ids)

        spiked_inputs = prepare_decoder(
            spiked, LINEAR, samples, laps, np.array([102]), params, time_bin_s=0.25
        )
        # The held-out lap's own template is unchanged by its own extra spikes...
        np.testing.assert_allclose(
            spiked_inputs.tuning_by_lap[target.index],
            inputs.tuning_by_lap[target.index],
        )
        # ...while every other lap's template, and the all-laps template, do change.
        assert not np.allclose(spiked_inputs.full_tuning, inputs.full_tuning)
        other = [lap for lap in inputs.laps if lap.index != target.index][0]
        assert not np.allclose(
            spiked_inputs.tuning_by_lap[other.index], inputs.tuning_by_lap[other.index]
        )

    def test_every_lap_is_decoded_exactly_once(self, tmp_path, synthetic_arrays):
        _, _, _, params, inputs = self._inputs(tmp_path, synthetic_arrays)
        result = decode(inputs, LINEAR, params)
        decoded_laps = set(np.unique(result.lap_of_bin))
        assert decoded_laps == {lap.index for lap in inputs.laps}


# ---------------------------------------------------------------------------
# Recovery on synthetic ensembles
# ---------------------------------------------------------------------------


def _ensemble_session(tmp_path, arrays, n_cells=12, seed=0, peak_hz=15.0):
    """Several place cells tiling the track, as one synthetic session."""
    rng = np.random.default_rng(seed)
    dt, n_laps = 0.0256, 24
    lengths = rng.integers(90, 150, n_laps)
    x = np.concatenate([
        np.linspace(0, 1.6, n) if i % 2 == 0 else np.linspace(1.6, 0, n)
        for i, n in enumerate(lengths)
    ])
    t = 100.0 + np.arange(len(x)) * dt
    centres = np.linspace(0.05, 1.55, n_cells)
    ids = 102 + np.arange(n_cells)

    times, cells = [], []
    for cid, centre in zip(ids, centres, strict=True):
        rate = 0.1 + peak_hz * np.exp(-0.5 * ((x - centre) / 0.15) ** 2)
        counts = rng.poisson(rate * dt)
        for i, c in enumerate(counts):
            if c:
                times.append(t[i] + rng.uniform(0, dt, c))
                cells.append(np.full(c, cid))
    times = np.concatenate(times)
    cells = np.concatenate(cells)
    order = np.argsort(times)

    import dataclasses

    from conftest import SYNTH_NAME, write_sessinfo
    from hc11.io import load_session

    arrays = dict(arrays)
    arrays.update(
        SpikeTimes=times[order][:, None], SpikeIDs=cells[order].astype(float)[:, None],
        PyrIDs=ids.astype(float)[:, None], IntIDs=np.array([[1.0]]),
        TimeStamps=t[None, :], TwoDLocation=np.column_stack([x, np.zeros_like(x)]),
        OneDLocation=x[:, None],
        MazeEpoch=np.array([[t[0], t[-1] + dt]]), PREEpoch=np.array([[0.0, t[0]]]),
        POSTEpoch=np.array([[t[-1] + dt, t[-1] + 100.0]]),
        sessDuration=np.array([[t[-1] + 100.0]]),
        Wake=np.array([[0.0, t[-1] + 100.0]]), Drowsy=np.array([[0.0, 1.0]]),
        NREM=np.array([[1.0, 2.0]]), Intermediate=np.array([[2.0, 3.0]]),
        REM=np.array([[3.0, 4.0]]),
    )
    session = load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))
    del dataclasses
    return session, ids


class TestRecovery:
    def test_decodes_synthetic_position_well(self, tmp_path, synthetic_arrays):
        session, ids = _ensemble_session(tmp_path, synthetic_arrays)
        laps = segment_laps(session, LINEAR, **LAP_KW)
        mask = np.ones(len(session.position_t), dtype=bool)
        samples = select_running_samples(session, LINEAR, laps, mask, "pos")
        params = config.load_params()
        inputs = prepare_decoder(session, LINEAR, samples, laps, ids, params, time_bin_s=0.25)
        result = decode(inputs, LINEAR, params)
        assert result.median_error_cm < 20.0, "12 tiling place cells should decode well"
        assert result.fraction_undecodable < 0.1

    def test_cross_validation_costs_accuracy(self, tmp_path, synthetic_arrays):
        """A template that has seen the decoded lap must do at least as well.

        The gap is small -- one lap of 24 barely moves a template -- so this is a
        directional check, not a magnitude one. Its value is as a tripwire: if
        cross-validation were silently bypassed the two would be identical.
        """
        session, ids = _ensemble_session(tmp_path, synthetic_arrays)
        laps = segment_laps(session, LINEAR, **LAP_KW)
        mask = np.ones(len(session.position_t), dtype=bool)
        samples = select_running_samples(session, LINEAR, laps, mask, "pos")
        params = config.load_params()
        inputs = prepare_decoder(session, LINEAR, samples, laps, ids, params, time_bin_s=0.25)
        honest = decode(inputs, LINEAR, params).median_error_cm
        leaky = decode(inputs, LINEAR, params, tuning_override=inputs.full_tuning)
        assert leaky.median_error_cm <= honest
        assert not np.allclose(inputs.full_tuning, inputs.tuning_by_lap[inputs.laps[0].index])

    def test_chance_baseline_is_far_worse(self, tmp_path, synthetic_arrays):
        session, ids = _ensemble_session(tmp_path, synthetic_arrays)
        laps = segment_laps(session, LINEAR, **LAP_KW)
        mask = np.ones(len(session.position_t), dtype=bool)
        samples = select_running_samples(session, LINEAR, laps, mask, "pos")
        params = config.load_params()
        inputs = prepare_decoder(session, LINEAR, samples, laps, ids, params, time_bin_s=0.25)
        real = decode(inputs, LINEAR, params).median_error_cm
        chance = chance_baseline(inputs, LINEAR, params, np.random.default_rng(1))
        assert np.median(chance) > 3 * real

    def test_error_falls_with_ensemble_size(self, tmp_path, synthetic_arrays):
        session, ids = _ensemble_session(tmp_path, synthetic_arrays)
        laps = segment_laps(session, LINEAR, **LAP_KW)
        mask = np.ones(len(session.position_t), dtype=bool)
        samples = select_running_samples(session, LINEAR, laps, mask, "pos")
        params = config.load_params()
        inputs = prepare_decoder(session, LINEAR, samples, laps, ids, params, time_bin_s=0.25)
        curve = error_vs_ensemble_size(inputs, LINEAR, params, np.random.default_rng(2))
        assert len(curve) >= 2
        assert curve[0]["median_error_cm"] > curve[-1]["median_error_cm"]

    def test_confusion_matrix_rows_sum_to_one(self, tmp_path, synthetic_arrays):
        session, ids = _ensemble_session(tmp_path, synthetic_arrays)
        laps = segment_laps(session, LINEAR, **LAP_KW)
        mask = np.ones(len(session.position_t), dtype=bool)
        samples = select_running_samples(session, LINEAR, laps, mask, "pos")
        params = config.load_params()
        inputs = prepare_decoder(session, LINEAR, samples, laps, ids, params, time_bin_s=0.25)
        matrix = decode(inputs, LINEAR, params).confusion_matrix(LINEAR)
        occupied = matrix.sum(axis=1) > 0
        np.testing.assert_allclose(matrix[occupied].sum(axis=1), 1.0)
        # Most mass on or near the diagonal.
        assert np.trace(matrix) / occupied.sum() > 0.4


def test_circular_error_wraps():
    """An error across the reward site must be short, not a full lap."""
    circle = Track(kind="circular", lo=0.0, hi=2.9, bin_size=0.1)
    assert circle.distance(0.05, 2.85) * 100 == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# Real data
# ---------------------------------------------------------------------------


@requires_data
def test_real_session_decodes_well_above_chance(real_sessions):
    from hc11.behavior import running_mask_from_config
    from hc11.fields import analyse_cell

    params = config.load_params()
    session = real_sessions("Achilles_10252013")
    track = Track.from_session(session, params["place_fields"]["bin_size_cm"])
    running = running_mask_from_config(session, params)
    laps = segment_laps_from_config(session, track, params)
    samples = select_running_samples(session, track, laps, running.mask, "pos")
    rng = np.random.default_rng(0)
    cells = [analyse_cell(session, track, samples, c, params, rng) for c in session.pyr_ids]
    ids = np.array(sorted({c.cluster_id for c in cells if c.is_place_cell}))

    inputs = prepare_decoder(session, track, samples, laps, ids, params, time_bin_s=0.25)
    result = decode(inputs, track, params)
    chance = chance_baseline(inputs, track, params, np.random.default_rng(3))
    assert result.median_error_cm < 20.0
    assert np.median(chance) > 5 * result.median_error_cm
    # Cross-validation really is costing something.
    leaky = decode(inputs, track, params, tuning_override=inputs.full_tuning)
    assert leaky.median_error_cm <= result.median_error_cm
