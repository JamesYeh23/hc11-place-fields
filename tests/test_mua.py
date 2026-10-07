"""Population-synchrony detection from spikes."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import SYNTH_NAME, requires_data, write_sessinfo
from hc11 import config
from hc11.events import participation
from hc11.io import load_session
from hc11.mua import detect_mua_events, population_rate, state_intervals_for_epoch


def _burst_session(tmp_path, arrays, burst_times, n_cells=20, burst_cells=15,
                   baseline_hz=1.0, seed=0, nrem=((10.0, 90.0),)):
    """A session with Poisson background plus population bursts at known times."""
    rng = np.random.default_rng(seed)
    ids = 102 + np.arange(n_cells)
    duration = 300.0
    times, cells = [], []
    for cid in ids:
        n = rng.poisson(baseline_hz * duration)
        times.append(rng.uniform(0, duration, n))
        cells.append(np.full(n, cid))
    for centre in burst_times:  # 50 ms burst, 6 spikes each from a subset
        for cid in rng.choice(ids, size=burst_cells, replace=False):
            times.append(centre + rng.uniform(-0.025, 0.025, 6))
            cells.append(np.full(6, cid))
    times = np.concatenate(times)
    cells = np.concatenate(cells)
    order = np.argsort(times)

    arrays = dict(arrays)
    arrays.update(
        SpikeTimes=times[order][:, None], SpikeIDs=cells[order].astype(float)[:, None],
        PyrIDs=ids.astype(float)[:, None], IntIDs=np.array([[1.0]]),
        PREEpoch=np.array([[0.0, 100.0]]), MazeEpoch=np.array([[100.0, 160.0]]),
        POSTEpoch=np.array([[160.0, duration]]), sessDuration=np.array([[duration]]),
        NREM=np.array(nrem, dtype=float),
        Wake=np.array([[100.0, 160.0]]), Drowsy=np.array([[0.0, 1.0]]),
        Intermediate=np.array([[1.0, 2.0]]), REM=np.array([[2.0, 3.0]]),
    )
    return load_session(write_sessinfo(tmp_path / f"{SYNTH_NAME}_sessInfo.mat", arrays))


class TestPopulationRate:
    def test_rate_matches_a_known_constant(self, tmp_path, synthetic_arrays):
        session = _burst_session(tmp_path, synthetic_arrays, [], n_cells=20, baseline_hz=2.0)
        _, rate = population_rate(session, np.array([10.0, 90.0]), session.pyr_ids,
                                  bin_s=0.001, sigma_s=0.05)
        assert np.mean(rate) == pytest.approx(40.0, rel=0.25)  # 20 cells * 2 Hz

    def test_grid_covers_the_interval(self, tmp_path, synthetic_arrays):
        session = _burst_session(tmp_path, synthetic_arrays, [])
        t, rate = population_rate(session, np.array([10.0, 20.0]), session.pyr_ids,
                                  bin_s=0.001, sigma_s=0.01)
        assert t.size == rate.size == 10_000
        assert 10.0 <= t[0] < 10.002 and 19.998 < t[-1] <= 20.0


class TestDetection:
    def test_finds_injected_bursts(self, tmp_path, synthetic_arrays):
        bursts = [20.0, 40.0, 60.0, 80.0]
        session = _burst_session(tmp_path, synthetic_arrays, bursts)
        events, _ = detect_mua_events(session, "PRE", config.load_params())
        assert len(events) >= len(bursts)
        for centre in bursts:
            assert np.min(np.abs(events.peak - centre)) < 0.03, f"missed burst at {centre}"

    def test_injected_bursts_stand_out_from_background_noise(self, tmp_path, synthetic_arrays):
        """What separates a real event from a threshold crossing is participation.

        A smoothed population-rate trace built from Poisson background alone
        still crosses mean + 3 SD from time to time, so event *count* is not by
        itself evidence of synchrony. Those crossings involve a handful of cells
        firing once each; an injected burst involves many. This contrast, not
        the count, is what the LFP cross-check in 6c is meant to confirm.
        """
        bursts = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
        session = _burst_session(tmp_path, synthetic_arrays, bursts, burst_cells=15)
        events, _ = detect_mua_events(session, "PRE", config.load_params())
        part = participation(session, events, session.pyr_ids)
        is_burst = np.min(
            np.abs(events.peak[:, None] - np.array(bursts)[None, :]), axis=1
        ) < 0.05
        assert is_burst.sum() >= 5, "most injected bursts should be found"
        if (~is_burst).any():
            assert part.n_active[is_burst].mean() > 2 * part.n_active[~is_burst].mean()

    def test_only_nrem_is_searched(self, tmp_path, synthetic_arrays):
        """A burst outside the scored state must not be detected."""
        session = _burst_session(
            tmp_path, synthetic_arrays, [20.0, 95.0], nrem=((10.0, 50.0),)
        )
        events, intervals = detect_mua_events(session, "PRE", config.load_params())
        np.testing.assert_allclose(intervals, [[10.0, 50.0]])
        assert np.min(np.abs(events.peak - 20.0)) < 0.03
        assert np.all(events.peak <= 50.0), "nothing may be found outside NREM"
        assert events.state_seconds == pytest.approx(40.0)

    def test_baseline_is_computed_within_the_state(self, tmp_path, synthetic_arrays):
        """Activity outside the searched state must not move the threshold."""
        params = config.load_params()
        quiet = _burst_session(tmp_path, synthetic_arrays, [20.0], nrem=((10.0, 50.0),))
        # Same NREM window, but a storm of extra bursts in MAZE, outside it.
        loud = _burst_session(
            tmp_path, synthetic_arrays, [20.0] + list(np.arange(101.0, 159.0, 0.5)),
            nrem=((10.0, 50.0),),
        )
        a, _ = detect_mua_events(quiet, "PRE", params)
        b, _ = detect_mua_events(loud, "PRE", params)
        assert a.baseline_mean == pytest.approx(b.baseline_mean, rel=1e-9)
        assert a.baseline_sd == pytest.approx(b.baseline_sd, rel=1e-9)

    def test_events_carry_participation(self, tmp_path, synthetic_arrays):
        session = _burst_session(tmp_path, synthetic_arrays, [20.0, 40.0], burst_cells=15)
        events, _ = detect_mua_events(session, "PRE", config.load_params())
        part = participation(session, events, session.pyr_ids)
        assert part.n_cells == 20
        assert part.n_active.max() >= 10, "an injected 15-cell burst should show up"

    def test_durations_respect_the_configured_bounds(self, tmp_path, synthetic_arrays):
        params = config.load_params()
        session = _burst_session(tmp_path, synthetic_arrays, list(np.arange(20.0, 80.0, 5.0)))
        events, _ = detect_mua_events(session, "PRE", params)
        cfg = params["events"]["mua"]
        assert events.duration.min() >= cfg["min_duration_ms"] / 1000 - 1e-9
        assert events.duration.max() <= cfg["max_duration_ms"] / 1000 + 1e-9

    def test_unknown_epoch_raises(self, tmp_path, synthetic_arrays):
        session = _burst_session(tmp_path, synthetic_arrays, [])
        with pytest.raises(ValueError, match="epoch"):
            state_intervals_for_epoch(session, "SIESTA", config.load_params())


@requires_data
@pytest.mark.parametrize("name", ["Achilles_10252013", "Achilles_11012013"])
def test_real_sessions_give_plausible_events(real_sessions, name):
    params = config.load_params()
    session = real_sessions(name)
    events, _ = detect_mua_events(session, "PRE", params)
    part = participation(session, events, session.pyr_ids)
    assert len(events) > 100
    # Durations in the published sharp-wave-ripple range.
    assert 40 < np.median(events.duration) * 1000 < 120
    # Participation around the 10-15% of the pyramidal population the paper describes.
    assert 0.05 < np.median(part.n_active) / len(session.pyr_ids) < 0.30
    # Events never leave the scored state.
    nrem = session.states["NREM"]
    pre = session.epoch("PRE")
    inside = np.zeros(len(events), dtype=bool)
    for lo, hi in nrem[(nrem[:, 1] > pre[0]) & (nrem[:, 0] < pre[1])]:
        inside |= (events.peak >= lo) & (events.peak <= hi)
    assert inside.all()
