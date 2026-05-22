# ABOUTME: Tests for QligFEP.charge_correction: parsing qdyn .corr logs, pooling replicates
# ABOUTME: per window, and trapezoidal lambda-integration of the charge-correction observables.

import numpy as np
import pytest

from QligFEP.charge_correction import (
    collect_leg_corr,
    integrate_windows,
    read_corr_file,
    reduce_windows,
)


def _write_corr(path, lam, frames, phi):
    """Write a .corr file in qdyn's format: step, lambdas..., u_obs..., phi_cog.

    lam: nstates lambda values; frames: list of u_obs rows (length nstates);
    phi: scalar phi_cog (constant across frames here, as in a real window).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for step, row in enumerate(frames, start=1):
        vals = list(lam) + list(row) + [phi]
        lines.append(f"{step * 10:10d}" + "".join(f" {v: .8E}" for v in vals))
    path.write_text("\n".join(lines) + "\n")
    return path


class TestReadCorrFile:
    def test_splits_lambdas_observables_and_phi(self, tmp_path):
        f = _write_corr(tmp_path / "w.corr", lam=[0.5, 0.5], frames=[[1.0, -2.0], [3.0, -4.0]], phi=-7.0)
        lambdas, u_obs, phi = read_corr_file(f)
        assert lambdas.tolist() == [0.5, 0.5]
        np.testing.assert_allclose(u_obs, [[1.0, -2.0], [3.0, -4.0]])
        np.testing.assert_allclose(phi, [-7.0, -7.0])

    def test_empty_file_raises(self, tmp_path):
        f = tmp_path / "empty.corr"
        f.write_text("")
        with pytest.raises(ValueError):
            read_corr_file(f)


class TestReduceWindows:
    def test_pools_replicates_by_window_basename(self, tmp_path):
        f1 = _write_corr(tmp_path / "1" / "md_0000_1000.corr", lam=[0.0, 1.0], frames=[[0.0, 2.0]], phi=-5.0)
        f2 = _write_corr(tmp_path / "2" / "md_0000_1000.corr", lam=[0.0, 1.0], frames=[[0.0, 4.0]], phi=-9.0)
        windows = reduce_windows([f1, f2])
        assert len(windows) == 1
        w = windows[0]
        assert w["lambda2"] == pytest.approx(1.0)
        assert w["du"] == pytest.approx(3.0)   # mean of (2-0) and (4-0)
        assert w["phi"] == pytest.approx(-7.0)  # mean of -5 and -9
        assert w["n_frames"] == 2

    def test_windows_sorted_by_lambda2(self, tmp_path):
        f_hi = _write_corr(tmp_path / "a" / "md_0000_1000.corr", lam=[0.0, 1.0], frames=[[0.0, 1.0]], phi=-1.0)
        f_lo = _write_corr(tmp_path / "a" / "md_1000_0000.corr", lam=[1.0, 0.0], frames=[[0.0, 1.0]], phi=-1.0)
        windows = reduce_windows([f_hi, f_lo])
        assert [w["lambda2"] for w in windows] == [0.0, 1.0]

    def test_rejects_non_two_state(self, tmp_path):
        # 3 states -> 1 + 3 + 3 + 1 = 8 columns; reduce should reject
        path = tmp_path / "w.corr"
        path.write_text(f"        10 {0.2:.8E} {0.3:.8E} {0.5:.8E} {1.0:.8E} {2.0:.8E} {3.0:.8E} {-1.0:.8E}\n")
        with pytest.raises(ValueError):
            reduce_windows([path])


class TestIntegrateWindows:
    def test_integrates_named_key(self):
        windows = [
            {"lambda2": 0.0, "du": 2.0, "phi": -3.0, "n_frames": 1},
            {"lambda2": 1.0, "du": 2.0, "phi": -3.0, "n_frames": 1},
        ]
        assert integrate_windows(windows, key="du") == pytest.approx(2.0)
        assert integrate_windows(windows, key="phi") == pytest.approx(-3.0)

    def test_linear_integrand(self):
        windows = [{"lambda2": x, "du": x, "phi": 0.0, "n_frames": 1} for x in np.linspace(0, 1, 11)]
        assert integrate_windows(windows, key="du") == pytest.approx(0.5, abs=1e-9)


class TestCollectLegCorr:
    def test_globs_and_integrates_both_observables(self, tmp_path):
        for rep in ("1", "2"):
            d = tmp_path / "FEP1" / "298" / rep
            _write_corr(d / "md_1000_0000.corr", lam=[1.0, 0.0], frames=[[0.0, 0.0]], phi=-2.0)
            _write_corr(d / "md_0000_1000.corr", lam=[0.0, 1.0], frames=[[0.0, 4.0]], phi=-6.0)
        run_dirs = [tmp_path / "FEP1" / "298" / "1", tmp_path / "FEP1" / "298" / "2"]
        windows, u_integral, phi_integral = collect_leg_corr(run_dirs)
        assert len(windows) == 2
        # du(0)=0, du(1)=4 -> 2.0 ; phi(0)=-2, phi(1)=-6 -> trapz = -4.0
        assert u_integral == pytest.approx(2.0)
        assert phi_integral == pytest.approx(-4.0)

    def test_returns_none_when_no_corr_files(self, tmp_path):
        d = tmp_path / "FEP1" / "298" / "1"
        d.mkdir(parents=True)
        windows, u_integral, phi_integral = collect_leg_corr([d])
        assert windows is None
        assert u_integral is None
        assert phi_integral is None
