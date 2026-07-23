"""Tests for the QligFEP.IO module.

Tests the run_qprep() function and error handling.
"""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from QligFEP.IO import read_slurm_diagnostics


class TestRunQprep:
    """Tests for the run_qprep() function."""

    def test_run_qprep_calls_subprocess_with_correct_args(self, tmp_path):
        """run_qprep should call subprocess.run with correct command."""
        # Import here to test that module loads correctly
        from QligFEP.IO import run_qprep

        # Create a fake qprep.out file that will be created by the mock
        output_file = tmp_path / "qprep.out"
        input_file = tmp_path / "qprep.inp"
        input_file.write_text("test input")

        # Create output file with no errors (empty output)
        output_file.write_text("")

        with patch("QligFEP.IO.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            run_qprep(
                qprep_path="/path/to/qprep",
                input_file=str(input_file),
                output_file=str(output_file),
                ff_name="AMBER14sb",
            )

            # Verify subprocess.run was called
            mock_run.assert_called_once()
            call_args = mock_run.call_args
            cmd = call_args[0][0]
            assert "/path/to/qprep" in cmd
            assert str(input_file) in cmd
            assert str(output_file) in cmd
            assert call_args[1]["shell"] is True

    def test_run_qprep_raises_on_general_error(self, tmp_path):
        """run_qprep should raise QprepError on general qprep errors."""
        from QligFEP.IO import QprepError, run_qprep

        output_file = tmp_path / "qprep.out"
        input_file = tmp_path / "qprep.inp"
        input_file.write_text("test input")

        # Create output file with error
        output_file.write_text("ERROR: Something went wrong\n")

        with patch("QligFEP.IO.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            with pytest.raises(QprepError):
                run_qprep(
                    qprep_path="/path/to/qprep",
                    input_file=str(input_file),
                    output_file=str(output_file),
                    ff_name="AMBER14sb",
                )

    def test_run_qprep_raises_on_atom_lib_missing(self, tmp_path):
        """run_qprep should raise QprepAtomLibMissingError on missing atom/lib errors."""
        from QligFEP.IO import QprepAtomLibMissingError, run_qprep

        output_file = tmp_path / "qprep.out"
        input_file = tmp_path / "qprep.inp"
        input_file.write_text("test input")

        # Create output file with atom/lib missing error
        output_file.write_text(">>> Atom CA in residue no.   42 not found in library entry for ALA\n")

        with patch("QligFEP.IO.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            with pytest.raises(QprepAtomLibMissingError):
                run_qprep(
                    qprep_path="/path/to/qprep",
                    input_file=str(input_file),
                    output_file=str(output_file),
                    ff_name="AMBER14sb",
                )

    def test_run_qprep_returns_completed_process(self, tmp_path):
        """run_qprep should return the CompletedProcess result."""
        from QligFEP.IO import run_qprep

        output_file = tmp_path / "qprep.out"
        input_file = tmp_path / "qprep.inp"
        input_file.write_text("test input")
        output_file.write_text("")  # No errors

        mock_result = MagicMock(spec=subprocess.CompletedProcess)
        mock_result.returncode = 0

        with patch("QligFEP.IO.subprocess.run", return_value=mock_result) as mock_run:
            result = run_qprep(
                qprep_path="/path/to/qprep",
                input_file=str(input_file),
                output_file=str(output_file),
                ff_name="AMBER14sb",
            )

            assert result is mock_result

    def test_run_qprep_uses_default_filenames(self, tmp_path):
        """run_qprep should use qprep.inp and qprep.out as defaults."""
        from QligFEP.IO import run_qprep

        # Change to tmp_path directory for this test
        import os

        original_cwd = os.getcwd()
        os.chdir(tmp_path)

        try:
            # Create default input and output files
            input_file = tmp_path / "qprep.inp"
            output_file = tmp_path / "qprep.out"
            input_file.write_text("test input")
            output_file.write_text("")

            with patch("QligFEP.IO.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0)

                # Call with just qprep_path - should use defaults
                run_qprep(qprep_path="/path/to/qprep")

                call_args = mock_run.call_args
                cmd = call_args[0][0]
                assert "qprep.inp" in cmd
                assert "qprep.out" in cmd
        finally:
            os.chdir(original_cwd)


class TestQprepErrorCheck:
    """Tests for the qprep_error_check() function."""

    def test_no_error_raises_nothing(self, tmp_path):
        """qprep_error_check should not raise when no errors in output."""
        from QligFEP.IO import qprep_error_check

        output_file = tmp_path / "qprep.out"
        output_file.write_text("Normal output\nNo errors here\n")

        # Should not raise
        qprep_error_check(output_file, "AMBER14sb")

    def test_general_error_raises_qprep_error(self, tmp_path):
        """qprep_error_check should raise QprepError on ERROR: lines."""
        from QligFEP.IO import QprepError, qprep_error_check

        output_file = tmp_path / "qprep.out"
        output_file.write_text("Some output\nERROR: Something failed\nMore output\n")

        with pytest.raises(QprepError):
            qprep_error_check(output_file, "AMBER14sb")

    def test_atom_missing_raises_atom_lib_error(self, tmp_path):
        """qprep_error_check should raise QprepAtomLibMissingError on atom/lib errors."""
        from QligFEP.IO import QprepAtomLibMissingError, qprep_error_check

        output_file = tmp_path / "qprep.out"
        output_file.write_text(
            "Some output\n>>> Atom CA in residue no.   42 not found in library entry for ALA\n"
        )

        with pytest.raises(QprepAtomLibMissingError):
            qprep_error_check(output_file, "AMBER14sb")

    def test_heavy_atom_missing_raises_atom_lib_error(self, tmp_path):
        """qprep_error_check should raise QprepAtomLibMissingError on heavy atom missing."""
        from QligFEP.IO import QprepAtomLibMissingError, qprep_error_check

        output_file = tmp_path / "qprep.out"
        output_file.write_text(">>> Heavy atom CA missing in residue    42\n")

        with pytest.raises(QprepAtomLibMissingError):
            qprep_error_check(output_file, "AMBER14sb")


class TestExceptionClasses:
    """Tests for exception classes."""

    def test_qprep_error_exists(self):
        """QprepError exception class should exist and be importable."""
        from QligFEP.IO import QprepError

        assert issubclass(QprepError, Exception)

    def test_qprep_atom_lib_missing_error_exists(self):
        """QprepAtomLibMissingError exception class should exist and be importable."""
        from QligFEP.IO import QprepAtomLibMissingError

        assert issubclass(QprepAtomLibMissingError, Exception)


class TestParseQprepTotalCharge:
    """Tests for parse_qprep_total_charge() function."""

    def test_parse_real_qprep_out(self):
        """Parse the real Tyk2 qprep.out and expect total charge = 7."""
        from QligFEP.IO import parse_qprep_total_charge

        qprep_out = Path(__file__).parents[2] / "tutorials" / "Tyk2" / "setupFEP" / "FEP_ejm_31_ejm_42" / "inputfiles" / "qprep.out"
        assert qprep_out.exists(), f"Test fixture not found: {qprep_out}"
        charge = parse_qprep_total_charge(qprep_out)
        assert charge == 7

    def test_parse_missing_charge_line_raises(self, tmp_path):
        """Raise ValueError when the charge line is not found."""
        from QligFEP.IO import parse_qprep_total_charge

        fake_out = tmp_path / "qprep.out"
        fake_out.write_text("some random output\nno charge info here\n")
        with pytest.raises(ValueError, match="total charge"):
            parse_qprep_total_charge(fake_out)

    def test_parse_negative_charge(self, tmp_path):
        """Parse a negative total charge correctly."""
        from QligFEP.IO import parse_qprep_total_charge

        fake_out = tmp_path / "qprep.out"
        fake_out.write_text("total charge of not excluded:  -3.00\n")
        charge = parse_qprep_total_charge(fake_out)
        assert charge == -3

    def test_parse_zero_charge(self, tmp_path):
        """Parse zero total charge correctly."""
        from QligFEP.IO import parse_qprep_total_charge

        fake_out = tmp_path / "qprep.out"
        fake_out.write_text("total charge of not excluded:   0.00\n")
        charge = parse_qprep_total_charge(fake_out)
        assert charge == 0


def _slurm(path, seed, replicate, runtime="1h:2m:3s", failure=None, exit_status="0"):
    """Write a run.sh-style slurm*.out footer (the block read_slurm_diagnostics parses).

    ``exit_status=None`` omits the line, reproducing logs written before the run scripts
    reported their own exit code.
    """
    lines = [failure] if failure else []
    lines += [
        f"#    Runtime: {runtime}",
        f"#    Random seed: {seed}",
        f"#    Replicate Number: {replicate}",
    ]
    if exit_status is not None:
        lines.append(f"#    Exit status: {exit_status}")
    path.write_text("\n".join(lines) + "\n")


class TestReadSlurmDiagnostics:
    """Tests for the shared slurm*.out footer parser used by both FEP analyzers."""

    def test_reads_runtime_seed_replicate_status(self, tmp_path):
        f = tmp_path / "slurm.run1.node.12345.out"
        _slurm(f, seed=42, replicate=1)
        info = read_slurm_diagnostics(f)
        assert (info.runtime, info.seed, info.replicate, info.status) == ("1h:2m:3s", "42", "1", "SUCCESS")

    def test_replicate_comes_from_footer_not_filename(self, tmp_path):
        # Multi-temperature array job: the slurm filename %a (5) differs from the true replicate
        # (2). The footer's "Replicate Number" must win, so the run is not mislabeled -- the bug
        # the old filename-regex approach had.
        f = tmp_path / "slurm.run5.node.12345.out"
        _slurm(f, seed=99, replicate=2)
        info = read_slurm_diagnostics(f)
        assert info.replicate == "2"
        assert info.seed == "99"

    def test_seed_from_footer_not_parameters_line(self, tmp_path):
        # A Parameters line carrying a trailing field (as run_neq.sh writes) must not feed the
        # seed; the "#    Random seed:" footer is the single source of truth.
        f = tmp_path / "slurm.run1.node.1.out"
        f.write_text(
            "Parameters T=298, replicate=1, seed=42, neq_reps=5\n"
            "#    Runtime: 0h:5m:0s\n"
            "#    Random seed: 777\n"
            "#    Replicate Number: 1\n"
        )
        info = read_slurm_diagnostics(f)
        assert info.seed == "777"
        assert info.runtime == "0h:5m:0s"
        assert info.status == "SUCCESS"

    @pytest.mark.parametrize(
        "marker,expected",
        [
            ("DUE TO TIME LIMIT", "TIMEOUT"),
            ("slurmstepd: error: job CANCELLED AT ...", "CANCELLED"),
            ("Out Of Memory", "OOM"),
            ("qdyn terminated abnormally", "CRASHED"),
            ("There are not enough slots available in the system", "MPI_LAUNCH_FAILED"),
            ("mpirun was unable to launch the specified application", "MPI_LAUNCH_FAILED"),
            ("A request was made to bind to that would result in binding more", "MPI_LAUNCH_FAILED"),
        ],
    )
    def test_failure_markers_map_to_status(self, tmp_path, marker, expected):
        f = tmp_path / "slurm.run1.node.1.out"
        _slurm(f, seed=1, replicate=1, failure=marker)
        assert read_slurm_diagnostics(f).status == expected

    def test_nonzero_exit_status_fails_a_run_with_no_known_marker(self, tmp_path):
        # An mpirun launch failure leaves a short but well-formed log: the footer parses, and
        # nothing in the body matches a marker. The reported exit code is what catches it.
        f = tmp_path / "slurm.run1.node.1.out"
        _slurm(f, seed=1, replicate=1, exit_status="1")
        assert read_slurm_diagnostics(f).status == "FAILED"

    def test_known_marker_outranks_the_exit_status(self, tmp_path):
        # Both signals present: keep the specific cause rather than flattening it to FAILED.
        f = tmp_path / "slurm.run1.node.1.out"
        _slurm(f, seed=1, replicate=1, failure="DUE TO TIME LIMIT", exit_status="1")
        assert read_slurm_diagnostics(f).status == "TIMEOUT"

    def test_absent_exit_status_still_reads_as_success(self, tmp_path):
        # Logs predating the Exit status footer line must not be reclassified as failures.
        f = tmp_path / "slurm.run1.node.1.out"
        _slurm(f, seed=1, replicate=1, exit_status=None)
        assert read_slurm_diagnostics(f).status == "SUCCESS"

    def test_missing_footer_defaults_gracefully(self, tmp_path):
        # A log killed before writing its footer yields empty/None fields and SUCCESS, so callers
        # degrade gracefully instead of crashing.
        f = tmp_path / "slurm.run1.node.1.out"
        f.write_text("Running job ...\nsome mid-run output\n")
        assert read_slurm_diagnostics(f) == ("", None, None, "SUCCESS")


