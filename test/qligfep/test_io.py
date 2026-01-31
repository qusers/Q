"""Tests for the QligFEP.IO module.

Tests the run_qprep() function and error handling.
"""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


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


