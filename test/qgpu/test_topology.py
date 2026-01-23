"""Tests for the Qgpu topology parser's vdW format handling."""
import os
import sys
import tempfile
import pytest

# Add Qgpu directory to path for imports
QGPU_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/Qgpu'))
sys.path.insert(0, QGPU_DIR)

import topology


# Test topology file paths
GEOMETRIC_TOPOLOGY = os.path.join(
    os.path.dirname(__file__), '../../test/q6/test1/nc20/lig_w.top'
)
ARITHMETIC_TOPOLOGY = os.path.join(
    os.path.dirname(__file__), 'amber_vdw_arithmetic/dualtop.top'
)


class TestVdwRuleParsing:
    """Tests for parsing vdW combination rule from topology files."""

    def test_geometric_vdw_rule_parsed_correctly(self):
        """Test that geometric format (vdw_rule=1) is parsed correctly."""
        reader = topology.Read_Topology(GEOMETRIC_TOPOLOGY)
        data = reader.Q()

        assert data['vdw_rule'] == '1', f"Expected vdw_rule='1', got '{data['vdw_rule']}'"

    def test_arithmetic_vdw_rule_parsed_correctly(self):
        """Test that arithmetic format (vdw_rule=2) is parsed correctly."""
        reader = topology.Read_Topology(ARITHMETIC_TOPOLOGY)
        data = reader.Q()

        assert data['vdw_rule'] == '2', f"Expected vdw_rule='2', got '{data['vdw_rule']}'"

    def test_catypes_populated_for_geometric(self):
        """Test that catypes dict is populated for geometric format."""
        reader = topology.Read_Topology(GEOMETRIC_TOPOLOGY)
        data = reader.Q()

        assert len(data['catypes']) > 0, "catypes should not be empty"
        # Each catype should have 7 values: mass, Aii_normal, Bii_normal,
        # Aii_polar, Bii_polar, Aii_14, Bii_14
        first_key = list(data['catypes'].keys())[0]
        assert len(data['catypes'][first_key]) == 7, "Each catype should have 7 values"

    def test_catypes_populated_for_arithmetic(self):
        """Test that catypes dict is populated for arithmetic format."""
        reader = topology.Read_Topology(ARITHMETIC_TOPOLOGY)
        data = reader.Q()

        assert len(data['catypes']) > 0, "catypes should not be empty"
        first_key = list(data['catypes'].keys())[0]
        assert len(data['catypes'][first_key]) == 7, "Each catype should have 7 values"


class TestVdwRuleValidation:
    """Tests for validation of vdW rule requirements."""

    def test_missing_vdw_rule_raises_error(self):
        """Test that missing vdw_rule causes an error."""
        # Create a minimal topology file without vdw_rule
        with tempfile.NamedTemporaryFile(mode='w', suffix='.top', delete=False) as f:
            f.write("""Q topology file
TITLE      test
DATE       20250123
VERSION    6.01
END        of header
    10    10 = Total no. of atoms, no. of solute atoms. Coordinates: (2*3 per line)
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
    10 = No. of integer atom codes
    1    1    1    1    1    1    1    1    1    1
    0     0 = No. of bonds, no. of solute bonds
    0 = No. of bond codes
    0     0 = No. of angles, no. of solute angles
    0 = No. of angle codes
    0     0 = No. of torsions, no. of solute torsions
    0 = No. of torsion codes
    0     0 = No. of impropers, no. of solute impropers
    0 = No. of improper codes
    10 = No. of atomic charges
  0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
     1     1 = No. of charge groups (total, solute)
     1     1
     1     2     3     4     5     6     7     8     9    10
Electrostatic 1-4 scaling factor and Coulomb constant
  0.5  332.0
Masses:
  12.0
sqrt (Aii) normal:
  1000.0
sqrt (Bii) normal:
  10.0
sqrt (Aii) polar:
  1000.0
sqrt (Bii) polar:
  10.0
sqrt (Aii) 1-4:
  1000.0
sqrt (Bii) 1-4:
  10.0
    0 = No. of type-2 vdW interactions
    0 = No. of 1-4 neighbours
    0 = No. of long 1-4 nbrs
    0 = No. of exclusions
    0 = No. of long exclusions
    1     1 = No. of residues (total, solute)
    1
UNK
Sequence:
    1
    1 = No. of separate molecules
    1 = No. of atom types
C
    0 = No. of SYBYL atom types
    0 = solvent type (0=SPC,1=3-atom,2=general)
  10.0   10.0 = Exclusion and solvent radii
   0.0    0.0    0.0 = Solute center
   0.0    0.0    0.0 = Solvent center
    0 = No. of excluded atoms
""")
            temp_path = f.name

        try:
            reader = topology.Read_Topology(temp_path)
            with pytest.raises(SystemExit):
                reader.Q()
        finally:
            os.unlink(temp_path)

    def test_invalid_vdw_rule_value_raises_error(self):
        """Test that invalid vdw_rule value (not 1 or 2) causes an error."""
        # Create a topology with invalid vdw_rule
        with tempfile.NamedTemporaryFile(mode='w', suffix='.top', delete=False) as f:
            f.write("""Q topology file
TITLE      test
DATE       20250123
VERSION    6.01
END        of header
    10    10 = Total no. of atoms, no. of solute atoms. Coordinates: (2*3 per line)
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
    10 = No. of integer atom codes
    1    1    1    1    1    1    1    1    1    1
    0     0 = No. of bonds, no. of solute bonds
    0 = No. of bond codes
    0     0 = No. of angles, no. of solute angles
    0 = No. of angle codes
    0     0 = No. of torsions, no. of solute torsions
    0 = No. of torsion codes
    0     0 = No. of impropers, no. of solute impropers
    0 = No. of improper codes
    10 = No. of atomic charges
  0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
     1     1 = No. of charge groups (total, solute)
     1     1
     1     2     3     4     5     6     7     8     9    10
       3 = vdW combination rule (1 = Geom. / 2 = Arit.)
Electrostatic 1-4 scaling factor and Coulomb constant
  0.5  332.0
Masses:
  12.0
sqrt (Aii) normal:
  1000.0
sqrt (Bii) normal:
  10.0
sqrt (Aii) polar:
  1000.0
sqrt (Bii) polar:
  10.0
sqrt (Aii) 1-4:
  1000.0
sqrt (Bii) 1-4:
  10.0
    0 = No. of type-2 vdW interactions
    0 = No. of 1-4 neighbours
    0 = No. of long 1-4 nbrs
    0 = No. of exclusions
    0 = No. of long exclusions
    1     1 = No. of residues (total, solute)
    1
UNK
Sequence:
    1
    1 = No. of separate molecules
    1 = No. of atom types
C
    0 = No. of SYBYL atom types
    0 = solvent type (0=SPC,1=3-atom,2=general)
  10.0   10.0 = Exclusion and solvent radii
   0.0    0.0    0.0 = Solute center
   0.0    0.0    0.0 = Solvent center
    0 = No. of excluded atoms
""")
            temp_path = f.name

        try:
            reader = topology.Read_Topology(temp_path)
            with pytest.raises(SystemExit):
                reader.Q()
        finally:
            os.unlink(temp_path)

    def test_vdw_rule_format_mismatch_geometric_with_arithmetic_sections(self):
        """Test that vdw_rule=1 with R*/epsilon sections causes an error."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.top', delete=False) as f:
            f.write("""Q topology file
TITLE      test
DATE       20250123
VERSION    6.01
END        of header
    10    10 = Total no. of atoms, no. of solute atoms. Coordinates: (2*3 per line)
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
    10 = No. of integer atom codes
    1    1    1    1    1    1    1    1    1    1
    0     0 = No. of bonds, no. of solute bonds
    0 = No. of bond codes
    0     0 = No. of angles, no. of solute angles
    0 = No. of angle codes
    0     0 = No. of torsions, no. of solute torsions
    0 = No. of torsion codes
    0     0 = No. of impropers, no. of solute impropers
    0 = No. of improper codes
    10 = No. of atomic charges
  0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
     1     1 = No. of charge groups (total, solute)
     1     1
     1     2     3     4     5     6     7     8     9    10
       1 = vdW combination rule (1 = Geom. / 2 = Arit.)
Electrostatic 1-4 scaling factor and Coulomb constant
  0.5  332.0
Masses:
  12.0
R* normal:
  1.9
epsilon normal:
  0.1
R* polar:
  1.9
epsilon polar:
  0.1
R* 1-4:
  1.9
epsilon 1-4:
  0.1
    0 = No. of type-2 vdW interactions
    0 = No. of 1-4 neighbours
    0 = No. of long 1-4 nbrs
    0 = No. of exclusions
    0 = No. of long exclusions
    1     1 = No. of residues (total, solute)
    1
UNK
Sequence:
    1
    1 = No. of separate molecules
    1 = No. of atom types
C
    0 = No. of SYBYL atom types
    0 = solvent type (0=SPC,1=3-atom,2=general)
  10.0   10.0 = Exclusion and solvent radii
   0.0    0.0    0.0 = Solute center
   0.0    0.0    0.0 = Solvent center
    0 = No. of excluded atoms
""")
            temp_path = f.name

        try:
            reader = topology.Read_Topology(temp_path)
            with pytest.raises(SystemExit):
                reader.Q()
        finally:
            os.unlink(temp_path)

    def test_vdw_rule_format_mismatch_arithmetic_with_geometric_sections(self):
        """Test that vdw_rule=2 with sqrt(Aii)/sqrt(Bii) sections causes an error."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.top', delete=False) as f:
            f.write("""Q topology file
TITLE      test
DATE       20250123
VERSION    6.01
END        of header
    10    10 = Total no. of atoms, no. of solute atoms. Coordinates: (2*3 per line)
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
   0.0   0.0   0.0   0.0   0.0   0.0
    10 = No. of integer atom codes
    1    1    1    1    1    1    1    1    1    1
    0     0 = No. of bonds, no. of solute bonds
    0 = No. of bond codes
    0     0 = No. of angles, no. of solute angles
    0 = No. of angle codes
    0     0 = No. of torsions, no. of solute torsions
    0 = No. of torsion codes
    0     0 = No. of impropers, no. of solute impropers
    0 = No. of improper codes
    10 = No. of atomic charges
  0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
     1     1 = No. of charge groups (total, solute)
     1     1
     1     2     3     4     5     6     7     8     9    10
       2 = vdW combination rule (1 = Geom. / 2 = Arit.)
Electrostatic 1-4 scaling factor and Coulomb constant
  0.5  332.0
Masses:
  12.0
sqrt (Aii) normal:
  1000.0
sqrt (Bii) normal:
  10.0
sqrt (Aii) polar:
  1000.0
sqrt (Bii) polar:
  10.0
sqrt (Aii) 1-4:
  1000.0
sqrt (Bii) 1-4:
  10.0
    0 = No. of type-2 vdW interactions
    0 = No. of 1-4 neighbours
    0 = No. of long 1-4 nbrs
    0 = No. of exclusions
    0 = No. of long exclusions
    1     1 = No. of residues (total, solute)
    1
UNK
Sequence:
    1
    1 = No. of separate molecules
    1 = No. of atom types
C
    0 = No. of SYBYL atom types
    0 = solvent type (0=SPC,1=3-atom,2=general)
  10.0   10.0 = Exclusion and solvent radii
   0.0    0.0    0.0 = Solute center
   0.0    0.0    0.0 = Solvent center
    0 = No. of excluded atoms
""")
            temp_path = f.name

        try:
            reader = topology.Read_Topology(temp_path)
            with pytest.raises(SystemExit):
                reader.Q()
        finally:
            os.unlink(temp_path)


class TestVdwCsvOutput:
    """Tests for vdW rule being correctly written to CSV output."""

    def test_vdw_rule_written_to_topo_csv_geometric(self):
        """Test that vdw_rule is correctly written to topo.csv line 9 for geometric."""
        reader = topology.Read_Topology(GEOMETRIC_TOPOLOGY)
        data = reader.Q()

        with tempfile.TemporaryDirectory() as tmpdir:
            writer = topology.Write_Topology(data)
            writer.CSV(tmpdir + '/')

            with open(os.path.join(tmpdir, 'topo.csv'), 'r') as f:
                lines = f.readlines()

            # Line 9 (0-indexed line 8) should be vdw_rule
            assert lines[8].strip() == '1', f"Expected '1', got '{lines[8].strip()}'"

    def test_vdw_rule_written_to_topo_csv_arithmetic(self):
        """Test that vdw_rule is correctly written to topo.csv line 9 for arithmetic."""
        reader = topology.Read_Topology(ARITHMETIC_TOPOLOGY)
        data = reader.Q()

        with tempfile.TemporaryDirectory() as tmpdir:
            writer = topology.Write_Topology(data)
            writer.CSV(tmpdir + '/')

            with open(os.path.join(tmpdir, 'topo.csv'), 'r') as f:
                lines = f.readlines()

            # Line 9 (0-indexed line 8) should be vdw_rule
            assert lines[8].strip() == '2', f"Expected '2', got '{lines[8].strip()}'"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
