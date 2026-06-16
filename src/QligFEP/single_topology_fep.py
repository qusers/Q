"""Single-topology single-Hamiltonian FEP builder.

A single-topology hybrid represents the shared core of two congeneric ligands once
(those atoms morph A->B) and only grows/shrinks the unique atoms, so it inherits the
small-difference convergence of a multistate edge while keeping the soft-core acting
on the appearing/vanishing atoms. `SingleTopologyFEP` subclasses `FEP` and overrides
only the topology-specific build steps; the shared machinery (lambda schedule, prm
merge, qprep run, counter-ion handling, ...) is reused from the base.

The hybrid = all of ligand A (shared core + vanishing atoms) plus ligand B's unique
(appearing) atoms. Shared atoms morph A->B (charge interpolates; LJ type morphs to the
X-prefixed ligand-B type); vanishing atoms go real->DUM; appearing atoms go DUM->real.
Ligand B's atom types clash with ligand A's (both use c1,c2,...), so every ligand-B
type is X-prefixed exactly as the dual-topology builder does -- only here the shared
atoms are collapsed to one copy.
"""
import os
import re
import shutil
import stat
from pathlib import Path

from rdkit import Chem

from .IO import replace, run_qprep
from .logger import logger
from .pdb_utils import pdb_parse_out, read_pdb_to_dataframe
from .qligfep import FEP
from .settings.settings import CLUSTER_DICT, CONFIGS
from .single_topology import map_single_topology
from .templates import get_equilibration_configs, get_production_config, render_md_input
from .templates.sections import format_sequence_restraint

#: Residue name of the single hybrid ligand (QligFEP convention, kept literal).
HYBRID_RESNAME = "{LIG}"

#: Solvent spec per leg. The protein leg merges the pre-equilibrated water.pdb;
#: the water leg lets qprep build a fresh HOH grid around the lone ligand.
SOLVENT_SPEC = {"protein": "4 water.pdb", "water": "1 HOH", "vacuum": "1 HOH"}

#: qprep "Missing <kind> type" -> prm section + leading-token count for that term.
_BOUNDARY_SECTION = {"bond": "[bonds]", "angle": "[angles]", "torsion": "[torsions]", "improper": "[impropers]"}
_BOUNDARY_NLEAD = {"bond": 2, "angle": 3, "torsion": 4, "improper": 4}


def _parse_lib(path):
    """Return atoms=[(name, type, charge)], bonds=[(n1, n2)], impropers=[(n1..n4)]."""
    atoms, bonds, impropers, block = [], [], [], None
    for line in Path(path).read_text().splitlines():
        s = line.split()
        if not s:
            continue
        if s[0].startswith("["):
            block = s[0]
            continue
        if block == "[atoms]" and s[0].isdigit():
            atoms.append((s[1], s[2], float(s[3])))
        elif block == "[bonds]":
            bonds.append((s[0], s[1]))
        elif block == "[impropers]":
            impropers.append(tuple(s[:4]))
    return atoms, bonds, impropers


def _parse_pdb_coords(path):
    """name -> (x, y, z) from ATOM/HETATM records."""
    coords = {}
    for line in Path(path).read_text().splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            coords[line[12:16].strip()] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            )
    return coords


def _element(name):
    """Element symbol = the alphabetic part of an atom name (C1 -> C, Cl11 -> Cl)."""
    return "".join(ch for ch in name if ch.isalpha())


class SingleTopologyFEP(FEP):
    """Build a single-topology hybrid for single-Hamiltonian FEP of two ligands."""

    def read_files(self):
        """Map the two ligands and build the hybrid atom table.

        Returns the same nested shape the orchestration expects:
        ([_, type_replacements], [change_charges, change_atoms], [n_hybrid, 0]).
        The X-prefix `type_replacements` feeds the inherited change_prm; the
        per-atom change_charges/change_atoms and the hybrid table feed
        write_FEP_file. The whole hybrid counts as one restraint group, so the
        second "ligand size" is 0.
        """
        atoms_a, bonds_a, impropers_a = _parse_lib(self.lig1 + ".lib")
        atoms_b, bonds_b, impropers_b = _parse_lib(self.lig2 + ".lib")
        coords_a = _parse_pdb_coords(self.lig1 + ".pdb")
        coords_b = _parse_pdb_coords(self.lig2 + ".pdb")

        mol_a = Chem.MolFromMolFile(self.lig1 + ".sdf", removeHs=False)
        mol_b = Chem.MolFromMolFile(self.lig2 + ".sdf", removeHs=False)
        if len(atoms_a) != mol_a.GetNumAtoms() or len(atoms_b) != mol_b.GetNumAtoms():
            raise ValueError("lib atom count does not match sdf; lib order must equal sdf order")

        self.st_map = map_single_topology(mol_a, mol_b, overlap_thresh=0.6)
        a_to_b = self.st_map.a_to_b
        shared_a = set(a_to_b)
        shared_b = set(a_to_b.values())

        # ligand-B atom index -> hybrid name: shared keep their A-partner's name,
        # appearing atoms are X-prefixed.
        name_b = {j: "X" + atoms_b[j][0] for j in self.st_map.appear}
        for i_a, j_b in a_to_b.items():
            name_b[j_b] = atoms_a[i_a][0]

        hyb = []
        for i, (name, ljtype, charge) in enumerate(atoms_a):
            if i in shared_a:
                j = a_to_b[i]
                state_b = ("X" + atoms_b[j][1], atoms_b[j][2])
                role = "shared"
            else:
                state_b = ("DUM", 0.0)
                role = "vanish"
            hyb.append({
                "name": name, "type": ljtype, "charge": charge, "xyz": coords_a[name],
                "role": role, "sA": (ljtype, charge), "sB": state_b, "elem": _element(name),
            })
        for j in self.st_map.appear:
            state_b = ("X" + atoms_b[j][1], atoms_b[j][2])
            hyb.append({
                "name": name_b[j], "type": state_b[0], "charge": atoms_b[j][2],
                "xyz": coords_b[atoms_b[j][0]], "role": "appear",
                "sA": ("DUM", 0.0), "sB": state_b, "elem": _element(atoms_b[j][0]),
            })

        self.hyb = hyb
        self.hyb_index = {h["name"]: k + 1 for k, h in enumerate(hyb)}
        self.bonds = self._merge_terms(bonds_a, bonds_b, atoms_b, name_b, shared_b)
        self.impropers = self._merge_terms(impropers_a, impropers_b, atoms_b, name_b, shared_b)

        type_replacements = {ljtype: "X" + ljtype for (_, ljtype, _) in atoms_b}
        change_charges = [[k + 1, h["sA"][1], h["sB"][1]] for k, h in enumerate(hyb)]
        change_atoms = [[k + 1, h["sA"][0], h["sB"][0]] for k, h in enumerate(hyb)]

        self.charge_lig1 = round(sum(h["sA"][1] for h in hyb))
        self.charge_lig2 = round(sum(h["sB"][1] for h in hyb))
        self.same_charge = self.charge_lig1 == self.charge_lig2

        return ([type_replacements, type_replacements],
                [change_charges, change_atoms],
                [len(hyb), 0])

    @staticmethod
    def _merge_terms(terms_a, terms_b, atoms_b, name_b, shared_b):
        """Bonds/impropers of ligand A as-is, plus ligand-B terms that touch at
        least one appearing atom (mapped to hybrid names)."""
        b_names = [a[0] for a in atoms_b]
        b_name_to_hyb = {atoms_b[j][0]: name_b[j] for j in range(len(atoms_b))}
        merged = list(terms_a)
        for term in terms_b:
            indices = [b_names.index(n) for n in term]
            if any(j not in shared_b for j in indices):
                merged.append(tuple(b_name_to_hyb[n] for n in term))
        return merged

    def change_lib(self, replacements, writedir):
        """Write the single hybrid library (replaces the dual two-library scheme)."""
        with open(writedir + "/singletop.lib", "w") as f:
            f.write(HYBRID_RESNAME + "\n\n[info]\n SYBYLtype RESIDUE\n\n[atoms]\n")
            for k, h in enumerate(self.hyb, 1):
                f.write(f"{k:5d}   {h['name']:<9s} {h['type']:<12s} {h['charge']:8.4f}\n")
            f.write("\n[bonds]\n")
            for x, y in self.bonds:
                f.write(f"{x:<9s} {y:<9s}\n")
            f.write("\n[impropers]\n")
            for improper in self.impropers:
                f.write("  ".join(f"{n:<8s}" for n in improper) + "\n")

    def _fep_atom_offset(self, writedir):
        """FEP-file topology offset. The water leg's only solute is the hybrid, so
        topology index == FEP index (offset 0). In the protein leg the hybrid is the
        last residue, so the offset is every preceding atom (protein + solvent)."""
        if self.system == "water":
            return 0
        total = read_pdb_to_dataframe(Path(writedir) / "top_p.pdb").shape[0]
        return total - len(self.hyb)

    def write_FEP_file(
        self, change_charges, change_atoms, FEP_vdw, writedir, lig_size1, lig_size2,
        softcore_method="gapsys", single_hamiltonian=True,
    ):
        """Write the single-Hamiltonian FEP file: per-atom charge/type morph with
        DUM-conditional soft-core (0 on real ends, 20 on a ghost end)."""
        offset = self._fep_atom_offset(writedir)
        with open(writedir + "/FEP1.fep", "w") as f:
            f.write(f"!info: {self.lig1} --> {self.lig2} (single topology)\n")
            f.write("[FEP]\nstates 2\nsoftcore_use_max_potential on\n")
            f.write(f"softcore_method {softcore_method}\n")
            if single_hamiltonian:
                f.write("single_hamiltonian on\n")

            f.write("\n[atoms]\n")
            for k in range(1, len(self.hyb) + 1):
                f.write(f"{k:<5d}{k + offset}\n")

            f.write("\n[change_charges]\n")
            for idx, charge_a, charge_b in change_charges:
                f.write(f"{idx:<5d}{charge_a:9.4f}{charge_b:9.4f}\n")

            f.write("\n[atom_types]\n")
            f.write("DUM       0.0000    0.0000    0    0    0.0000    0.0000    1.0080\n")
            for line in FEP_vdw:
                f.write(line + "\n")

            f.write("\n[change_atoms]\n")
            for idx, type_a, type_b in change_atoms:
                f.write(f"{idx:<5d}{type_a:<9s}{type_b:<9s}\n")

            f.write("\n[softcore]\n")
            for k, h in enumerate(self.hyb, 1):
                sc_a = 20 if h["sA"][0] == "DUM" else 0
                sc_b = 20 if h["sB"][0] == "DUM" else 0
                f.write(f"{k:<5d}{sc_a:<6d}{sc_b}\n")

    def merge_pdbs(self, writedir):
        """Write the merged PDB (self.pdb_fname): the prepared protein (protein leg
        only) followed by the single hybrid ligand as one LIG residue placed last.
        Single topology has one solute molecule, so unlike dual topology there is no
        LID partner and no coordinate nudge to separate overlapping copies."""
        lig_resnr = self.residueoffset + 1
        serial = self.atomoffset
        lines = []
        if self.system == "protein":
            prot = Path("protein.pdb").read_text()
            lines.append(prot if prot.endswith("\n") else prot + "\n")
        for h in self.hyb:
            serial += 1
            x, y, z = h["xyz"]
            entry = ["HETATM", serial, h["name"], "", "LIG", "", lig_resnr, "",
                     float(x), float(y), float(z), 1.0, 0.0, h["elem"].upper(), ""]
            lines.append(pdb_parse_out(entry) + "\n")
        Path(writedir, self.pdb_fname).write_text("".join(lines))

    def set_restraints(self, writedir, restraint_method, strict_check: bool = True) -> list:
        """Single topology is one molecule, so there is no LIG<->LID distance tether.
        The hybrid is instead held by a sequence (position) restraint applied per
        stage in write_md_files; no distance-restraint pairs are returned."""
        return []

    def write_qprep(self, writedir):
        """Write qprep.inp for the single hybrid: one force-field library + the single
        singletop.lib + the parameter file. The shared-core/appearing-R-group junction
        introduces bonded terms the merged prm lacks; self.qprep() adds them iteratively
        into singletop_full.prm, which this file references."""
        center, density, cysbond_str = self._prepare_qprep_geometry(writedir)
        solvate = self.system != "vacuum"
        solvent = SOLVENT_SPEC[self.system]
        lines = [
            f"rl {self.lib_file}",
            "rl singletop.lib",
            "rprm singletop_full.prm",
            f"rp {self.pdb_fname}",
            "set solvent_pack 2.3",
            f"set solute_density {density:.5f}",
            f"boundary 1 {center} {self.sphereradius}",
        ]
        lines.append(
            f"solvate {center} {self.sphereradius} {solvent}"
            if solvate
            else f"!solvate {center} {self.sphereradius} {solvent}"
        )
        if cysbond_str:
            lines.append(cysbond_str.rstrip())
        lines += [f"maketop singletop_{self.system}", "writetop singletop.top", "wp top_p.pdb y", "q"]
        Path(writedir, "qprep.inp").write_text("\n".join(lines) + "\n")

    def qprep(self, writedir, max_iterations: int = 8):
        """Run qprep, adding the missing junction boundary parameters between passes
        (qprep reveals bonds, then angles, then torsions/impropers as each level builds)
        until none remain or max_iterations is reached. The merged prm is copied to
        singletop_full.prm, which accumulates the boundary parameters."""
        start_cwd = os.getcwd()
        os.chdir(writedir)
        try:
            qprep_path = CLUSTER_DICT[self.cluster]["QPREP"]
            merged_prm = f"{self.FF}_{self.lig1}_{self.lig2}_merged.prm"
            full_prm = "singletop_full.prm"
            shutil.copy(merged_prm, full_prm)
            seen = set()
            it = 0
            for it in range(1, max_iterations + 1):
                run_qprep(qprep_path, "qprep.inp", "qprep.out", self.FF)
                missing = self._parse_missing_params("qprep.out")
                new = [m for m in missing if m not in seen]
                if not new:
                    break
                seen.update(new)
                self._add_boundary_params(full_prm, new)
            remaining = self._parse_missing_params("qprep.out")
            if remaining:
                logger.warning(
                    f"{len(remaining)} missing parameter(s) remain after {it} qprep pass(es) for "
                    f"{self.lig1}->{self.lig2} ({self.system}); first: {remaining[0]}"
                )
            if not Path("singletop.top").exists():
                logger.error(
                    f"qprep did not produce singletop.top for {self.lig1}->{self.lig2} ({self.system})"
                )
        finally:
            os.chdir(start_cwd)

    @staticmethod
    def _parse_missing_params(out_path):
        """Return [(kind, (atom_types,...)), ...] from qprep's 'Missing ... type' lines."""
        pat = re.compile(r"Missing (bond|angle|torsion|improper) type for atom types (.+)")
        missing = []
        for line in Path(out_path).read_text().splitlines():
            m = pat.search(line)
            if m:
                missing.append((m.group(1), tuple(m.group(2).split())))
        return missing

    @staticmethod
    def _zero_boundary_line(kind, types):
        """A vanishing<->appearing cross-term couples states that never coexist, so it
        gets a zero force constant (with a placeholder geometry Q still parses)."""
        tail = {
            "bond": "      0.0     1.50",
            "angle": "      0.0     120.0",
            "torsion": "      0.000     1.000     0.000    1",
            "improper": "      0.0     0.0",
        }[kind]
        return "   ".join(types) + tail

    def _add_boundary_params(self, full_prm_path, missing):
        """Append the boundary parameters qprep reports missing at the shared-core /
        appearing-R-group junction. A term mixing a shared-core type with an appearing
        type is copied from its all-X (pure ligand-B) equivalent, mapping shared types
        to their morph partner; a term coupling a vanishing type to an appearing type is
        zeroed. Extras are flushed at the end of their section so qprep re-reads them."""
        shared_to_X = {
            h["sA"][0]: h["sB"][0] for h in self.hyb if h["sA"][0] != "DUM" and h["sB"][0] != "DUM"
        }
        vanish = {h["sA"][0] for h in self.hyb if h["sB"][0] == "DUM" and h["sA"][0] != "DUM"}

        text = Path(full_prm_path).read_text()
        section_lines, cur = {}, None
        for line in text.splitlines():
            s = line.strip()
            if s.startswith("[") and s.endswith("]"):
                cur = s
                section_lines.setdefault(cur, [])
            elif cur is not None and s and not s.startswith("!"):
                section_lines[cur].append(line)

        extras = {sec: [] for sec in _BOUNDARY_SECTION.values()}
        for kind, types in missing:
            sec, n = _BOUNDARY_SECTION[kind], _BOUNDARY_NLEAD[kind]
            if any(t in vanish for t in types):
                extras[sec].append(self._zero_boundary_line(kind, types))
                continue
            allx = [shared_to_X.get(t, t) for t in types]
            hit = None
            for ln in section_lines.get(sec, []):
                tok = ln.split()
                if tok[:n] == allx or tok[:n] == allx[::-1]:
                    ordered = list(types) if tok[:n] == allx else list(types)[::-1]
                    hit = "   ".join(ordered) + "   " + "   ".join(tok[n:])
                    break
            if hit:
                extras[sec].append(hit)
            else:
                logger.warning(f"No all-X match for missing {kind} {types} (looked up {allx})")

        out, cur = [], None
        for line in text.splitlines():
            s = line.strip()
            if s.startswith("[") and s.endswith("]"):
                if cur in extras and extras[cur]:
                    out.extend(extras[cur])
                    extras[cur] = []
                cur = s
            out.append(line)
        if cur in extras and extras[cur]:
            out.extend(extras[cur])
        Path(full_prm_path).write_text("\n".join(out) + "\n")

    def write_md_files(
        self, lambdas, writedir, lig_size1, lig_size2, overlapping_atoms, wall_restraints: str = ""
    ):
        """Write equilibration (eq1-eq5) and the single-Hamiltonian production chain.

        There is no LIG<->LID tether; the hybrid is held by a sequence (position)
        restraint on its atom range in every stage. The production windows form one
        continuous chain restarting from the equilibrated end, each carrying a
        [foreign_lambdas] schedule so the engine writes the true U(lambda') energies
        (the .en.fl file) the single-Hamiltonian BAR/MBAR analysis consumes.

        Leg-specific protocol (validated): the protein leg equilibrates at lambda_c=0.5
        and the chain runs lambda_c 0->1; the water leg equilibrates at lambda_c=1 (the
        grown end, so dense water has no DUM cavity to clash into) and the chain runs
        1->0, restarting from that equilibrated end.

        Returns:
            list: [eq_files, chain_files, []] -- chain_files in execution order.
        """
        offset = self._fep_atom_offset(writedir)
        first, last = offset + 1, offset + len(self.hyb)
        n_windows = len(lambdas)
        sched = [round(i / (n_windows - 1), 4) for i in range(n_windows)]
        foreign = " ".join(f"{lc:.4f}" for lc in sched)
        eq_lc = 1.0 if self.system == "water" else 0.5

        eq_files = []
        eq_configs = get_equilibration_configs(self.timestep, int(self.sphereradius))
        for i, eq_config in enumerate(eq_configs):
            eq_config.params.topology = "singletop.top"
            # eq1-eq4 use the graded restraint force; eq5 keeps the hybrid pinned
            # (force 10) since no tether takes over for production.
            force = eq_config.sequence_restraint_force if i < 4 else 10.0
            content = render_md_input(
                params=eq_config.params,
                lambda1=f"{1 - eq_lc:.3f}",
                lambda2=f"{eq_lc:.3f}",
                trajectory_file=f"{eq_config.name}.dcd",
                final_file=f"{eq_config.name}.re",
                restart_file=f"eq{i}.re" if i > 0 else None,
                sequence_restraints=format_sequence_restraint(first, last, force),
                is_eq1=(i == 0),
            )
            (Path(writedir) / f"{eq_config.name}.inp").write_text(content)
            eq_files.append(f"{eq_config.name}.inp")

        prod_config = get_production_config(self.timestep, int(self.sphereradius))
        prod_config.params.topology = "singletop.top"
        seq10 = format_sequence_restraint(first, last, 10.0)
        chain_order = list(reversed(sched)) if self.system == "water" else list(sched)

        chain_files = []
        prev = "eq5.re"
        for lc in chain_order:
            tag = f"md_{int(round(lc * 1000)):04d}"
            content = render_md_input(
                params=prod_config.params,
                lambda1=f"{1 - lc:.4f}",
                lambda2=f"{lc:.4f}",
                trajectory_file=f"{tag}.dcd",
                final_file=f"{tag}.re",
                restart_file=prev,
                energy_file=f"{tag}.en",
                sequence_restraints=seq10,
                foreign_lambdas=foreign,
            )
            (Path(writedir) / f"{tag}.inp").write_text(content)
            chain_files.append(f"{tag}.inp")
            prev = f"{tag}.re"

        return [eq_files, chain_files, []]

    def write_runfile(self, writedir, file_list):
        """Write the SLURM run script. Reuses the standard template (replicate/seed/
        temperature handling) but lists eq + the single-Hamiltonian production chain
        in execution order and drops the qfep step: single-Hamiltonian energies are
        non-linear in lambda, so ddG comes from a post-hoc foreign-lambda BAR/MBAR over
        the .en.fl files, not qfep's linear reconstruction."""
        eq_files, chain_files = file_list[0], file_list[1]
        src = CONFIGS["INPUT_DIR"] + "/run.sh"
        tgt = writedir + "/run" + self.cluster + ".sh"
        replacements = dict(CLUSTER_DICT[self.cluster])
        replacements["FEPS"] = "FEP1.fep"
        with open(src) as infile, open(tgt, "w") as outfile:
            for line in infile:
                stripped = line.strip()
                # qfep is invalid for single-Hamiltonian; ddG is computed post-hoc.
                if stripped.startswith("cp $inputfiles/qfep.inp") or stripped.startswith("timeout 3m QFEP"):
                    continue
                if stripped == "#SBATCH --array=1-TOTAL_JOBS":
                    replacements["TOTAL_JOBS"] = str(self.replicates)
                if stripped == "temperatures=(TEMP_VAR)":
                    replacements["TEMP_VAR"] = str(self.temperature)
                if stripped == "seeds=(RANDOM_SEEDS)":
                    replacements["RANDOM_SEEDS"] = " ".join(str(s) for s in self.seeds)
                if stripped == "#SBATCH -A ACCOUNT" and "ACCOUNT" not in replacements:
                    line = ""
                if stripped == "#SBATCH -J JOBNAME":
                    if self.cluster == "SNELLIUS":
                        outfile.write("#SBATCH -p rome\n")
                    elif self.cluster == "DARDEL":
                        outfile.write("#SBATCH -p shared\n")
                    jobname = {"water": "w_", "protein": "p_", "vacuum": "v_"}[self.system]
                    replacements["JOBNAME"] = jobname + self.lig1 + "_" + self.lig2
                outline = replace(line, replacements)
                outfile.write(outline)
                if stripped == "#EQ_FILES":
                    for f in eq_files:
                        base = Path(f).stem
                        outfile.write(
                            f"time mpirun -n $SLURM_NTASKS --bind-to core $qdyn {base}.inp > {base}.log\n"
                        )
                if stripped == "#RUN_FILES":
                    for f in chain_files:
                        base = Path(f).stem
                        outfile.write(
                            f"time mpirun -n $SLURM_NTASKS --bind-to core $qdyn {base}.inp > {base}.log\n"
                        )
                if stripped == "#CLEANUP" and self.to_clean:
                    outfile.write("rm -r " + " ".join("*" + x for x in self.to_clean) + "\n")
        try:
            st = os.stat(tgt)
            os.chmod(tgt, st.st_mode | stat.S_IEXEC)
        except Exception:
            logger.warning(f"Could not change permission for {tgt}")
