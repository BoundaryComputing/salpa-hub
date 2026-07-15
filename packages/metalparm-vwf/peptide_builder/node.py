"""Peptide Builder Node — assemble a standalone peptide for downstream fusion.

Two input modes:
  - ``sequence``: tleap builds the peptide from a space-separated 3-letter
    sequence (e.g., ``HIS TYR TYR LEU ALA GLU ALA``) under the chosen
    forcefield (ff19SB or ff14SB), with optional ACE/NME caps.
  - ``pdb``: user supplies a peptide PDB directly. Common cases — a
    sequence with mutations, multiple chains, post-translational
    modifications, or a peptide pre-relaxed by an external MD step.
    The PDB must use standard residue/atom names.

Both modes write three artifacts:
  - ``peptide.pdb``  — coordinates
  - ``peptide.lib``  — parameterized AMBER OFF unit (saveoff from tleap)
  - ``peptide.frcmod`` — placeholder (empty for pure ff19SB; kept symmetric
    with the metal side which always ships a paired .lib + .frcmod)

Forwarded data:
  ``output_pdb``, ``output_lib``, ``output_frcmod``, ``forcefield``,
  ``peptide_residues``. ``ep_fragment_fuse_topology`` consumes ``output_lib`` +
  ``output_frcmod`` via the standard 3-tier predecessor pattern.

Pure-string helpers (sequence normalization, tleap script emission, PDB
residue counting) live in ``core.py`` so they're unit-testable without
bocoflow_core or a tleap binary.
"""

from __future__ import annotations

import os
import shutil
import subprocess

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FolderParameter, SelectParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import (
        PEPTIDE_FRCMOD_PLACEHOLDER, build_peptide_from_pdb_tleap_script,
        build_peptide_tleap_script, count_pdb_residues,
        peptide_pdb_preprocess, peptide_residue_count, rename_lib_unit,
        validate_user_pdb,
    )
except ImportError:
    try:
        from core import (  # type: ignore
            PEPTIDE_FRCMOD_PLACEHOLDER, build_peptide_from_pdb_tleap_script,
            build_peptide_tleap_script, count_pdb_residues,
            peptide_pdb_preprocess, peptide_residue_count, rename_lib_unit,
            validate_user_pdb,
        )
    except ImportError:  # server-side introspection (no sys.path tweak yet)
        PEPTIDE_FRCMOD_PLACEHOLDER = ""
        build_peptide_from_pdb_tleap_script = None
        build_peptide_tleap_script = None
        count_pdb_residues = None
        peptide_pdb_preprocess = None
        peptide_residue_count = None
        rename_lib_unit = None
        validate_user_pdb = None


def _get_from_predecessors(predecessor_data, key):
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


class PeptideBuilder(Node):
    """Build (or accept) a standalone peptide and emit its parameterized lib."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="peptide"),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for peptide.pdb / peptide.lib / peptide.frcmod",
        ),
        "peptide_mode": SelectParameter(
            "Peptide Mode",
            default="sequence",
            options=["sequence", "pdb"],
            docstring="'sequence' builds via tleap; 'pdb' loads a user-supplied PDB.",
        ),
        "peptide_sequence": StringParameter(
            "Peptide Sequence",
            default="HIS TYR TYR LEU ALA GLU ALA",
            docstring=(
                "Space-separated 3-letter codes (e.g., HIS TYR TYR LEU ALA GLU ALA). "
                "ACE/NME caps added automatically per n_term/c_term. "
                "Used only when peptide_mode=sequence."
            ),
        ),
        "peptide_pdb": FileParameterEdit(
            "Peptide PDB",
            default="",
            docstring=(
                "User-supplied peptide PDB (peptide_mode=pdb). Must use standard "
                "ff19SB/ff14SB residue + atom names. Leave empty to auto-discover "
                "``output_peptide_pdb`` from a snp_builder predecessor."
            ),
            optional=True,
        ),
        "n_term": SelectParameter(
            "N-terminal cap",
            default="ACE",
            options=["ACE", "charged", "neutral"],
            docstring="Sequence-mode only.",
        ),
        "c_term": SelectParameter(
            "C-terminal cap",
            default="NME",
            options=["NME", "charged", "neutral"],
            docstring="Sequence-mode only.",
        ),
        "forcefield": SelectParameter(
            "Peptide Force Field",
            default="ff19SB",
            options=["ff19SB", "ff14SB"],
            docstring="Forwarded to ep_fragment_fuse_topology so both nodes use the same protein FF.",
        ),
        "save_topology": BooleanParameter(
            "Save peptide.prmtop / peptide.rst7",
            default=False,
            docstring=(
                "Also write a peptide-only AMBER topology (useful for pre-relaxing "
                "the peptide before fuse)."
            ),
        ),
        # PDB-mode preprocessing (v1.8.0) — let users feed AlphaFold /
        # ProteinMPNN / experimental PDBs directly without external tooling.
        "chain_filter": StringParameter(
            "Chain filter (PDB mode)",
            default="",
            docstring=(
                "Comma-separated list of chains to keep (e.g., 'A' or 'A,B'). "
                "Empty / 'all' keeps every chain. PDB mode only."
            ),
        ),
        "residue_range": StringParameter(
            "Residue range (PDB mode)",
            default="",
            docstring=(
                "Residue-sequence range to keep (e.g., '5-30'). Empty / 'all' "
                "keeps everything. PDB mode only — useful for trimming an "
                "AlphaFold output to just a binding loop."
            ),
        ),
        "drop_heteroatoms": BooleanParameter(
            "Drop heteroatoms (PDB mode)",
            default=True,
            docstring=(
                "Drop HETATM records and ATOM records with non-standard "
                "residue names (waters, ligands, post-translationally "
                "modified residues that ff19SB doesn't have templates for). "
                "PDB mode only."
            ),
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting peptide build...", node_id=self.node_id, progress=0)

        if build_peptide_tleap_script is None:
            raise NodeException("setup",
                "peptide_builder core.py could not be imported — run this "
                "node in a pixi env with ambertools available.")

        try:
            result = NodeResult()
            case_name = flow_vars["case_name"].get_value() or "peptide"
            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)

            mode = flow_vars["peptide_mode"].get_value() or "sequence"
            forcefield = flow_vars["forcefield"].get_value() or "ff19SB"
            n_term = flow_vars["n_term"].get_value() or "ACE"
            c_term = flow_vars["c_term"].get_value() or "NME"
            sequence = flow_vars["peptide_sequence"].get_value() or ""
            user_pdb = self.resolve_path(flow_vars["peptide_pdb"].get_value()) or ""
            save_topology = bool(flow_vars["save_topology"].get_value())

            # In PDB mode: 3-tier resolution — explicit > predecessor > error
            if mode == "pdb" and not user_pdb:
                ref = _get_from_predecessors(predecessor_data, "output_peptide_pdb")
                if ref:
                    user_pdb = self.resolve_path(ref) or ""

            peptide_pdb = os.path.join(out_dir, "peptide.pdb")
            peptide_lib = os.path.join(out_dir, "peptide.lib")
            peptide_frcmod = os.path.join(out_dir, "peptide.frcmod")
            tleap_path = os.path.join(out_dir, "peptide.tleap")

            if mode == "pdb":
                stream_log(f"Validating user-supplied PDB: {user_pdb}",
                           node_id=self.node_id, progress=10)
                try:
                    validate_user_pdb(user_pdb)
                except ValueError as ve:
                    raise NodeException("setup", str(ve))

                chain_filter = flow_vars["chain_filter"].get_value() or ""
                residue_range = flow_vars["residue_range"].get_value() or ""
                drop_het = bool(flow_vars["drop_heteroatoms"].get_value())

                # Stage + preprocess the user PDB into work_dir. The
                # preprocessor handles MODEL 1 selection, chain + residue
                # filtering, heteroatom drop, YASARA atom rename, and HIS
                # tautomer inference — so AlphaFold / ProteinMPNN /
                # YASARA / experimental PDBs all work without external
                # cleanup.
                staged_input = os.path.join(out_dir, "peptide_input.pdb")
                stream_log(
                    f"Preprocessing PDB (chain={chain_filter or 'all'}, "
                    f"range={residue_range or 'all'}, drop_het={drop_het})...",
                    node_id=self.node_id, progress=20,
                )
                try:
                    stats = peptide_pdb_preprocess(
                        input_pdb=user_pdb,
                        output_pdb=staged_input,
                        chain_filter=chain_filter,
                        residue_range=residue_range,
                        drop_heteroatoms=drop_het,
                    )
                except ValueError as ve:
                    raise NodeException("setup", str(ve))
                pep_size = stats["residues"]
                if pep_size <= 0:
                    raise NodeException(
                        "setup",
                        "After preprocessing, no peptide residues remain. "
                        f"Stats: {stats}",
                    )
                tleap_script = build_peptide_from_pdb_tleap_script(
                    forcefield=forcefield,
                    peptide_pdb_basename="peptide_input.pdb",
                    save_topology=save_topology,
                )
            else:
                try:
                    tleap_script = build_peptide_tleap_script(
                        forcefield=forcefield,
                        peptide_sequence=sequence,
                        n_term=n_term,
                        c_term=c_term,
                        save_topology=save_topology,
                    )
                except ValueError as ve:
                    raise NodeException("setup", str(ve))
                pep_size = peptide_residue_count(sequence, n_term, c_term)
                if pep_size <= 0:
                    raise NodeException("setup", "Peptide residue count is zero")

            with open(tleap_path, "w") as f:
                f.write(tleap_script)

            # tleap's saveoff APPENDS to existing libs rather than overwriting.
            # Clean up stale outputs from prior runs in the same work_dir so we
            # don't end up with two units in peptide.lib.
            for stale in (peptide_lib, peptide_pdb, peptide_frcmod,
                          os.path.join(out_dir, "peptide.prmtop"),
                          os.path.join(out_dir, "peptide.rst7"),
                          os.path.join(out_dir, "leap.log")):
                if os.path.isfile(stale):
                    os.unlink(stale)

            stream_log(
                f"Running tleap to build {pep_size}-residue peptide ({mode} mode)...",
                node_id=self.node_id, progress=40,
            )
            proc = subprocess.run(
                ["tleap", "-f", "peptide.tleap"],
                cwd=out_dir, capture_output=True, text=True,
            )
            leap_log = os.path.join(out_dir, "leap.log")
            if os.path.isfile(leap_log):
                shutil.copy2(leap_log, os.path.join(out_dir, "peptide.leap.log"))

            if not (os.path.isfile(peptide_pdb) and os.path.getsize(peptide_pdb) > 0):
                msg = proc.stderr or proc.stdout or "tleap failed"
                raise NodeException(
                    "execution",
                    f"tleap did not produce peptide.pdb. stderr tail:\n{msg[-2000:]}",
                )
            if not (os.path.isfile(peptide_lib) and os.path.getsize(peptide_lib) > 0):
                msg = proc.stderr or proc.stdout or "tleap failed"
                raise NodeException(
                    "execution",
                    f"tleap did not produce peptide.lib. stderr tail:\n{msg[-2000:]}",
                )

            # Re-count residues from the PDB tleap actually wrote (PDB mode may
            # have had the count from validate_user_pdb; sequence mode trusts
            # peptide_residue_count). Either way, PDB-derived count is canonical.
            actual_residues = count_pdb_residues(peptide_pdb)
            if actual_residues > 0:
                pep_size = actual_residues

            # Rename the unit inside peptide.lib from tleap's default "mol"
            # to "pep" so it doesn't collide with fragment.lib's "mol" when
            # fuse does `loadoff peptide.lib + loadoff fragment.lib`.
            try:
                rename_lib_unit(peptide_lib, "pep")
            except (ValueError, OSError) as ex:
                stream_log(f"Warning: failed to rename peptide.lib unit: {ex}",
                           node_id=self.node_id, progress=85)

            # Ship a placeholder frcmod so the fuse boundary is symmetric
            with open(peptide_frcmod, "w") as f:
                f.write(PEPTIDE_FRCMOD_PLACEHOLDER)

            stream_log(
                f"Peptide ready: {pep_size} residues",
                node_id=self.node_id, progress=100,
            )

            output = {
                "case_name": case_name,
                "working_path": self.format_output_path(out_dir),
                "output_pdb": self.format_output_path(peptide_pdb),
                # Role-specific keys — these are what ep_fragment_fuse_topology and
                # ep_apply_coords resolve for the peptide side (the generic
                # output_pdb / output_lib / output_frcmod below would collide
                # with fragment_align's same-named generic keys if both nodes
                # were upstream of the same consumer).
                "output_peptide_pdb": self.format_output_path(peptide_pdb),
                "output_peptide_lib": self.format_output_path(peptide_lib),
                "output_peptide_frcmod": self.format_output_path(peptide_frcmod),
                # Generic keys kept for backward compat with consumers that
                # don't yet know about the role-specific names.
                "output_lib": self.format_output_path(peptide_lib),
                "output_frcmod": self.format_output_path(peptide_frcmod),
                "peptide_residues": pep_size,
                "forcefield": forcefield,
            }
            if save_topology:
                prmtop = os.path.join(out_dir, "peptide.prmtop")
                rst7 = os.path.join(out_dir, "peptide.rst7")
                if os.path.isfile(prmtop):
                    output["output_prmtop"] = self.format_output_path(prmtop)
                if os.path.isfile(rst7):
                    output["output_rst7"] = self.format_output_path(rst7)

            result.data = output
            result.files["output"] = {
                "pdb": output["output_pdb"],
                "lib": output["output_lib"],
                "frcmod": output["output_frcmod"],
            }
            if "output_prmtop" in output:
                result.files["output"]["prmtop"] = output["output_prmtop"]
                result.files["output"]["rst7"] = output["output_rst7"]
            result.success = True
            result.message = (
                f"Built peptide ({pep_size} residues, {forcefield}, mode={mode}); "
                f"emitted .pdb + .lib + .frcmod"
            )
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("peptide build", str(e))
