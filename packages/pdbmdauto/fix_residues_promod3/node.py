"""
fix-residues-promod3 — BoCoFlow node wrapper.

Fills missing residues in PDB structures using ProMod3 (SWISS-MODEL engine).
Reads the merged PDB and per-chain alignment files from upstream nodes,
converts to ProMod3 format, and runs `pm build-model`.

Input: predecessor data from merge_pdb_chains (output_pdb, chain_types)
Output: Repaired PDB with missing residues filled
"""

import os
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter,
    StringParameter,
)

try:
    from .core import process_fix_residues
except ImportError:
    process_fix_residues = None


class FixResiduesPromod3(Node):
    """
    Fix missing residues in PDB structures using ProMod3.

    Uses the ProMod3 homology modeling engine (SWISS-MODEL) to fill
    gaps in protein structures. ProMod3 searches a fragment database
    (~19.7M fragments), performs loop modeling, rebuilds sidechains,
    and runs energy minimization.

    Requires:
    - ProMod3 installed (`pixi add promod3` from bioconda channel)
    - Per-chain .ali alignment files (from gen_ali)
    - Merged PDB structure (from merge_pdb_chains)

    Only protein chains are modeled — DNA/RNA chains are preserved as-is.

    Input: predecessor data (working_path, output_pdb, chain_types)
    Output: Repaired PDB file in Merge/ subdirectory
    """

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            default="",
            docstring="Case identifier. Leave empty to use predecessor data.",
        ),
        "model_termini": BooleanParameter(
            "Model Terminal Extensions",
            default=False,
            docstring="If true, also model terminal extensions (N/C-terminal residues beyond the template).",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute missing residue repair with ProMod3."""
        try:
            from bocoflow_core.stream_logger import stream_log
        except ImportError:
            stream_log = lambda msg, **kw: log_message(msg)

        try:
            result = NodeResult()
            stream_log("Starting missing residue repair (ProMod3)", node_id=self.node_id, progress=0)

            # --- Parse predecessor data (flat dict) ---
            if not predecessor_data or not predecessor_data[0]:
                raise NodeException("fix_residues_promod3", "No predecessor data.")
            input_data = predecessor_data[0]

            # --- Resolve parameters ---
            case_name = (
                flow_vars["case_name"].get_value()
                or input_data.get("case_name")
                or "protein"
            )
            model_termini = flow_vars["model_termini"].get_value()

            # Resolve output directory
            working_path = input_data.get("working_path", "")
            output_dir = self.resolve_path(working_path) if working_path else ""
            if not output_dir:
                raise NodeException("fix_residues_promod3", "No output directory from predecessor.")

            # Find merged PDB. Prefer the explicit path from the upstream
            # merge_pdb_chains node (input_data["output_pdb"]) — scanning the
            # merge_dir for *.pdb files breaks on re-runs because this node's
            # own output (fixed.pdb) also lives in merge_dir and alphabetises
            # before merge.pdb. If that happens, pdb_stem becomes "fixed"
            # instead of "merge", ali_to_fasta writes alignment headers
            # `>fixed.A/B`, and fix_chain_residues loads fixed.pdb as the
            # template — ProMod3 sees the previously-filled loops in the
            # template positions the alignment expects to be gaps, emitting
            # "Alignment-structure mismatch at pos N" for every position
            # that differs from the gap-filled loop residue. See pdbmdauto
            # force-rerun bug discovered 2026-04-23.
            pdb_path = None
            merge_folder = input_data.get("merge_folder", "")
            merge_dir = self.resolve_path(merge_folder) if merge_folder \
                else os.path.join(output_dir, "Merge")

            output_pdb_hint = input_data.get("output_pdb", "")
            if output_pdb_hint:
                candidate = self.resolve_path(output_pdb_hint)
                if os.path.exists(candidate):
                    pdb_path = candidate

            # Fallback scan — exclude this node's own output so a stale
            # fixed.pdb from a prior run can never be picked up as input.
            SELF_OUTPUT_NAMES = {"fixed.pdb"}
            if not pdb_path and os.path.isdir(merge_dir):
                for f in sorted(os.listdir(merge_dir)):
                    if f.endswith(".pdb") and f not in SELF_OUTPUT_NAMES:
                        pdb_path = os.path.join(merge_dir, f)
                        break
            if not pdb_path:
                # Last-resort scan in output_dir with the same exclusion.
                for f in sorted(os.listdir(output_dir)):
                    if f.endswith(".pdb") and f not in SELF_OUTPUT_NAMES:
                        pdb_path = os.path.join(output_dir, f)
                        break
            if not pdb_path:
                raise NodeException("fix_residues_promod3", "No PDB file found.")

            log_message(f"Case: {case_name}, PDB: {pdb_path}")
            stream_log("Converting alignments to ProMod3 format", node_id=self.node_id, progress=10)

            # Get chain info from predecessor
            chain_types = input_data.get("chain_types", {})
            chain_list = input_data.get("pdb_chain_list", input_data.get("selected_pdb_chain_list", []))
            protein_chains = [c for c, t in chain_types.items() if t == "P1"]

            if not protein_chains:
                protein_chains = chain_list  # Fallback: treat all as protein

            # Warn about DNA/RNA chains (ProMod3 is protein-only)
            dna_rna_chains = [c for c, t in chain_types.items() if t == "DL"]
            if dna_rna_chains:
                warning = (
                    f"DNA/RNA chains detected ({', '.join(dna_rna_chains)}) but not supported "
                    f"by ProMod3. Only protein chains will be modeled. "
                    f"DNA/RNA chains will be absent from the output PDB."
                )
                log_message(f"WARNING: {warning}")
                stream_log(f"Warning: {warning}", node_id=self.node_id, progress=15)

            log_message(f"Protein chains to model: {protein_chains}")
            stream_log("Running ProMod3 (per-chain gap filling)", node_id=self.node_id, progress=20)

            # --- Process ---
            fix_result = process_fix_residues(
                pdb_path=pdb_path,
                ali_dir=output_dir,
                output_dir=output_dir,
                case_name=case_name,
                chain_ids=chain_list,
                protein_chains=protein_chains,
                model_termini=model_termini,
            )

            if not fix_result.success:
                # Stream the per-chain failure log to the UI so the actual
                # ProMod3 error (caught inside core.process_fix_residues)
                # surfaces without needing to dig into the file log.
                log_message(f"ProMod3 log:\n{fix_result.promod3_log}")
                stream_log(
                    f"ProMod3 build-model failed:\n{fix_result.promod3_log}",
                    node_id=self.node_id,
                    progress=100,
                    level="error",
                )
                # Include the tail of the promod3_log in the NodeException
                # message so the /api/logs error record carries it too.
                tail = (fix_result.promod3_log or "").strip()
                if len(tail) > 1500:
                    tail = tail[-1500:]
                raise NodeException(
                    "fix_residues_promod3",
                    f"ProMod3 build-model failed:\n{tail}" if tail
                    else "ProMod3 build-model failed. Check logs for details."
                )

            stream_log("ProMod3 completed successfully", node_id=self.node_id, progress=90)

            # --- Build result ---
            formatted_output_dir = self.format_output_path(output_dir)

            result.files["output"]["fixed_pdb"] = self.format_output_path(fix_result.output_pdb)

            result.data.update({
                "case_name": case_name,
                "working_path": formatted_output_dir,
                "output_pdb": self.format_output_path(fix_result.output_pdb),
                "pdb_chain_list": chain_list,
                "chain_types": chain_types,
                "chain_info": input_data.get("chain_info", {}),
                "num_chains_modeled": fix_result.num_chains_processed,
                "total_residues_added": fix_result.total_residues_added,
                "chain_details": fix_result.chain_details,
            })

            result.metadata.update({
                "case_name": case_name,
                "output_dir": formatted_output_dir,
                "execution_time": datetime.now().isoformat(),
            })

            result.success = True
            result.message = (
                f"Fixed missing residues in {fix_result.num_chains_processed} protein chain(s): "
                f"+{fix_result.total_residues_added} residues added by ProMod3"
            )

            log_message(result.message)
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            log_message(f"Error in FixResiduesPromod3: {str(e)}")
            raise NodeException("fix_residues_promod3", str(e))
