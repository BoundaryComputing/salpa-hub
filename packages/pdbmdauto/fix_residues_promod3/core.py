"""
fix-residues-promod3 core — pure Python logic, no BoCoFlow dependencies.

Fills missing residues in PDB structures using ProMod3, the engine behind
SWISS-MODEL. Uses the ProMod3 Python API for full control over the pipeline:

    BuildRawModel → FillLoopsByDatabase → CloseGaps → ModelTermini → BuildSidechains

Energy minimization (MinimizeModelEnergy) is skipped because it crashes on
macOS ARM64 due to an OpenMM issue. The downstream GROMACS pipeline performs
its own energy minimization anyway.

ProMod3 replaces MODELLER — Apache 2.0 license, no registration required.
"""

import json
import os
import sys
from dataclasses import dataclass, field

# Beside fixed.pdb: which model residues ProMod3 rebuilt, chain by chain. Original Atom
# Groups reads it, because the model's numbers are not the deposited ones.
REBUILT_MAP_NAME = "rebuilt_residues.json"


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class FixResiduesResult:
    """Result of missing residue repair."""

    output_pdb: str
    num_chains_processed: int = 0
    total_residues_added: int = 0
    chain_details: dict = field(
        default_factory=dict
    )  # chain_id -> {before, after, added}
    promod3_log: str = ""
    success: bool = False
    rebuilt_map: str = ""


# ---------------------------------------------------------------------------
# Alignment format conversion
# ---------------------------------------------------------------------------


def ali_to_fasta(
    ali_path: str,
    output_path: str,
    chain_id: str = "",
    pdb_stem: str = "template",
) -> str:
    """Convert a MODELLER .ali file to ProMod3-compatible FASTA alignment.

    ProMod3 convention: template header is `<PDB_FILENAME_STEM>.<CHAIN_ID>`.
    Target must be named "target" or "trg" or be the first sequence.
    """
    with open(ali_path, "r") as f:
        content = f.read()

    blocks = [b.strip() for b in content.split("\n\n") if b.strip()]
    if len(blocks) != 2:
        raise ValueError(f"Expected 2 blocks in .ali file, got {len(blocks)}")

    def extract_sequence(block: str) -> str:
        lines = block.strip().split("\n")
        seq_lines = lines[2:]
        seq = "".join(seq_lines).rstrip("*")
        return seq

    template_seq = extract_sequence(blocks[0])
    target_seq = extract_sequence(blocks[1])

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    chain_suffix = f".{chain_id}" if chain_id else ""
    with open(output_path, "w") as f:
        f.write(f">target\n{target_seq}\n")
        f.write(f">{pdb_stem}{chain_suffix}\n{template_seq}\n")

    return output_path


# ---------------------------------------------------------------------------
# The record of rebuilt residues
# ---------------------------------------------------------------------------


def rebuilt_positions(fasta_path: str) -> list:
    """Target positions (1-based) where the template has no residue.

    ProMod3 numbers the model by position in the target sequence, so these are the
    numbers the rebuilt residues carry in fixed.pdb. A column where the target has a
    gap is not a target position and shifts nothing.
    """
    sequences, current = [], None
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current = []
                sequences.append(current)
            elif current is not None:
                current.append(line)
    if len(sequences) < 2:
        raise ValueError(
            f"{os.path.basename(fasta_path)} holds fewer than two sequences"
        )
    target, template = "".join(sequences[0]), "".join(sequences[1])
    if len(target) != len(template):
        raise ValueError(
            f"{os.path.basename(fasta_path)}: the two rows differ in length"
        )
    position, rebuilt = 0, []
    for t, m in zip(target, template):
        if t == "-":
            continue
        position += 1
        if m == "-":
            rebuilt.append(position)
    return rebuilt


def model_residue_numbers(pdb_path: str) -> list:
    """[(chain, [residue numbers in order])], chains in the order the file lists them."""
    order, numbers = [], {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain, resid = line[21], int(line[22:26])
            if chain not in numbers:
                order.append(chain)
                numbers[chain] = []
            if not numbers[chain] or numbers[chain][-1] != resid:
                numbers[chain].append(resid)
    return [(chain, numbers[chain]) for chain in order]


def write_rebuilt_map(fixed_pdb: str, fasta_by_chain: dict, out_path: str) -> str:
    """Write which residues of fixed.pdb ProMod3 rebuilt, chain by chain in file order.

    A rebuilt residue is one at a target position the template lacked AND present in
    the model: with terminal modelling off, a missing tag is absent, not rebuilt.
    """
    chains = []
    for chain, residues in model_residue_numbers(fixed_pdb):
        fasta = fasta_by_chain.get(chain)
        gaps = set(rebuilt_positions(fasta)) if fasta else set()
        chains.append(
            {
                "chain": chain,
                "residues": residues,
                "rebuilt": [n for n in residues if n in gaps],
            }
        )
    record = {
        "numbering": "model",
        "note": (
            "Residue numbers as in fixed.pdb, where ProMod3 numbers every chain from 1 "
            "by position in the target sequence. 'rebuilt' lists the residues ProMod3 "
            "built because the template had no coordinates for them."
        ),
        "chains": chains,
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2)
        f.write("\n")
    return out_path


# ---------------------------------------------------------------------------
# ProMod3 Python API pipeline
# ---------------------------------------------------------------------------


def fix_chain_residues(
    pdb_path: str,
    fasta_path: str,
    chain_id: str,
) -> dict:
    """Fix missing residues in a single protein chain using ProMod3 Python API.

    Pipeline: BuildRawModel → FillLoopsByDatabase → CloseGaps → ModelTermini → BuildSidechains

    Args:
        pdb_path: Path to template PDB.
        fasta_path: Path to ProMod3-format FASTA alignment.
        chain_id: Chain ID to process.

    Returns:
        dict with keys: model (ost Entity), residues_before, residues_after, gaps_filled
    """
    from ost import io, mol, seq
    from promod3 import loop, modelling

    # Load template, restrict to this chain
    tpl = io.LoadPDB(pdb_path, restrict_chains=chain_id)
    mol.alg.AssignSecStruct(tpl)

    # Load alignment
    aln = io.LoadAlignment(fasta_path)
    aln.AttachView(1, tpl.CreateFullView())

    # Build raw model
    mhandle = modelling.BuildRawModel(aln)
    residues_before = mhandle.model.GetResidueCount()
    n_gaps = len(mhandle.gaps)

    # Load fragment databases
    frag_db = loop.LoadFragDB()
    struct_db = loop.LoadStructureDB()
    torsion_sampler = loop.LoadTorsionSamplerCoil()

    # Fill internal loops from fragment database
    modelling.FillLoopsByDatabase(mhandle, frag_db, struct_db, torsion_sampler)

    # Close remaining internal gaps (Monte Carlo fallback)
    modelling.CloseGaps(
        mhandle,
        fragment_db=frag_db,
        structure_db=struct_db,
        torsion_sampler=torsion_sampler,
    )

    # Model terminal extensions (N/C-terminal missing residues)
    modelling.ModelTermini(mhandle, torsion_sampler)

    # Rebuild sidechains
    modelling.BuildSidechains(mhandle)

    # Skip MinimizeModelEnergy — crashes on macOS ARM64 (OpenMM issue).
    # Downstream GROMACS performs its own energy minimization.

    residues_after = mhandle.model.GetResidueCount()

    return {
        "model": mhandle.model,
        "residues_before": residues_before,
        "residues_after": residues_after,
        "gaps_initial": n_gaps,
        "gaps_remaining": len(mhandle.gaps),
    }


# ---------------------------------------------------------------------------
# High-level processing
# ---------------------------------------------------------------------------


def process_fix_residues(
    pdb_path: str,
    ali_dir: str,
    output_dir: str,
    case_name: str,
    chain_ids: list,
    protein_chains: list = None,
    model_termini: bool = True,
    merge_folder_name: str = "Merge",
) -> FixResiduesResult:
    """Fix missing residues in a PDB structure using ProMod3.

    Processes each protein chain independently (ProMod3 requires single-chain
    alignment). Combines results into a single output PDB.

    Args:
        pdb_path: Path to merged PDB (from merge_pdb_chains).
        ali_dir: Directory with per-chain .ali files (from gen_ali).
        output_dir: Base output directory.
        case_name: Case identifier.
        chain_ids: All chain IDs in the structure.
        protein_chains: Protein-only chain IDs. If None, uses all chain_ids.
        model_termini: Whether to model terminal extensions (default True).
        merge_folder_name: Subfolder for output.

    Returns:
        FixResiduesResult with repaired PDB path and statistics.
    """
    chains_to_model = protein_chains if protein_chains else chain_ids

    merge_dir = os.path.join(output_dir, merge_folder_name)
    os.makedirs(merge_dir, exist_ok=True)

    # PDB filename stem for alignment headers
    pdb_stem = os.path.splitext(os.path.basename(pdb_path))[0]

    # Process each protein chain
    chain_models = []
    chain_details = {}
    alignment_files = {}
    failed_chains = []
    total_added = 0
    log_lines = []

    for chain_id in sorted(chains_to_model):
        ali_path = os.path.join(ali_dir, chain_id, "homology.ali")
        if not os.path.exists(ali_path):
            log_lines.append(f"Chain {chain_id}: FAILED — no alignment at {ali_path}")
            chain_details[chain_id] = {"error": "no alignment"}
            failed_chains.append(chain_id)
            continue

        # Convert .ali to ProMod3 FASTA
        fasta_path = os.path.join(merge_dir, f"alignment_{chain_id}.fasta")
        ali_to_fasta(ali_path, fasta_path, chain_id=chain_id, pdb_stem=pdb_stem)
        alignment_files[chain_id] = fasta_path

        # Fix residues for this chain
        try:
            result = fix_chain_residues(pdb_path, fasta_path, chain_id)
            added = result["residues_after"] - result["residues_before"]
            total_added += added
            chain_models.append((chain_id, result["model"]))
            chain_details[chain_id] = {
                "residues_before": result["residues_before"],
                "residues_after": result["residues_after"],
                "added": added,
                "gaps_initial": result["gaps_initial"],
                "gaps_remaining": result["gaps_remaining"],
            }
            log_lines.append(
                f"Chain {chain_id}: {result['residues_before']} → {result['residues_after']} "
                f"residues (+{added}), {result['gaps_initial']} gaps → {result['gaps_remaining']}"
            )
        except Exception as e:
            log_lines.append(f"Chain {chain_id}: FAILED — {e}")
            chain_details[chain_id] = {"error": str(e)}
            failed_chains.append(chain_id)

    # Combine chain models into a single PDB
    output_pdb = os.path.join(merge_dir, "fixed.pdb")

    # A protein chain left out here would be missing from every later step, which
    # then runs on a different system from the one deposited. Up to 1.2.5 the step
    # still reported success; now nothing is written and the step fails.
    if failed_chains:
        log_lines.append(
            f"\nfixed.pdb not written: chain(s) {', '.join(failed_chains)} could not be "
            f"modelled, and a structure without them is not the deposited system."
        )
        return FixResiduesResult(
            output_pdb=output_pdb,
            num_chains_processed=len(chain_models),
            total_residues_added=total_added,
            chain_details=chain_details,
            promod3_log="\n".join(log_lines),
            success=False,
        )

    if chain_models:
        from ost import io, mol

        # Save each chain model separately, then combine via simple file concat
        # (OST entity API for combining chains is complex; PDB concat is simpler)
        chain_pdbs = []
        for chain_id, model in chain_models:
            chain_pdb = os.path.join(merge_dir, f"_chain_{chain_id}_fixed.pdb")
            # Rename chain from 'A' (ProMod3 default) to actual chain ID
            for chain in model.GetChainList():
                edi = model.EditXCS()
                edi.RenameChain(chain, chain_id)
                break
            io.SavePDB(model, chain_pdb)
            chain_pdbs.append(chain_pdb)

        # Combine by reading all ATOM/HETATM lines
        with open(output_pdb, "w") as out:
            for cpdb in chain_pdbs:
                with open(cpdb) as f:
                    for line in f:
                        if line.startswith(("ATOM", "HETATM", "TER")):
                            out.write(line)
            out.write("END\n")

        # Clean up temp files
        for cpdb in chain_pdbs:
            os.remove(cpdb)

        log_lines.append(f"\nCombined output: {output_pdb}")
        rebuilt_map = write_rebuilt_map(
            output_pdb, alignment_files, os.path.join(merge_dir, REBUILT_MAP_NAME)
        )
        log_lines.append(f"Rebuilt residues recorded in {rebuilt_map}")
        success = True
    else:
        log_lines.append("\nNo chains were successfully processed")
        rebuilt_map = ""
        success = False

    return FixResiduesResult(
        output_pdb=output_pdb,
        num_chains_processed=len(chain_models),
        total_residues_added=total_added,
        chain_details=chain_details,
        promod3_log="\n".join(log_lines),
        success=success,
        rebuilt_map=rebuilt_map,
    )
