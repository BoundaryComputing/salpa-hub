"""SnP Builder Node — Zn->Sn swap, add 2 axial OMe, cap aniline NH.

Reads a YASARA-style PDB that contains a ``UNK`` residue (Zn-porphyrin + aniline
linker) and emits a capped Sn(IV)(OMe)2-porphyrin fragment as XYZ + PDB, ready
for the existing ``easyparm-vwf`` 5-node chain.

Outputs explicit ``output_xyz``/``output_pdb``/``output_cap_atoms`` keys so
downstream nodes (``ep_bond_detection`` for QM prep, ``ep_fragment_fuse_topology`` for
the final tleap merge) can consume them via the standard 3-tier predecessor
data-flow pattern.
"""

import os

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FloatParameter, FolderParameter,
    SelectParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import build_snp_fragment, extract_peptide, write_outputs
except ImportError:  # direct-path import (node_runner adds node_dir to sys.path)
    try:
        from core import build_snp_fragment, extract_peptide, write_outputs
    except ImportError:  # server-side (no heavy deps); OPTIONS still introspectable
        build_snp_fragment = None
        extract_peptide = None
        write_outputs = None


def _parse_residue_range(raw):
    """Parse a 'lo-hi' residue-range string into an inclusive ``(int, int)``
    tuple. Empty / whitespace-only input returns ``None`` (whole chain).

    Raises ``ValueError`` on malformed input so the node can surface a
    clear setup error.
    """
    s = (raw or "").strip()
    if not s:
        return None
    parts = s.split("-")
    if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
        raise ValueError(f"expected 'lo-hi' (e.g. '1-7'), got {raw!r}")
    try:
        lo, hi = int(parts[0]), int(parts[1])
    except ValueError:
        raise ValueError(f"residue numbers must be integers, got {raw!r}")
    if lo > hi:
        raise ValueError(f"lo ({lo}) must be <= hi ({hi})")
    return (lo, hi)


class SnpBuilder(Node):
    """Convert YASARA-style Zn-porphyrin UNK residue into a capped Sn(OMe)2
    porphyrin fragment ready for EasyParm parameterization."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="snp"),
        "input_pdb": FileParameterEdit(
            "Input PDB",
            docstring="YASARA-style PDB with a UNK residue (Zn-porphyrin + aniline linker)",
        ),
        "unk_resname": StringParameter(
            "UNK Residue Name",
            default="UNK",
            docstring="Residue name of the placeholder porphyrin fragment",
        ),
        "metal_in": StringParameter(
            "Placeholder Metal",
            default="ZN",
            docstring="Atom name of the metal to replace",
        ),
        "metal_out": StringParameter(
            "Target Metal",
            default="SN",
            docstring="Atom name of the new metal center",
        ),
        "axial_ligand": SelectParameter(
            "Axial Ligand",
            default="OMe",
            options=["OMe", "OH", "Cl", "none"],
            docstring="Ligand placed above and below the porphyrin plane",
        ),
        "axial_bond_len": FloatParameter(
            "Axial Bond Length (Å)",
            default=2.00,
            docstring="Metal-ligand distance for axial ligands",
        ),
        "cap_style": SelectParameter(
            "Aniline Cap",
            default="ACE",
            options=["ACE", "H", "NHMe"],
            docstring="Cap on the aniline NH2 side; ACE gives correct amide-N electronics for QM",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for snp_frag.xyz and snp_frag.pdb",
        ),
        "extract_peptide": BooleanParameter(
            "Extract Peptide",
            default=True,
            docstring=(
                "Also write the peptide chain (everything except the UNK residue) "
                "as ``peptide_from_pdb.pdb``. Forwarded as ``output_peptide_pdb`` so "
                "peptide_builder can pick it up — gives the fused complex consistent "
                "coordinates for free, no minimization needed."
            ),
        ),
        "peptide_residue_range": StringParameter(
            "Peptide Residue Range",
            default="",
            docstring=(
                "Optional. When 'Extract Peptide' is on, restrict the extracted "
                "chain to this inclusive residue-number span, written 'lo-hi' "
                "(e.g. '1-7' to carve the Case 1 heptapeptide out of the "
                "23-residue snpp.pdb). The truncated chain is re-terminated with "
                "a synthetic TER so tleap loadpdb caps it cleanly. "
                "Empty = whole chain."
            ),
        ),
        "cap_peptide_termini": BooleanParameter(
            "Cap Peptide Termini (ACE/NME)",
            default=False,
            docstring=(
                "When 'Extract Peptide' is on, geometrically add an ACE N-cap "
                "and NME C-cap to the extracted chain (Ac-…-NH-CH3) and drop "
                "the charged-terminus atoms. Use for short peptides whose "
                "folding is the scientific question — charged termini suppress "
                "helix propensity. Default off (verbatim termini)."
            ),
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting SnP fragment build...", node_id=self.node_id, progress=0)
        if build_snp_fragment is None:
            raise NodeException("setup",
                "core.py could not be imported — run this node in a pixi env "
                "with numpy + biopython available.")

        try:
            result = NodeResult()
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}

            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "snp")
            input_pdb = self.resolve_path(flow_vars["input_pdb"].get_value())
            if not input_pdb or not os.path.isfile(input_pdb):
                raise NodeException("setup",
                    f"Input PDB not found: {input_pdb}. Set the Input PDB parameter.")

            unk = flow_vars["unk_resname"].get_value() or "UNK"
            metal_in = flow_vars["metal_in"].get_value() or "ZN"
            metal_out = flow_vars["metal_out"].get_value() or "SN"
            axial = flow_vars["axial_ligand"].get_value() or "OMe"
            bond_len = float(flow_vars["axial_bond_len"].get_value() or 2.00)
            cap = flow_vars["cap_style"].get_value() or "ACE"
            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)

            stream_log(f"Parsing {os.path.basename(input_pdb)}...",
                       node_id=self.node_id, progress=20)
            build = build_snp_fragment(
                input_pdb,
                unk_resname=unk,
                metal_in=metal_in,
                metal_out=metal_out,
                axial_ligand=axial,
                axial_bond_len=bond_len,
                cap_style=cap,
            )

            stream_log(f"Writing outputs ({len(build.atoms)} atoms)...",
                       node_id=self.node_id, progress=70)
            paths = write_outputs(build, out_dir, basename="snp_frag", resname="SNP")

            do_extract = bool(flow_vars["extract_peptide"].get_value())
            try:
                residue_range = _parse_residue_range(
                    flow_vars["peptide_residue_range"].get_value())
            except ValueError as ex:
                raise NodeException("setup",
                    f"Invalid Peptide Residue Range: {ex}")
            cap_termini = bool(flow_vars["cap_peptide_termini"].get_value())
            peptide_pdb_path = None
            n_pep_residues = 0
            if do_extract:
                peptide_pdb_path = os.path.join(out_dir, "peptide_from_pdb.pdb")
                try:
                    n_pep_residues = extract_peptide(
                        input_pdb, peptide_pdb_path, unk_resname=unk,
                        residue_range=residue_range,
                        cap_termini=cap_termini,
                    )
                except Exception as ex:
                    raise NodeException("execution",
                        f"extract_peptide failed: {ex}")
                range_note = (f" (residues {residue_range[0]}-{residue_range[1]})"
                              if residue_range else "")
                cap_note = " +ACE/NME caps" if cap_termini else ""
                stream_log(
                    f"Extracted peptide: {n_pep_residues} residues{range_note}"
                    f"{cap_note} → {peptide_pdb_path}",
                    node_id=self.node_id, progress=85,
                )

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(out_dir),
                "output_xyz": self.format_output_path(paths["xyz"]),
                "output_pdb": self.format_output_path(paths["pdb"]),
                "output_cap_atoms": paths["cap_atoms"],
                "n_atoms": len(build.atoms),
                "cap_style": cap,
                "axial_ligand": axial,
                "metal": metal_out,
            }
            if peptide_pdb_path and os.path.isfile(peptide_pdb_path):
                result.data["output_peptide_pdb"] = self.format_output_path(peptide_pdb_path)
                result.data["peptide_residues"] = n_pep_residues
            result.files["input"] = {"pdb": self.format_output_path(input_pdb)}
            result.files["output"] = {
                "xyz": self.format_output_path(paths["xyz"]),
                "pdb": self.format_output_path(paths["pdb"]),
            }
            if peptide_pdb_path and os.path.isfile(peptide_pdb_path):
                result.files["output"]["peptide_pdb"] = self.format_output_path(peptide_pdb_path)
            result.success = True
            extra = f", + peptide ({n_pep_residues} res)" if do_extract else ""
            result.message = (f"Built SnP fragment: {len(build.atoms)} atoms, "
                              f"{cap} cap, {axial} axial{extra}")

            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("snp builder", str(e))
