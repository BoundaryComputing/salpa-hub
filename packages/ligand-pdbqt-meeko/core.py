"""Pure-Python core for ligand-pdbqt-meeko.

The science lives here; node.py is a thin Salpa wrapper around it. No Salpa
imports, so this file is unit-testable on its own; RDKit and Meeko are imported
inside the functions that need them so the pure helpers stay importable.

What it does: read every record of an SDF, make sure it carries explicit
hydrogens with coordinates, hand it to Meeko's MoleculePreparation (Gasteiger
charges, AutoDock atom types, the rotatable-bond tree) and write one PDBQT per
record, named by the record's title. Meeko merges non-polar hydrogens into their
heavy atoms, which is why a PDBQT has fewer atoms than the SDF record.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional


def safe_name(title: str, index: int) -> str:
    """A file name from an SDF title: letters, digits, `-`, `_`; positional otherwise."""
    cleaned = re.sub(r"[^A-Za-z0-9_-]+", "_", (title or "").strip()).strip("_")
    return cleaned or f"ligand{index + 1}"


def read_sdf(sdf_path: str) -> List:
    """Every record of the SDF as an RDKit Mol, hydrogens kept. A record RDKit cannot
    read is a ValueError naming its position: a docking run on three of four requested
    ligands is a wrong answer, not a partial one."""
    from rdkit import Chem

    mols = []
    for i, mol in enumerate(Chem.SDMolSupplier(sdf_path, removeHs=False)):
        if mol is None:
            raise ValueError(f"record {i + 1} of {os.path.basename(sdf_path)} could not be read as a molecule")
        mols.append(mol)
    if not mols:
        raise ValueError(f"{os.path.basename(sdf_path)} holds no molecules")
    return mols


def ensure_hydrogens(mol) -> tuple:
    """Return (mol with explicit hydrogens and 3D coordinates, hydrogens_were_added).

    Meeko needs explicit hydrogens to place polar ones and merge the rest. An SDF from
    a 2D sketch or one stripped of hydrogens gets them added with coordinates; the
    caller reports that, because coordinates invented for hydrogens are a modelling
    decision the user should know about.
    """
    from rdkit import Chem

    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        raise ValueError(f"{mol.GetProp('_Name') if mol.HasProp('_Name') else 'a record'} has no 3D coordinates — run it through a conformer generator first")
    heavy_h = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
    implicit = sum(a.GetTotalNumHs(includeNeighbors=False) for a in mol.GetAtoms())
    if implicit == 0 and heavy_h > 0:
        return mol, False
    return Chem.AddHs(mol, addCoords=True), True


def prepare_pdbqt(mol, rigid_macrocycles: bool = False, hydrate: bool = False,
                  flexible_amides: bool = False) -> Dict:
    """One record → its PDBQT string and the numbers a docking user checks.

    Returns: pdbqt (str), n_atoms_pdbqt, n_atoms_input, torsions (TORSDOF),
    total_charge, hydrogens_added. Raises ValueError with Meeko's reason when a
    record cannot be prepared.
    """
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol, added = ensure_hydrogens(mol)
    prep = MoleculePreparation(
        rigid_macrocycles=rigid_macrocycles, hydrate=hydrate, flexible_amides=flexible_amides,
    )
    setups = prep.prepare(mol)
    if not setups:
        raise ValueError("Meeko produced no setup for this record")
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"Meeko could not write a PDBQT: {err}")
    atoms = [ln for ln in pdbqt.splitlines() if ln.startswith(("ATOM", "HETATM"))]
    tors = [ln for ln in pdbqt.splitlines() if ln.startswith("TORSDOF")]
    charge = sum(float(ln[70:76]) for ln in atoms)
    return {
        "pdbqt": pdbqt,
        "n_atoms_pdbqt": len(atoms),
        "n_atoms_input": mol.GetNumAtoms(),
        "torsions": int(tors[0].split()[1]) if tors else 0,
        "total_charge": round(charge, 3),
        "hydrogens_added": added,
    }


def process(sdf_path: str, output_dir: str, rigid_macrocycles: bool = False,
            hydrate: bool = False, flexible_amides: bool = False) -> Dict:
    """Prepare every record of `sdf_path` and write `<name>.pdbqt` files into `output_dir`.

    Returns a plain dict — the node wrapper turns it into a NodeResult:
        ligands: one entry per record, in order — name, file, n_atoms_pdbqt,
                 n_atoms_input, torsions, total_charge, hydrogens_added
        n_in / n_out: records read / files written
        pdbqt_files: the paths, in order
    """
    mols = read_sdf(sdf_path)
    os.makedirs(output_dir, exist_ok=True)
    ligands = []
    seen: Dict[str, int] = {}
    for i, mol in enumerate(mols):
        title = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        name = safe_name(title, i)
        if name in seen:                      # two records with one title still get two files
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 1
        try:
            out = prepare_pdbqt(mol, rigid_macrocycles, hydrate, flexible_amides)
        except ValueError as e:
            raise ValueError(f"record {i + 1} ({title or 'untitled'}): {e}")
        path = os.path.join(output_dir, f"{name}.pdbqt")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(out["pdbqt"])
        ligands.append({
            "name": name, "file": path,
            "n_atoms_pdbqt": out["n_atoms_pdbqt"], "n_atoms_input": out["n_atoms_input"],
            "torsions": out["torsions"], "total_charge": out["total_charge"],
            "hydrogens_added": out["hydrogens_added"],
        })
    return {
        "ligands": ligands,
        "n_in": len(mols),
        "n_out": len(ligands),
        "pdbqt_files": [lig["file"] for lig in ligands],
    }
