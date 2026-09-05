"""Pure-Python core for smiles-to-3d.

The science lives here; node.py is a thin Salpa wrapper around it. This file has
no Salpa imports, so it is unit-testable on its own, and RDKit is imported inside
the functions that need it so the pure helpers stay importable without it.

What it does: for each SMILES, add hydrogens, embed one or more 3D conformers
with RDKit's ETKDG (a fixed random seed, so the same input gives the same
coordinates), optimise each with the MMFF94 force field, and write every
conformer to one SDF file with its MMFF energy as a property.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

#: A SMILES line is "SMILES [name]". Lines starting with # are comments.
Record = Tuple[str, str]  # (smiles, name)


def read_smiles(text: str = "", file_path: str = "") -> List[Record]:
    """Collect (smiles, name) records from a text block and/or a file.

    One molecule per line, an optional name after whitespace (the .smi convention).
    A molecule with no name is named by its 1-based position. Blank lines and
    `#` comments are ignored. File reads pass encoding="utf-8": bare open() defaults
    to cp1252 on Windows and silently corrupts non-ASCII names.
    """
    lines: List[str] = []
    if text:
        lines += text.splitlines()
    if file_path:
        with open(file_path, encoding="utf-8") as fh:
            lines += fh.read().splitlines()
    records: List[Record] = []
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        smiles = parts[0]
        name = parts[1].strip() if len(parts) > 1 else f"mol{len(records) + 1}"
        records.append((smiles, name))
    return records


def embed_conformers(smiles: str, num_conformers: int = 1, seed: int = 42,
                     optimize: bool = True):
    """Return an RDKit Mol with `num_conformers` embedded, optimised 3D conformers
    and a list of their MMFF94 energies (kcal/mol), lowest first is NOT assumed —
    the order is the embedding order, and conformer ids match that order.

    Raises ValueError for a SMILES RDKit cannot parse, or one it cannot embed
    (which happens for a few exotic cases); the caller decides whether that is
    fatal. It is: a node must not silently skip a molecule it was asked for.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=max(1, num_conformers), params=params))
    if not conf_ids:
        raise ValueError(f"RDKit could not embed a 3D conformer for {smiles!r}")
    energies: List[float] = []
    props = AllChem.MMFFGetMoleculeProperties(mol)
    for cid in conf_ids:
        if optimize and props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
            ff.Minimize(maxIts=2000)
            energies.append(float(ff.CalcEnergy()))
        elif props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
            energies.append(float(ff.CalcEnergy()))
        else:
            # MMFF has no parameters for this molecule (rare: some metals, radicals).
            energies.append(float("nan"))
    return mol, energies


def process(records: List[Record], output_file: str, num_conformers: int = 1,
            seed: int = 42, optimize: bool = True) -> Dict:
    """Embed every record and write all conformers to `output_file` (SDF).

    Returns a plain dict — the node wrapper turns it into a NodeResult:

        molecules: one entry per input, in order — name, smiles, n_atoms (with H),
                   n_heavy_atoms, n_conformers, lowest_energy_kcal_mol
        n_in / n_out: records read / molecules written
        output_file: the SDF path as given

    Each SDF record carries the properties SMILES, conformer_id and
    MMFF94_energy_kcal_mol; its title is the molecule's name.
    """
    from rdkit import Chem

    summaries = []
    writer = Chem.SDWriter(output_file)
    try:
        for smiles, name in records:
            mol, energies = embed_conformers(smiles, num_conformers, seed, optimize)
            mol.SetProp("_Name", name)
            mol.SetProp("SMILES", smiles)
            for cid, energy in zip(range(mol.GetNumConformers()), energies):
                mol.SetProp("conformer_id", str(cid))
                mol.SetProp("MMFF94_energy_kcal_mol", f"{energy:.4f}")
                writer.write(mol, confId=cid)
            finite = [e for e in energies if e == e]  # drop NaN
            summaries.append({
                "name": name,
                "smiles": smiles,
                "n_atoms": mol.GetNumAtoms(),
                "n_heavy_atoms": mol.GetNumHeavyAtoms(),
                "n_conformers": mol.GetNumConformers(),
                "lowest_energy_kcal_mol": round(min(finite), 4) if finite else None,
            })
    finally:
        writer.close()
    return {
        "molecules": summaries,
        "n_in": len(records),
        "n_out": len(summaries),
        "output_file": output_file,
    }


def conformer_coordinates(smiles: str, seed: int = 42) -> List[Tuple[float, float, float]]:
    """The optimised coordinates of the single conformer for `smiles` — what the
    determinism test compares across two runs with the same seed."""
    mol, _ = embed_conformers(smiles, 1, seed, optimize=True)
    conf = mol.GetConformer(0)
    return [tuple(round(v, 6) for v in conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
