# Water Box Demo Data

A minimal water box simulation case for testing the GROMACS MD Run (Local) node.

## Contents

| File | Description |
|------|-------------|
| `water_box.gro` | Coordinates for 27 TIP3P water molecules in a 1.5nm cube |
| `topol.top` | Topology file using OPLS-AA force field |
| `em.mdp` | Energy minimization parameters (steepest descent) |
| `nvt.mdp` | NVT equilibration parameters (10 ps at 300K) |

## Quick Start

### Using the BoCoFlow Node

1. Add the "GROMACS MD Run (Local)" node to your workflow
2. Configure the node with these demo files:
   - **Input GRO file**: `water_box.gro`
   - **Input TOP file**: `topol.top`
   - **Input MDP file**: `em.mdp` (for energy minimization)
3. Execute the node

### Using Command Line

```bash
# Energy minimization
gmx grompp -f em.mdp -c water_box.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v

# NVT equilibration (after minimization)
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v
```

## Expected Results

### Energy Minimization
- Should complete in < 1000 steps
- Final potential energy: ~ -1000 to -2000 kJ/mol
- Maximum force: < 1000 kJ/mol/nm

### NVT Equilibration
- Temperature should stabilize around 300K
- Total energy should be relatively constant
- Runtime: ~10 seconds on modern hardware

## System Details

- **Water model**: TIP3P (OPLS-AA force field)
- **Box size**: 1.5 nm × 1.5 nm × 1.5 nm
- **Number of molecules**: 27 water molecules (81 atoms)
- **Periodic boundaries**: XYZ

## Notes

- This is a minimal demo system - real simulations typically use larger boxes
- The OPLS-AA force field is included with GROMACS by default
- For larger systems, use `gmx solvate` to generate water boxes

## References

- [GROMACS Documentation](https://manual.gromacs.org/)
- [GROMACS Tutorials](https://tutorials.gromacs.org/)
- [TIP3P Water Model](https://en.wikipedia.org/wiki/TIP3P_water_model)
