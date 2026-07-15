# GROMACS MD Run (Local)

A simple BoCoFlow node for running GROMACS molecular dynamics simulations locally. This node is designed for demos, tutorials, and testing purposes.

## Features

- Run complete MD simulations (grompp + mdrun) in one step
- Easy configuration via BoCoFlow UI
- Supports all standard GROMACS simulation types:
  - Energy minimization
  - NVT equilibration
  - NPT equilibration
  - Production runs

## Demo Data

This package includes ready-to-use demo data for testing:

```
demo_data/water_box/
├── water_box.gro    # 27 TIP3P water molecules in 1.5nm cube
├── topol.top        # Topology (OPLS-AA force field)
├── em.mdp           # Energy minimization parameters
├── nvt.mdp          # NVT equilibration (10ps at 300K)
└── README.md        # Usage instructions
```

### Quick Test

1. Add the node to your workflow
2. Set input files to the demo_data/water_box files
3. Run energy minimization (~5 seconds)

See [demo_data/water_box/README.md](demo_data/water_box/README.md) for details.

## Requirements

This package uses [pixi](https://pixi.sh) to manage dependencies:

- **GROMACS** >= 2023 (installed automatically via conda-forge)
- **Python** >= 3.9

## Installation

```bash
cd gmx_mdrun_local
pixi install
```

## Usage

### In BoCoFlow GUI

1. Add the "GROMACS MD Run (Local)" node to your workflow
2. Configure the input files:
   - **Topology File (.top)**: Your GROMACS topology file
   - **Structure File (.gro)**: Input structure/coordinates
   - **Parameters File (.mdp)**: MD parameters
   - **Index File (.ndx)**: Optional index file
3. Set execution options:
   - **Run Label**: Output file prefix (e.g., "nvt", "npt", "md")
   - **Number of Threads**: OpenMP threads (0 = auto-detect)
   - **Max Warnings**: Maximum grompp warnings allowed
4. Execute the node

### Standalone Usage (Level 1 API)

The core functions can be used independently without BoCoFlow:

```python
from core import run_md_simulation, check_gromacs_available

# Check GROMACS installation
if check_gromacs_available():
    result = run_md_simulation(
        top_file="topol.top",
        gro_file="conf.gro",
        mdp_file="em.mdp",
        working_dir="/path/to/simulation",
        run_label="em"
    )

    if result.success:
        print(f"Simulation completed: {result.message}")
        print(f"Output structure: {result.gro_file}")
    else:
        print(f"Simulation failed: {result.message}")
```

## Output Files

The node generates standard GROMACS output files:

| File | Description |
|------|-------------|
| `{run_label}.tpr` | Portable run input file |
| `{run_label}.gro` | Output structure |
| `{run_label}.xtc` | Trajectory (if trajectory output enabled) |
| `{run_label}.edr` | Energy file |
| `{run_label}.log` | Log file |

## Architecture

This package follows the **Node Wrapper Mechanism** pattern:

```
gmx-mdrun-local/
├── core.py       # Level 1: Pure Python functions (BoCoFlow-independent)
├── node.py       # Level 2: BoCoFlow wrapper
├── meta.toml     # Package metadata
├── pixi.toml     # Environment specification
├── demo_data/    # Example data for testing
│   └── water_box/
└── tests/        # Test suite
    ├── test_core.py   # Level 1 tests (no BoCoFlow required)
    └── test_node.py   # Level 2 tests (requires BoCoFlow)
```

### Level 1: Core Functions (`core.py`)

Pure Python functions that can be:
- Tested with standard pytest
- Used from command line
- Called from Jupyter notebooks
- Imported into other projects

### Level 2: Node Wrapper (`node.py`)

BoCoFlow integration that handles:
- Parameter extraction from UI
- Path resolution (abs:/rel: prefixes)
- Result formatting for workflow

## Testing

```bash
# Run all tests
pixi run test

# Or directly with pytest
pytest tests/ -v
```

## For HPC/SLURM Support

For running simulations on HPC clusters with SLURM, use the full `gmx-mdrun` node instead, which supports both local and remote execution modes.

## Citation

If you use this node in your research, please cite GROMACS:

> Abraham, M. J. et al. (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. SoftwareX, 1-2, 19-25. https://doi.org/10.1016/j.softx.2015.06.001

## License

MIT License
