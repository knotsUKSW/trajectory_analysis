# Folding Analysis Rust Library

This is a Rust implementation of the trajectory analysis functionality, designed to be faster than the Python version for processing large trajectory files.

## Current Features

- **PDB File Reading**: Extract CA atom coordinates from PDB trajectory files
- **Trajectory Trait**: Trait-based design for extensibility

## Usage

```rust
use folding_analysis_rs::{PdbTrajectory, Trajectory};

let trajectory = PdbTrajectory::new("path/to/trajectory.pdb");
let frames_data = trajectory.read_pdb(None)?; // Read all frames
// or
let frames_data = trajectory.read_pdb(Some(100))?; // Read first 100 frames

// Access frame data
for (frame_num, frame_data) in frames_data {
    println!("Frame {} has {} residues", frame_num, frame_data.len());
    for (residue_num, coord) in frame_data {
        println!("  Residue {}: ({}, {}, {})", residue_num, coord.x, coord.y, coord.z);
    }
}
```

## Future: Python Bindings

This library is designed to work with Python bindings using PyO3 and maturin. To add Python bindings:

1. Add PyO3 to `Cargo.toml`:
```toml
[dependencies]
pyo3 = { version = "0.21", features = ["extension-module"] }
```

2. Use `maturin` to build Python wheels:
```bash
maturin develop  # For development
maturin build   # For distribution
```

3. The library can then be imported in Python:
```python
import folding_analysis_rs

trajectory = folding_analysis_rs.PdbTrajectory("path/to/trajectory.pdb")
frames_data = trajectory.read_pdb(max_frames=None)
```

## Building

```bash
cargo build --release
```

## Testing

```bash
cargo test
```

