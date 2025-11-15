# Folding Analysis

A Python toolkit for analyzing protein folding trajectories and contact formation.

## Features

- **Trajectory Parsing**: Parse PDB trajectory files and calculate native contact formation
- **Contact Clustering**: Identify and cluster protein contacts using adjacency-based algorithms
- **Visualization**: Create contact maps, probability matrices, folding pathway plots, and trajectory analysis plots
- **Analysis**: Calculate contact probabilities, cluster formation, and folding metrics
- **Interactive CLI**: Rich-formatted command-line interface with progress bars and detailed output

## Installation

This project uses `uv` for dependency management. Activate the virtual environment:

```bash
source .venv/bin/activate
```

## Usage

### Trajectory Analysis

The main CLI tool for analyzing trajectory files. **Note**: Before analyzing trajectories, you need to parse and cluster contacts using `contacts parse` to create a CSV file with columns `i`, `j`, `r`, and `cluster`.

```bash
python run.py trajectory analyze <INPUT_FILE_PATH> [OPTIONS]
```

#### Basic Usage

```bash
# Analyze a trajectory file with default settings
python run.py trajectory analyze data/unfolded_1/133.9/105/traj.pdb
```

#### With Custom Contacts CSV File

```bash
# Specify a custom native contacts CSV file (must contain columns: i, j, r, cluster)
python run.py trajectory analyze traj.pdb --contacts data/BDD_clustered.csv
```

#### Debug Mode (Limit Frames)

```bash
# Process only first 100 frames for faster testing
python run.py trajectory analyze traj.pdb --max-frames 100
```

#### Custom Parameters

```bash
# Use custom cutoff distance and output path
python run.py trajectory analyze traj.pdb \
    --cutoff-distance 1.5 \
    --output-csv results/my_results.csv
```

#### Without Saving CSV

```bash
# Analyze without saving results to CSV
python run.py trajectory analyze traj.pdb --no-save-csv
```

#### Complete Example

```bash
# Full example with all options
python run.py trajectory analyze \
    data/unfolded_1/133.9/105/traj.pdb \
    --contacts data/BDD_clustered.csv \
    --cutoff-distance 1.2 \
    --max-frames 1000 \
    --output-csv results/trajectory_analysis.csv
```

### Trajectory Visualization

Create two-panel plots showing native contact formation and cluster filling over time:

```bash
python run.py trajectory draw <INPUT_FILE_PATH> [OPTIONS]
```

#### Basic Usage

```bash
# Create visualization by parsing PDB file
python run.py trajectory draw data/unfolded_1/133.9/105/traj.pdb
```

#### Using Pre-parsed CSV Data

```bash
# Use existing CSV file (faster if trajectory was already analyzed)
python run.py trajectory draw data/unfolded_1/133.9/105/traj.pdb \
    --csv data/unfolded_1/133.9/105/traj_contacts.csv
```

#### Custom Window Size for Running Average

```bash
# Use larger window for smoother running average
python run.py trajectory draw traj.pdb --window-size 200
```

#### Different Output Format

```bash
# Save as SVG instead of PNG
python run.py trajectory draw traj.pdb --output-format svg
```

#### Without Displaying Plot

```bash
# Save plot without displaying (useful for scripts)
python run.py trajectory draw traj.pdb --no-show
```

#### Complete Example

```bash
# Full example with all options
python run.py trajectory draw \
    data/unfolded_1/133.9/105/traj.pdb \
    --contacts data/BDD_clustered.csv \
    --csv data/unfolded_1/133.9/105/traj_contacts.csv \
    --window-size 100 \
    --output-format png \
    --no-show
```

### Contact Map Parsing

Parse contact files and perform clustering to create a CSV file with cluster assignments:

```bash
python run.py contacts parse <INPUT_FILE> [OPTIONS]
```

#### Basic Usage

```bash
# Parse contact file with default clustering (10 clusters)
python run.py contacts parse data/BDD.cont
```

#### Custom Number of Clusters

```bash
# Create 5 clusters
python run.py contacts parse data/BDD.cont --clusters 5
```

#### Custom Cutoff Size

```bash
# Use larger minimum cluster size
python run.py contacts parse data/BDD.cont --clusters 10 --cutoff-size 8
```

#### Custom Output Path

```bash
# Specify custom output CSV path
python run.py contacts parse data/BDD.cont \
    --clusters 10 \
    --output-csv results/contacts_clustered.csv
```

#### Complete Example

```bash
# Full example with all options
python run.py contacts parse data/BDD.cont \
    --clusters 10 \
    --cutoff-size 5 \
    --output-csv results/BDD_clustered.csv
```

### Contact Map Visualization

Create a contact map visualization showing contacts as colored rectangles with cluster assignments:

```bash
python run.py contacts draw <INPUT_FILE> [OPTIONS]
```

#### Basic Usage

```bash
# Create visualization by clustering and drawing
python run.py contacts draw data/BDD.cont --clusters 10
```

#### Using Pre-clustered CSV Data

```bash
# Use existing clustered CSV file (faster if contacts were already parsed)
python run.py contacts draw data/BDD.cont \
    --csv data/BDD_clustered.csv
```

#### Custom Output Format

```bash
# Save as SVG instead of PNG
python run.py contacts draw data/BDD.cont --output-format svg
```

#### Without Displaying Plot

```bash
# Save plot without displaying (useful for scripts)
python run.py contacts draw data/BDD.cont --no-show
```

#### Complete Example

```bash
# Full example with all options
python run.py contacts draw data/BDD.cont \
    --clusters 10 \
    --cutoff-size 5 \
    --output-format png \
    --no-show
```

### Command Options

#### `trajectory analyze` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process |
| `--save-csv` / `--no-save-csv` | | `True` | Save results to CSV file |
| `--output-csv` | | Auto-generated | Custom path for output CSV file |

#### `trajectory draw` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--csv` | | `None` | Path to CSV file with parsed trajectory data (if not provided, will parse PDB) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff (only used if parsing PDB) |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process (only used if parsing PDB) |
| `--window-size` | | `100` | Window size for running average |
| `--output-format` | | `png` | Output format: `png` or `svg` |
| `--no-show` | | `False` | Do not display the plot (only save it) |

#### `contacts parse` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--clusters` | | `10` | Number of clusters to create |
| `--cutoff-size` | | `5` | Minimum cluster size to retain |
| `--output-csv` | | Auto-generated | Custom path for output CSV file |

#### `contacts draw` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--csv` | | `None` | Path to CSV file with clustered contacts (if not provided, will cluster from input file) |
| `--clusters` | | `10` | Number of clusters to create (only used if not using --csv) |
| `--cutoff-size` | | `5` | Minimum cluster size to retain (only used if not using --csv) |
| `--output-format` | | `png` | Output format: `png` or `svg` |
| `--no-show` | | `False` | Do not display the plot (only save it) |

### Help

Get help for any command:

```bash
# Main help
python run.py --help

# Trajectory group help
python run.py trajectory --help

# Analyze command help
python run.py trajectory analyze --help

# Draw command help
python run.py trajectory draw --help

# Contacts group help
python run.py contacts --help

# Parse command help
python run.py contacts parse --help

# Draw command help
python run.py contacts draw --help
```

## Project Structure

```
folding_analysis/
├── data/                    # Input data files
│   ├── BDD.cont            # Native contacts file (raw)
│   ├── BDD_clustered.csv   # Clustered contacts CSV (created by contacts parse)
│   └── unfolded_*/          # Trajectory directories
├── results/                 # Output directory
│   └── plots/              # Generated plots
├── src/                    # Source code
│   ├── trajectory.py      # Trajectory parsing and analysis
│   ├── contacts.py        # Contact map loading and clustering
│   ├── clustering.py       # Contact clustering algorithms
│   ├── drawing.py         # Visualization functions
│   └── trajectory_parsing.py  # Trajectory data parsing
└── run.py                  # CLI entry point
```

## Output

### `trajectory analyze` Output

The `trajectory analyze` command generates:

1. **CSV File**: Results saved as CSV with columns:
   - `frame`: Frame number
   - `contacts`: Number of existing contacts
   - `q`: Fraction of native contacts existing
   - `contact_list`: List of existing contact pairs
   - `clusters_filling`: Dictionary of cluster formation fractions

2. **Console Output**: 
   - Progress bars for file reading and processing
   - Summary statistics (Q range, average contacts, etc.)
   - Preview of results

### `trajectory draw` Output

The `trajectory draw` command generates:

1. **Plot File**: Two-panel visualization saved to `results/plots/`:
   - **Left Panel**: Fraction of native contacts (q) vs frame number
     - Shows raw data (light blue) and running average (bold blue)
   - **Right Panel**: Cluster filling fraction vs frame number
     - One curve per cluster showing raw data and running average
   - Both panels include grid, legends, and axis labels

2. **Console Output**:
   - Configuration summary
   - Progress bars (if parsing PDB file)
   - Success message with output file path

### `contacts parse` Output

The `contacts parse` command generates:

1. **CSV File**: Contact map with clustering saved as CSV with columns:
   - `i`: First residue index
   - `j`: Second residue index
   - `r6`: r6 parameter
   - `r12`: r12 parameter
   - `r`: Calculated native distance
   - `cluster`: Cluster number (0 = unassigned, 1+ = assigned clusters)

2. **Console Output**:
   - Configuration summary
   - Contact loading confirmation
   - Clustering progress
   - Cluster distribution statistics table
   - Preview of results

### `contacts draw` Output

The `contacts draw` command generates:

1. **Plot File**: Contact map visualization saved to `results/plots/`:
   - Contacts shown as colored rectangles below the diagonal
   - Each cluster has a distinct color (using tab20 colormap)
   - Contacts with cluster=0 (unassigned) are shown in black
   - Includes grid, legend, diagonal line, and statistics text box

2. **Console Output**:
   - Configuration summary
   - Contact loading/clustering progress
   - Cluster information summary
   - Success message with output file path

## Examples

### Example 1: Quick Analysis

```bash
python run.py trajectory analyze data/unfolded_1/133.9/105/test_frame.pdb
```

### Example 2: Debug Mode with Limited Frames

```bash
python run.py trajectory analyze data/unfolded_1/133.9/105/traj.pdb \
    --max-frames 50 \
    --cutoff-distance 1.3
```

### Example 3: Custom Output Location

```bash
python run.py trajectory analyze data/unfolded_1/133.9/105/traj.pdb \
    --output-csv results/analysis_133.9_105.csv \
    --contacts data/BDD_clustered.csv
```

### Example 4: Create Visualization from PDB File

```bash
# Parse and visualize in one step
python run.py trajectory draw data/unfolded_1/133.9/105/traj.pdb \
    --window-size 100
```

### Example 5: Create Visualization from CSV (Faster)

```bash
# Use pre-parsed CSV data for faster visualization
python run.py trajectory draw data/unfolded_1/133.9/105/traj.pdb \
    --csv data/unfolded_1/133.9/105/traj_contacts.csv \
    --window-size 50 \
    --output-format svg
```

### Example 6: Workflow - Analyze then Visualize

```bash
# Step 1: Analyze trajectory and save to CSV
python run.py trajectory analyze data/unfolded_1/133.9/105/traj.pdb \
    --output-csv results/traj_contacts.csv

# Step 2: Create visualization from CSV (much faster)
python run.py trajectory draw data/unfolded_1/133.9/105/traj.pdb \
    --csv results/traj_contacts.csv \
    --window-size 100 \
    --no-show
```

### Example 7: Parse and Cluster Contacts

```bash
# Parse contact file and create 10 clusters
python run.py contacts parse data/BDD.cont --clusters 10
```

### Example 8: Custom Contact Clustering

```bash
# Create 5 clusters with custom cutoff size
python run.py contacts parse data/BDD.cont \
    --clusters 5 \
    --cutoff-size 8 \
    --output-csv results/contacts_5clusters.csv
```

### Example 9: Visualize Contact Map

```bash
# Create and visualize contact map
python run.py contacts draw data/BDD.cont --clusters 10
```

### Example 10: Visualize from Clustered CSV

```bash
# Use pre-clustered CSV for faster visualization
python run.py contacts draw data/BDD.cont \
    --csv data/BDD_clustered.csv \
    --output-format svg \
    --no-show
```

### Example 11: Workflow - Parse then Visualize

```bash
# Step 1: Parse and cluster contacts
python run.py contacts parse data/BDD.cont \
    --clusters 10 \
    --output-csv results/BDD_clustered.csv

# Step 2: Create visualization from CSV (faster)
python run.py contacts draw data/BDD.cont \
    --csv results/BDD_clustered.csv \
    --output-format png \
    --no-show
```

### Example 12: Complete Workflow - Contacts to Trajectory Analysis

```bash
# Step 1: Parse and cluster contacts (required before trajectory analysis)
python run.py contacts parse data/BDD.cont \
    --clusters 10 \
    --output-csv data/BDD_clustered.csv

# Step 2: Analyze trajectory using the clustered contacts CSV
python run.py trajectory analyze data/unfolded_1/133.9/105/traj.pdb \
    --contacts data/BDD_clustered.csv \
    --cutoff-distance 1.2 \
    --output-csv results/traj_analysis.csv

# Step 3: Visualize trajectory results
python run.py trajectory draw data/unfolded_1/133.9/105/traj.pdb \
    --contacts data/BDD_clustered.csv \
    --csv results/traj_analysis.csv \
    --window-size 100
```

## Requirements

- Python 3.9+
- pandas
- numpy
- matplotlib
- seaborn
- rich
- click
- scikit-learn

## License

[Add your license here]

