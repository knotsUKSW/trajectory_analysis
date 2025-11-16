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

### Contact Map Parsing

The primary method for parsing contact files and performing clustering. This command loads contacts, performs clustering, saves the results to CSV, and by default creates a contact map visualization.

```bash
python run.py contacts parse <INPUT_FILE> [OPTIONS]
```

**Note**: This is the recommended first step before analyzing trajectories. The output CSV file (with columns `i`, `j`, `r`, and `cluster`) is required for trajectory analysis.

#### Basic Usage

```bash
# Parse contact file with default settings (10 clusters, creates visualization)
python run.py contacts parse data/BDD.cont
```

This will:
- Load contacts from `data/BDD.cont`
- Perform clustering (10 clusters by default)
- Save CSV to `results/contacts_11_clusters.csv` (number may vary)
- Create and save contact map plot to `results/plots/contacts_11_clusters.png`

#### Custom Number of Clusters

```bash
# Create 5 clusters instead of default 10
python run.py contacts parse data/BDD.cont --clusters 5
```

#### Custom Cutoff Size

```bash
# Use larger minimum cluster size (filters out smaller clusters)
python run.py contacts parse data/BDD.cont \
    --clusters 10 \
    --cutoff-size 8
```

#### Custom Output Name

```bash
# Specify custom base name for output files
python run.py contacts parse data/BDD.cont \
    --clusters 10 \
    --output-name my_contacts
```

This will create:
- `results/my_contacts.csv`
- `results/plots/my_contacts.png`

#### Disable Visualization

```bash
# Parse without creating the contact map plot
python run.py contacts parse data/BDD.cont --no-draw
```

#### Custom Plot Format

```bash
# Save plot as SVG instead of PNG
python run.py contacts parse data/BDD.cont \
    --output-format svg
```

#### Without Displaying Plot

```bash
# Save plot without displaying (useful for scripts)
python run.py contacts parse data/BDD.cont --no-show
```

#### Complete Example

```bash
# Full example with all options
python run.py contacts parse data/BDD.cont \
    --clusters 10 \
    --cutoff-size 5 \
    --output-name BDD_clustered \
    --output-format png \
    --no-show
```

#### Contact Map Visualization

**`draw`** - Creates a contact map visualization from an existing clustered CSV file.

```bash
python run.py contacts draw <CSV_FILE> [OPTIONS]
```

Parameters: `--output-name` (default: auto-generated from CSV filename), `--output-format` (default: `png`), `--no-show` (flag, default: `False`)

### Trajectory Analysis

The main CLI tool for analyzing trajectory files. **Note**: Before analyzing trajectories, you must first parse and cluster contacts using `contacts parse` (see Contact Map Parsing section above) to create a CSV file with columns `i`, `j`, `r`, and `cluster`.

#### Main Method: `parse`

The `parse` command orchestrates the complete trajectory analysis pipeline, performing all steps in sequence and saving all output files in the same directory as the input trajectory file.

```bash
python run.py trajectory parse <INPUT_FILE_PATH> [OPTIONS]
```

**What it does:**

The `parse` command performs the following steps in sequence:

1. **Reads trajectory**: Parses the PDB trajectory file and calculates which native contacts are formed in each frame, saving results as `{base}_parsed.csv`
2. **Creates visualization**: Generates a two-panel plot showing fraction of native contacts (q) and cluster filling over time, saving as `{base}_trajectory_analysis.{format}`
3. **Summarizes (float mode)**: Calculates mean cluster filling fractions in windows, saving as `{base}_summarized.csv` and plotting as `{base}_summarized.{format}`
4. **Summarizes (binary mode)**: Repeats summarization with binary conversion (0 or 1 based on cutoff), saving as `{base}_summarized_binary.csv` and plotting as `{base}_summarized_binary.{format}`
5. **Classifies trajectory**: Determines cluster formation order by reverse-tracking cluster breaks, saving as `{base}_class.txt`
6. **Optionally animates**: Creates an animated contact map showing contact formation over time, saving as `{base}_animation.gif` (if `--animate` flag is used)

All output files are saved in the same directory as the input trajectory file.

**Parameters:**

- `INPUT_FILE_PATH` (required): Path to the PDB trajectory file to analyze
- `--contacts`, `-c` (default: `data/BDD_clustered.csv`): Path to native contacts CSV file (must contain columns: i, j, r, cluster)
- `--cutoff-distance` (default: `1.2`): Multiplier for native distance cutoff
- `--max-frames` (default: `None`, all frames): Maximum number of frames to process (useful for debugging)
- `--window-size` (default: `10000`): Window size for summarization (number of frames per window)
- `--cutoff` (default: `0.5`): Cutoff for binary classification (values >= cutoff become 1, values < cutoff become 0)
- `--output-format` (default: `png`): Output format for plots - `png` or `svg`
- `--show-plot` (flag, default: `False`): Display plots (default is to only save them)
- `--animate` (flag, default: `False`): Create animation showing contact formation over time

**Basic Usage:**

```bash
# Run complete analysis pipeline with default settings
python run.py trajectory parse data/unfolded_1/133.9/105/traj.pdb
```

This will create all output files in the same directory as `traj.pdb`:
- `traj_parsed.csv`
- `traj_trajectory_analysis.png`
- `traj_summarized.csv`
- `traj_summarized.png`
- `traj_summarized_binary.csv`
- `traj_summarized_binary.png`
- `traj_class.txt`

**With Animation:**

```bash
# Include animation in the pipeline
python run.py trajectory parse traj.pdb --animate
```

This will also create:
- `traj_animation.gif`

**Custom Parameters:**

```bash
# Use custom cutoff distance, window size, and binary cutoff
python run.py trajectory parse traj.pdb \
    --cutoff-distance 1.5 \
    --window-size 5000 \
    --cutoff 0.6 \
    --output-format svg
```

**Debug Mode:**

```bash
# Process only first 100 frames for faster testing
python run.py trajectory parse traj.pdb --max-frames 100
```

**Complete Example:**

```bash
# Full example with all options
python run.py trajectory parse \
    data/unfolded_1/133.9/105/traj.pdb \
    --contacts data/BDD_clustered.csv \
    --cutoff-distance 1.2 \
    --max-frames 10000 \
    --window-size 10000 \
    --cutoff 0.5 \
    --output-format png \
    --animate
```

#### Other Methods

**`read`** - Reads trajectory file and calculates native contact formation for each frame.

```bash
python run.py trajectory read <INPUT_FILE_PATH> [OPTIONS]
```

Parameters: `--contacts`, `-c` (default: `data/BDD_clustered.csv`), `--cutoff-distance` (default: `1.2`), `--max-frames` (default: `None`), `--save-csv` / `--no-save-csv` (default: `True`), `--output-csv` (default: auto-generated as `{base}_parsed.csv`)

**`draw`** - Creates two-panel plots showing native contact formation (q) and cluster filling over time.

```bash
python run.py trajectory draw <INPUT_FILE_PATH> [OPTIONS]
```

Parameters: `--contacts`, `-c` (default: `data/BDD_clustered.csv`), `--csv` (default: `None`), `--cutoff-distance` (default: `1.2`, only used if reading PDB), `--max-frames` (default: `None`, only used if reading PDB), `--window-size` (default: `100`), `--output-format` (default: `png`), `--no-show` (flag, default: `False`)

**`summarize`** - Summarizes trajectory by calculating mean cluster filling fractions in windows.

```bash
python run.py trajectory summarize <INPUT_FILE_PATH> [OPTIONS]
```

Parameters: `--contacts`, `-c` (default: `data/BDD_clustered.csv`), `--csv` (default: `None`), `--cutoff-distance` (default: `1.2`, only used if reading PDB), `--max-frames` (default: `None`, only used if reading PDB), `--window-size` (default: `10000`), `--output-csv` (default: auto-generated), `--cutoff` (default: `None`), `--plot` / `--no-plot` (default: `False`), `--output-format` (default: `png`, only used with `--plot`), `--no-show` (flag, default: `False`, only used with `--plot`)

**`plot-summary`** - Plots cluster formation time series from trajectory summary.

```bash
python run.py trajectory plot-summary <INPUT_FILE_PATH> [OPTIONS]
```

Parameters: `--contacts`, `-c` (default: `data/BDD_clustered.csv`), `--summary-csv` (default: `None`), `--csv` (default: `None`, only used if `summary-csv` is not provided), `--cutoff-distance` (default: `1.2`, only used if generating summary), `--max-frames` (default: `None`, only used if generating summary), `--window-size` (default: `10000`, only used if generating summary), `--cutoff` (default: `None`, only used if generating summary), `--output-format` (default: `png`), `--no-show` (flag, default: `False`)

**`classify`** - Classifies trajectory by determining cluster formation order through reverse-tracking cluster breaks.

```bash
python run.py trajectory classify <INPUT_FILE_PATH> [OPTIONS]
```

Parameters: `--contacts`, `-c` (default: `data/BDD_clustered.csv`), `--summary-csv` (default: `None`), `--csv` (default: `None`, only used if `summary-csv` is not provided), `--cutoff-distance` (default: `1.2`, only used if generating summary), `--max-frames` (default: `None`, only used if generating summary), `--window-size` (default: `10000`, only used if generating summary), `--cutoff` (default: `0.5`, only used if generating summary), `--output` (default: auto-generated as `{base}_class.txt`)

**`animate`** - Creates an animated contact map showing which contacts are formed in each frame.

```bash
python run.py trajectory animate <INPUT_FILE_PATH> [OPTIONS]
```

Parameters: `--contacts`, `-c` (default: `data/BDD_clustered.csv`), `--csv` (default: `None`), `--cutoff-distance` (default: `1.2`, only used if reading PDB), `--max-frames` (default: `None`, only used if reading PDB), `--interval` (default: `50`), `--output-format` (default: `gif`), `--fps` (default: `20`, only used for `mp4`)


### Command Options

#### `trajectory parse` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process (useful for debugging) |
| `--window-size` | | `10000` | Window size for summarization (number of frames per window) |
| `--cutoff` | | `0.5` | Cutoff for binary classification (values >= cutoff become 1, values < cutoff become 0) |
| `--output-format` | | `png` | Output format for plots: `png` or `svg` |
| `--show-plot` | | `False` | Display plots (default is to only save them) |
| `--animate` | | `False` | Create animation showing contact formation over time |

#### `trajectory read` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process |
| `--save-csv` / `--no-save-csv` | | `True` | Save results to CSV file |
| `--output-csv` | | Auto-generated | Custom path for output CSV file (default: `{base}_parsed.csv`) |

#### `trajectory draw` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--csv` | | `None` | Path to CSV file with parsed trajectory data (if not provided, will read PDB) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff (only used if reading PDB) |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process (only used if reading PDB) |
| `--window-size` | | `100` | Window size for running average |
| `--output-format` | | `png` | Output format: `png` or `svg` |
| `--no-show` | | `False` | Do not display the plot (only save it) |

#### `trajectory summarize` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--csv` | | `None` | Path to CSV file with parsed trajectory data (if not provided, will read PDB) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff (only used if reading PDB) |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process (only used if reading PDB) |
| `--window-size` | | `10000` | Number of frames per window |
| `--output-csv` | | Auto-generated | Custom path for output summary CSV file |
| `--cutoff` | | `None` | Convert probabilities to binary (0 or 1). Values >= cutoff become 1, values < cutoff become 0. If not specified, returns float probabilities |
| `--plot` / `--no-plot` | | `False` | Create a plot showing each cluster as a time series |
| `--output-format` | | `png` | Output format for the plot: `png` or `svg` (only used with --plot) |
| `--no-show` | | `False` | Do not display the plot (only save it, only used with --plot) |

#### `trajectory plot-summary` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--summary-csv` | | `None` | Path to summary CSV file (if not provided, will generate summary from trajectory data) |
| `--csv` | | `None` | Path to CSV file with parsed trajectory data (only used if summary-csv is not provided) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff (only used if generating summary) |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process (only used if generating summary) |
| `--window-size` | | `10000` | Number of frames per window (only used if generating summary) |
| `--cutoff` | | `None` | Convert probabilities to binary (0 or 1) before plotting. Values >= cutoff become 1, values < cutoff become 0. If not specified, plots float probabilities (only used if generating summary) |
| `--output-format` | | `png` | Output format: `png` or `svg` |
| `--no-show` | | `False` | Do not display the plot (only save it) |

#### `trajectory classify` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--summary-csv` | | `None` | Path to binary summary CSV file (if not provided, will generate from trajectory) |
| `--csv` | | `None` | Path to CSV file with parsed trajectory data (only used if summary-csv is not provided) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff (only used if generating summary) |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process (only used if generating summary) |
| `--window-size` | | `10000` | Number of frames per window (only used if generating summary) |
| `--cutoff` | | `0.5` | Cutoff for binary classification (only used if generating summary) |
| `--output` | | Auto-generated | Path for output classification file (default: `{base}_class.txt`) |

#### `trajectory animate` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--contacts` | `-c` | `data/BDD_clustered.csv` | Path to native contacts CSV file (must contain columns: i, j, r, cluster) |
| `--csv` | | `None` | Path to CSV file with parsed trajectory data (if not provided, will read PDB) |
| `--cutoff-distance` | | `1.2` | Multiplier for native distance cutoff (only used if reading PDB) |
| `--max-frames` | | `None` (all frames) | Maximum number of frames to process (only used if reading PDB) |
| `--interval` | | `50` | Delay between frames in milliseconds |
| `--output-format` | | `gif` | Output format: `gif` or `mp4` |
| `--fps` | | `20` | Frames per second for saved animation (only used for mp4) |

#### `contacts parse` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--clusters` | | `10` | Number of clusters to create |
| `--cutoff-size` | | `5` | Minimum cluster size to retain |
| `--output-name` | | Auto-generated | Base name for output files (default: `contacts_{n}_clusters` where n is number of clusters) |
| `--no-draw` | | `False` | Do not create contact map visualization (default: visualization is created) |
| `--output-format` | | `png` | Output format for the plot: `png` or `svg` (only used if drawing) |
| `--no-show` | | `False` | Do not display the plot (only save it, only used if drawing) |

#### `contacts draw` Options

| Option | Short | Default | Description |
|--------|-------|---------|------------|
| `--output-name` | | Auto-generated | Base name for output plot file (default: from CSV file name) |
| `--output-format` | | `png` | Output format: `png` or `svg` |
| `--no-show` | | `False` | Do not display the plot (only save it) |

## Project Structure

```
folding_analysis/
├── src/                    # Source code
│   ├── trajectory.py      # Trajectory parsing and analysis
│   ├── contacts.py        # Contact map loading and clustering
│   └── utils.py           # Utility functions
└── run.py                  # CLI entry point
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

