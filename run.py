#!/usr/bin/env python3
"""
CLI script for trajectory analysis using Click.
"""
import click
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from trajectory import Trajectory
from contacts import ContactMap
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
import pandas as pd
import ast


@click.group()
def cli():
    """Trajectory analysis CLI tool."""
    pass


@cli.group()
def trajectory():
    """Trajectory analysis commands."""
    pass


@trajectory.command('parse')
@click.argument('input_file_path', type=click.Path(exists=True, readable=True))
@click.option('--contacts', '-c', 
              type=click.Path(exists=True, readable=True),
              default='data/BDD_clustered.csv',
              help='Path to native contacts CSV file (must contain columns: i, j, r, cluster). Default: data/BDD_clustered.csv')
@click.option('--cutoff-distance', 
              type=float,
              default=1.2,
              help='Multiplier for native distance cutoff (default: 1.2)')
@click.option('--max-frames',
              type=int,
              default=None,
              help='Maximum number of frames to process (for debugging, default: all frames)')
@click.option('--window-size',
              type=int,
              default=10000,
              help='Window size for summarization (default: 10000)')
@click.option('--cutoff',
              type=float,
              default=0.5,
              help='Cutoff for binary classification (default: 0.5)')
@click.option('--output-format',
              type=click.Choice(['png', 'svg'], case_sensitive=False),
              default='png',
              help='Output format for plots (default: png)')
@click.option('--animate',
              is_flag=True,
              default=False,
              help='Create animation (default: False)')
def parse_trajectory(input_file_path, contacts, cutoff_distance, max_frames, window_size, cutoff, output_format, animate):
    """
    Main method for parsing trajectory: reads, draws, summarizes, and classifies.
    
    This command orchestrates the complete trajectory analysis pipeline:
    1. Reads trajectory and saves as _parsed.csv
    2. Draws trajectory visualization
    3. Summarizes trajectory (float mode) and saves as _summarized.csv
    4. Plots summary and saves as _summarized.{format}
    5. Summarizes trajectory (binary mode) and saves as _summarized_binary.csv
    6. Plots binary summary and saves as _summarized_binary.{format}
    7. Classifies trajectory and saves as _class.txt
    8. Optionally animates trajectory and saves as .gif
    
    All output files are saved in the same directory as the trajectory file.
    
    INPUT_FILE_PATH: Path to the PDB trajectory file to analyze
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]TRAJECTORY PARSING (COMPLETE PIPELINE)[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    config_table.add_row("📁 Input PDB file:", f"[cyan]{input_file_path}[/cyan]")
    config_table.add_row("📁 Contacts file:", f"[cyan]{contacts}[/cyan]")
    config_table.add_row("⚙️  Cutoff distance:", f"[yellow]{cutoff_distance}[/yellow]")
    if max_frames:
        config_table.add_row("🔢 Max frames:", f"[yellow]{max_frames}[/yellow] [dim](debug mode)[/dim]")
    else:
        config_table.add_row("🔢 Max frames:", "[green]All[/green]")
    config_table.add_row("📊 Window size:", f"[yellow]{window_size}[/yellow]")
    config_table.add_row("🔢 Binary cutoff:", f"[yellow]{cutoff}[/yellow]")
    config_table.add_row("🖼️  Output format:", f"[cyan]{output_format}[/cyan]")
    config_table.add_row("🎬 Animation:", f"[green]{'Yes' if animate else 'No'}[/green]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Validate contacts CSV file
        if not os.path.exists(contacts):
            raise FileNotFoundError(f"Contacts CSV file not found: {contacts}")
        
        # Check if it's a CSV file
        if not contacts.lower().endswith('.csv'):
            console.print(f"[yellow]Warning:[/yellow] File '{contacts}' does not have .csv extension. "
                        f"Trajectory class expects a CSV file with columns: i, j, r, cluster")
        
        # Create Trajectory object
        console.print("📊 [bold]Loading native contacts from CSV...[/bold]")
        trajectory = Trajectory(input_file_path, contacts)
        console.print(f"✅ [green]Loaded {trajectory.total_contacts} native contacts[/green]")
        
        # Run complete parsing pipeline
        console.print(f"\n🚀 [bold]Starting complete trajectory analysis pipeline...[/bold]")
        trajectory.parse(
            cutoff_distance=cutoff_distance,
            max_frames=max_frames,
            window_size=window_size,
            cutoff=cutoff,
            output_format=output_format,
            animate=animate
        )
        
        # Get output directory
        traj_dir = os.path.dirname(input_file_path) or '.'
        base_name = os.path.splitext(os.path.basename(input_file_path))[0]
        
        console.print(f"\n✅ [bold green]Complete pipeline finished![/bold green]")
        console.print(f"\n📁 [bold]Output files (in {traj_dir}):[/bold]")
        console.print(f"  📄 {base_name}_parsed.csv")
        console.print(f"  🖼️  {base_name}_trajectory_analysis.{output_format}")
        console.print(f"  📄 {base_name}_summarized.csv")
        console.print(f"  🖼️  {base_name}_summarized.{output_format}")
        console.print(f"  📄 {base_name}_summarized_binary.csv")
        console.print(f"  🖼️  {base_name}_summarized_binary.{output_format}")
        console.print(f"  📝 {base_name}_class.txt")
        if animate:
            console.print(f"  🎬 {base_name}_animation.gif")
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@trajectory.command('read')
@click.argument('input_file_path', type=click.Path(exists=True, readable=True))
@click.option('--contacts', '-c', 
              type=click.Path(exists=True, readable=True),
              default='data/BDD_clustered.csv',
              help='Path to native contacts CSV file (must contain columns: i, j, r, cluster). Default: data/BDD_clustered.csv')
@click.option('--cutoff-distance', 
              type=float,
              default=1.2,
              help='Multiplier for native distance cutoff (default: 1.2)')
@click.option('--max-frames',
              type=int,
              default=None,
              help='Maximum number of frames to process (for debugging, default: all frames)')
def read_trajectory(input_file_path, contacts, cutoff_distance, max_frames):
    """
    Read trajectory file and calculate native contact formation.
    
    Reads the PDB trajectory file and calculates which native contacts are formed
    in each frame. Saves results as {base}_parsed.csv in the trajectory directory.
    
    INPUT_FILE_PATH: Path to the PDB trajectory file to read
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]TRAJECTORY READING[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    config_table.add_row("📁 Input PDB file:", f"[cyan]{input_file_path}[/cyan]")
    config_table.add_row("📁 Contacts file:", f"[cyan]{contacts}[/cyan]")
    config_table.add_row("⚙️  Cutoff distance:", f"[yellow]{cutoff_distance}[/yellow]")
    if max_frames:
        config_table.add_row("🔢 Max frames:", f"[yellow]{max_frames}[/yellow] [dim](debug mode)[/dim]")
    else:
        config_table.add_row("🔢 Max frames:", "[green]All[/green]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Validate contacts CSV file
        if not os.path.exists(contacts):
            raise FileNotFoundError(f"Contacts CSV file not found: {contacts}")
        
        # Check if it's a CSV file
        if not contacts.lower().endswith('.csv'):
            console.print(f"[yellow]Warning:[/yellow] File '{contacts}' does not have .csv extension. "
                        f"Trajectory class expects a CSV file with columns: i, j, r, cluster")
        
        # Create Trajectory object
        console.print("📊 [bold]Loading native contacts from CSV...[/bold]")
        trajectory = Trajectory(input_file_path, contacts)
        console.print(f"✅ [green]Loaded {trajectory.total_contacts} native contacts[/green]")
        
        # Read trajectory (always saves to CSV)
        console.print(f"\n🔍 [bold]Reading trajectory...[/bold]")
        trajectory.read_trajectory(
            cutoff_distance=cutoff_distance,
            max_frames=max_frames
        )
        
        # Load the saved CSV to get statistics
        base_name = os.path.splitext(os.path.basename(input_file_path))[0]
        output_dir = os.path.dirname(input_file_path) or '.'
        output_path = os.path.join(output_dir, f"{base_name}_parsed.csv")
        
        if not os.path.exists(output_path):
            raise ValueError("Failed to read trajectory - output file not found")
        
        df = pd.read_csv(output_path)
        
        if df.empty:
            raise ValueError("Failed to read trajectory or trajectory is empty")
        
        console.print(f"\n✅ [bold green]Reading complete![/bold green]")
        console.print(f"📊 Total frames processed: [cyan]{len(df)}[/cyan]")
        console.print(f"📈 Q range: [cyan]{df['q'].min():.3f} - {df['q'].max():.3f}[/cyan]")
        console.print(f"📊 Average contacts per frame: [cyan]{df['contacts'].mean():.1f}[/cyan]")
        
        console.print(f"\n💾 [bold]Output file:[/bold]")
        console.print(f"📄 CSV: [cyan]{output_path}[/cyan]")
        
        # Display first few rows
        console.print(f"\n📄 [bold]First few rows:[/bold]")
        console.print(df.head().to_string())
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@trajectory.command('draw')
@click.argument('input_file_path', type=click.Path(exists=True, readable=True))
@click.option('--contacts', '-c',
              type=click.Path(exists=True, readable=True),
              default='data/BDD_clustered.csv',
              help='Path to native contacts CSV file (must contain columns: i, j, r, cluster). Default: data/BDD_clustered.csv')
@click.option('--csv', 
              type=click.Path(exists=True, readable=True),
              default=None,
              help='Path to CSV file with parsed trajectory data (if not provided, will read PDB)')
@click.option('--cutoff-distance', 
              type=float,
              default=1.2,
              help='Multiplier for native distance cutoff (default: 1.2, only used if reading PDB)')
@click.option('--max-frames',
              type=int,
              default=None,
              help='Maximum number of frames to process (only used if reading PDB, default: all frames)')
@click.option('--window-size',
              type=int,
              default=100,
              help='Window size for running average (default: 100)')
@click.option('--output-format',
              type=click.Choice(['png', 'svg'], case_sensitive=False),
              default='png',
              help='Output format for the plot (default: png)')
@click.option('--no-show',
              is_flag=True,
              default=False,
              help='Do not display the plot (only save it)')
def draw_trajectory(input_file_path, contacts, csv, cutoff_distance, max_frames, window_size, output_format, no_show):
    """
    Create visualization plots for trajectory analysis.
    
    Creates two-panel plots showing:
    - Left: Fraction of native contacts (q) vs frame number
    - Right: Cluster filling fraction vs frame number
    
    INPUT_FILE_PATH: Path to the PDB trajectory file (or base path if using --csv)
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]TRAJECTORY VISUALIZATION[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    if csv:
        config_table.add_row("📄 CSV file:", f"[cyan]{csv}[/cyan]")
        config_table.add_row("📁 PDB file (reference):", f"[cyan]{input_file_path}[/cyan]")
    else:
        config_table.add_row("📁 Input PDB file:", f"[cyan]{input_file_path}[/cyan]")
    config_table.add_row("📁 Contacts file:", f"[cyan]{contacts}[/cyan]")
    if not csv:
        config_table.add_row("⚙️  Cutoff distance:", f"[yellow]{cutoff_distance}[/yellow]")
        if max_frames:
            config_table.add_row("🔢 Max frames:", f"[yellow]{max_frames}[/yellow] [dim](debug mode)[/dim]")
        else:
            config_table.add_row("🔢 Max frames:", "[green]All[/green]")
    config_table.add_row("📊 Window size:", f"[yellow]{window_size}[/yellow]")
    config_table.add_row("🖼️  Output format:", f"[cyan]{output_format}[/cyan]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Validate contacts CSV file
        if not os.path.exists(contacts):
            raise FileNotFoundError(f"Contacts CSV file not found: {contacts}")
        
        # Check if it's a CSV file
        if not contacts.lower().endswith('.csv'):
            console.print(f"[yellow]Warning:[/yellow] File '{contacts}' does not have .csv extension. "
                        f"Trajectory class expects a CSV file with columns: i, j, r, cluster")
        
        # Create Trajectory object (it will load the CSV internally)
        console.print("📊 [bold]Loading native contacts from CSV...[/bold]")
        trajectory = Trajectory(input_file_path, contacts)
        console.print(f"✅ [green]Loaded {trajectory.total_contacts} native contacts[/green]")
        
        # Read trajectory if needed
        if not csv:
            console.print(f"\n🔍 [bold]Reading trajectory...[/bold]")
            trajectory.read_trajectory(
                cutoff_distance=cutoff_distance,
                max_frames=max_frames
            )
            console.print(f"✅ [green]Trajectory read and saved[/green]")
        
        # Create visualization
        console.print(f"\n🎨 [bold]Creating visualization...[/bold]")
        
        trajectory.draw(
            window_size=window_size,
            output_format=output_format
        )
        
        # Determine output path
        base_name = os.path.splitext(os.path.basename(input_file_path))[0]
        output_path = f'results/plots/{base_name}_trajectory_analysis.{output_format}'
        
        console.print(f"\n✅ [bold green]Visualization complete![/bold green]")
        console.print(f"💾 Plot saved to: [cyan]{output_path}[/cyan]")
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@trajectory.command('summarize')
@click.argument('input_file_path', type=click.Path(exists=True, readable=True))
@click.option('--contacts', '-c',
              type=click.Path(exists=True, readable=True),
              default='data/BDD_clustered.csv',
              help='Path to native contacts CSV file (must contain columns: i, j, r, cluster). Default: data/BDD_clustered.csv')
@click.option('--csv',
              type=click.Path(exists=True, readable=True),
              default=None,
              help='Path to CSV file with parsed trajectory data (if not provided, will parse PDB file)')
@click.option('--cutoff-distance',
              type=float,
              default=1.2,
              help='Multiplier for native distance cutoff (default: 1.2, only used if parsing PDB)')
@click.option('--max-frames',
              type=int,
              default=None,
              help='Maximum number of frames to process (only used if parsing PDB, default: all frames)')
@click.option('--window-size',
              type=int,
              default=10000,
              help='Number of frames per window (default: 10000)')
@click.option('--output-csv',
              type=click.Path(),
              default=None,
              help='Path for output summary CSV file (default: auto-generated from input path)')
@click.option('--cutoff',
              type=float,
              default=None,
              help='Convert probabilities to binary (0 or 1). Values >= cutoff become 1, values < cutoff become 0. If not specified, returns float probabilities.')
@click.option('--plot/--no-plot',
              default=False,
              help='Create a plot showing each cluster as a time series (default: False)')
@click.option('--output-format',
              type=click.Choice(['png', 'svg'], case_sensitive=False),
              default='png',
              help='Output format for the plot (default: png, only used with --plot)')
@click.option('--no-show',
              is_flag=True,
              default=False,
              help='Do not display the plot (only save it, only used with --plot)')
def summarize_trajectory(input_file_path, contacts, csv, cutoff_distance, max_frames, window_size, output_csv, cutoff, plot, output_format, no_show):
    """
    Summarize trajectory by calculating mean cluster filling fractions in windows.
    
    Groups frames into windows and calculates the mean fraction of contacts
    established for each cluster in each window. Output is a CSV file with
    windows as rows and clusters as columns.
    
    INPUT_FILE_PATH: Path to the PDB trajectory file to analyze
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]TRAJECTORY SUMMARY[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    if csv:
        config_table.add_row("📄 CSV file:", f"[cyan]{csv}[/cyan]")
        config_table.add_row("📁 PDB file (reference):", f"[cyan]{input_file_path}[/cyan]")
    else:
        config_table.add_row("📁 Input PDB file:", f"[cyan]{input_file_path}[/cyan]")
        config_table.add_row("⚙️  Cutoff distance:", f"[yellow]{cutoff_distance}[/yellow]")
        if max_frames:
            config_table.add_row("🔢 Max frames:", f"[yellow]{max_frames}[/yellow] [dim](debug mode)[/dim]")
        else:
            config_table.add_row("🔢 Max frames:", "[green]All[/green]")
    config_table.add_row("📁 Contacts file:", f"[cyan]{contacts}[/cyan]")
    config_table.add_row("📊 Window size:", f"[yellow]{window_size}[/yellow]")
    if cutoff is not None:
        config_table.add_row("🔢 Binary cutoff:", f"[yellow]{cutoff}[/yellow] [dim](binary mode)[/dim]")
    else:
        config_table.add_row("🔢 Binary cutoff:", "[green]None[/green] [dim](float mode)[/dim]")
    if output_csv:
        config_table.add_row("📄 Output CSV:", f"[cyan]{output_csv}[/cyan]")
    config_table.add_row("🎨 Create plot:", f"[green]{'Yes' if plot else 'No'}[/green]")
    if plot:
        config_table.add_row("🖼️  Output format:", f"[cyan]{output_format}[/cyan]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Validate contacts CSV file
        if not os.path.exists(contacts):
            raise FileNotFoundError(f"Contacts CSV file not found: {contacts}")
        
        # Check if it's a CSV file
        if not contacts.lower().endswith('.csv'):
            console.print(f"[yellow]Warning:[/yellow] File '{contacts}' does not have .csv extension. "
                        f"Trajectory class expects a CSV file with columns: i, j, r, cluster")
        
        # Create Trajectory object (it will load the CSV internally)
        console.print("📊 [bold]Loading native contacts from CSV...[/bold]")
        trajectory = Trajectory(input_file_path, contacts)
        console.print(f"✅ [green]Loaded {trajectory.total_contacts} native contacts[/green]")
        
        # Read trajectory if needed
        if not csv:
            console.print(f"\n🔍 [bold]Reading trajectory...[/bold]")
            trajectory.read_trajectory(
                cutoff_distance=cutoff_distance,
                max_frames=max_frames
            )
            console.print(f"✅ [green]Trajectory read and saved[/green]")
        
        # Create summary
        console.print(f"\n📊 [bold]Creating trajectory summary...[/bold]")
        if cutoff is not None:
            console.print(f"[dim]Using binary mode with cutoff={cutoff}[/dim]")
        trajectory.summarize_trajectory(
            window_size=window_size,
            cutoff=cutoff
        )
        
        # Determine output path and load for display
        base_name = os.path.splitext(os.path.basename(input_file_path))[0]
        output_dir = os.path.dirname(input_file_path) or '.'
        if cutoff is not None:
            actual_output_csv = os.path.join(output_dir, f"{base_name}_summary_binary_{cutoff}.csv")
        else:
            actual_output_csv = os.path.join(output_dir, f"{base_name}_summary.csv")
        
        # Load summary for display
        summary_df = pd.read_csv(actual_output_csv, index_col=0)
        
        console.print(f"\n✅ [bold green]Summary complete![/bold green]")
        console.print(f"📊 Summary shape: [cyan]{summary_df.shape}[/cyan] (windows x clusters)")
        console.print(f"📈 Number of windows: [cyan]{len(summary_df)}[/cyan]")
        console.print(f"🔢 Number of clusters: [cyan]{len(summary_df.columns)}[/cyan]")
        console.print(f"💾 Summary saved to: [cyan]{actual_output_csv}[/cyan]")
        
        # Display preview
        console.print(f"\n📋 [bold]Summary preview:[/bold]")
        console.print(summary_df.head().to_string())
        
        # Create plot if requested
        if plot:
            console.print(f"\n🎨 [bold]Creating summary plot...[/bold]")
            trajectory.plot_summary(
                cutoff=cutoff,
                output_format=output_format
            )
            
            # Determine plot output path
            base_name = os.path.splitext(os.path.basename(input_file_path))[0]
            plot_path = f'results/plots/{base_name}_summary.{output_format}'
            console.print(f"✅ [green]Plot saved to:[/green] [cyan]{plot_path}[/cyan]")
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@trajectory.command('plot-summary')
@click.argument('input_file_path', type=click.Path(exists=True, readable=True))
@click.option('--contacts', '-c',
              type=click.Path(exists=True, readable=True),
              default='data/BDD_clustered.csv',
              help='Path to native contacts CSV file (must contain columns: i, j, r, cluster). Default: data/BDD_clustered.csv')
@click.option('--summary-csv',
              type=click.Path(exists=True, readable=True),
              default=None,
              help='Path to summary CSV file (if not provided, will generate summary from trajectory data)')
@click.option('--csv',
              type=click.Path(exists=True, readable=True),
              default=None,
              help='Path to CSV file with parsed trajectory data (only used if summary-csv is not provided)')
@click.option('--cutoff-distance',
              type=float,
              default=1.2,
              help='Multiplier for native distance cutoff (default: 1.2, only used if generating summary)')
@click.option('--max-frames',
              type=int,
              default=None,
              help='Maximum number of frames to process (only used if generating summary, default: all frames)')
@click.option('--window-size',
              type=int,
              default=10000,
              help='Number of frames per window (default: 10000, only used if generating summary)')
@click.option('--cutoff',
              type=float,
              default=None,
              help='Convert probabilities to binary (0 or 1) before plotting. Values >= cutoff become 1, values < cutoff become 0. If not specified, plots float probabilities. Only used if generating summary.')
@click.option('--output-format',
              type=click.Choice(['png', 'svg'], case_sensitive=False),
              default='png',
              help='Output format for the plot (default: png)')
@click.option('--no-show',
              is_flag=True,
              default=False,
              help='Do not display the plot (only save it)')
def plot_summary(input_file_path, contacts, summary_csv, csv, cutoff_distance, max_frames, window_size, cutoff, output_format, no_show):
    """
    Plot cluster formation time series from trajectory summary.
    
    Creates a plot showing each cluster column from the summary DataFrame as a time series.
    Can use an existing summary CSV file or generate a summary on the fly.
    
    INPUT_FILE_PATH: Path to the PDB trajectory file (used as reference or for generating summary)
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]TRAJECTORY SUMMARY PLOT[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    config_table.add_row("📁 Input PDB file:", f"[cyan]{input_file_path}[/cyan]")
    if summary_csv:
        config_table.add_row("📄 Summary CSV:", f"[cyan]{summary_csv}[/cyan]")
    else:
        if csv:
            config_table.add_row("📄 Trajectory CSV:", f"[cyan]{csv}[/cyan]")
        else:
            config_table.add_row("⚙️  Cutoff distance:", f"[yellow]{cutoff_distance}[/yellow]")
            if max_frames:
                config_table.add_row("🔢 Max frames:", f"[yellow]{max_frames}[/yellow] [dim](debug mode)[/dim]")
            else:
                config_table.add_row("🔢 Max frames:", "[green]All[/green]")
            config_table.add_row("📊 Window size:", f"[yellow]{window_size}[/yellow]")
    config_table.add_row("📁 Contacts file:", f"[cyan]{contacts}[/cyan]")
    if cutoff is not None:
        config_table.add_row("🔢 Binary cutoff:", f"[yellow]{cutoff}[/yellow] [dim](binary mode)[/dim]")
    else:
        config_table.add_row("🔢 Binary cutoff:", "[green]None[/green] [dim](float mode)[/dim]")
    config_table.add_row("🖼️  Output format:", f"[cyan]{output_format}[/cyan]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Validate contacts CSV file
        if not os.path.exists(contacts):
            raise FileNotFoundError(f"Contacts CSV file not found: {contacts}")
        
        # Check if it's a CSV file
        if not contacts.lower().endswith('.csv'):
            console.print(f"[yellow]Warning:[/yellow] File '{contacts}' does not have .csv extension. "
                        f"Trajectory class expects a CSV file with columns: i, j, r, cluster")
        
        # Create Trajectory object (it will load the CSV internally)
        console.print("📊 [bold]Loading native contacts from CSV...[/bold]")
        trajectory = Trajectory(input_file_path, contacts)
        console.print(f"✅ [green]Loaded {trajectory.total_contacts} native contacts[/green]")
        
        # Load or generate summary
        summary_df = None
        if summary_csv:
            # Load existing summary CSV
            console.print(f"\n📄 [bold]Loading summary from CSV...[/bold]")
            if not os.path.exists(summary_csv):
                raise FileNotFoundError(f"Summary CSV file not found: {summary_csv}")
            summary_df = pd.read_csv(summary_csv, index_col=0)
            # Ensure index is integer (frame numbers)
            summary_df.index = summary_df.index.astype(int)
            # Apply binary conversion if cutoff is provided (even for loaded CSV)
            if cutoff is not None:
                console.print(f"[dim]Applying binary conversion with cutoff={cutoff} to loaded summary[/dim]")
                cluster_columns = [col for col in summary_df.columns if col.startswith('cluster_')]
                for col in cluster_columns:
                    summary_df[col] = (summary_df[col] >= cutoff).astype(int)
            console.print(f"✅ [green]Loaded summary with {len(summary_df)} windows and {len(summary_df.columns)} clusters[/green]")
        else:
            # Read trajectory if needed
            if not csv:
                console.print(f"\n🔍 [bold]Reading trajectory...[/bold]")
                trajectory.read_trajectory(
                    cutoff_distance=cutoff_distance,
                    max_frames=max_frames
                )
                console.print(f"✅ [green]Trajectory read and saved[/green]")
            
            # Generate summary
            console.print(f"\n📊 [bold]Generating trajectory summary...[/bold]")
            if cutoff is not None:
                console.print(f"[dim]Using binary mode with cutoff={cutoff}[/dim]")
            trajectory.summarize_trajectory(
                window_size=window_size,
                cutoff=cutoff
            )
            
            # Load summary for display
            base_name = os.path.splitext(os.path.basename(input_file_path))[0]
            output_dir = os.path.dirname(input_file_path) or '.'
            if cutoff is not None:
                summary_csv_path = os.path.join(output_dir, f"{base_name}_summary_binary_{cutoff}.csv")
            else:
                summary_csv_path = os.path.join(output_dir, f"{base_name}_summary.csv")
            summary_df = pd.read_csv(summary_csv_path, index_col=0)
            console.print(f"✅ [green]Generated summary with {len(summary_df)} windows[/green]")
        
        # Create plot
        console.print(f"\n🎨 [bold]Creating summary plot...[/bold]")
        trajectory.plot_summary(
            cutoff=cutoff,
            output_format=output_format
        )
        
        # Determine plot output path
        base_name = os.path.splitext(os.path.basename(input_file_path))[0]
        plot_path = f'results/plots/{base_name}_summary.{output_format}'
        
        console.print(f"\n✅ [bold green]Plot complete![/bold green]")
        console.print(f"💾 Plot saved to: [cyan]{plot_path}[/cyan]")
        console.print(f"📊 Summary shape: [cyan]{summary_df.shape}[/cyan] (windows x clusters)")
        console.print(f"📈 Number of windows: [cyan]{len(summary_df)}[/cyan]")
        console.print(f"🔢 Number of clusters: [cyan]{len(summary_df.columns)}[/cyan]")
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@trajectory.command('classify')
@click.argument('input_file_path', type=click.Path(exists=True, readable=True))
@click.option('--contacts', '-c',
              type=click.Path(exists=True, readable=True),
              default='data/BDD_clustered.csv',
              help='Path to native contacts CSV file (must contain columns: i, j, r, cluster). Default: data/BDD_clustered.csv')
@click.option('--summary-csv',
              type=click.Path(exists=True, readable=True),
              default=None,
              help='Path to binary summary CSV file (if not provided, will generate from trajectory)')
@click.option('--csv',
              type=click.Path(exists=True, readable=True),
              default=None,
              help='Path to CSV file with parsed trajectory data (only used if summary-csv is not provided)')
@click.option('--cutoff-distance',
              type=float,
              default=1.2,
              help='Multiplier for native distance cutoff (default: 1.2, only used if generating summary)')
@click.option('--max-frames',
              type=int,
              default=None,
              help='Maximum number of frames to process (only used if generating summary, default: all frames)')
@click.option('--window-size',
              type=int,
              default=10000,
              help='Number of frames per window (default: 10000, only used if generating summary)')
@click.option('--cutoff',
              type=float,
              default=0.5,
              help='Cutoff for binary classification (default: 0.5, only used if generating summary)')
@click.option('--output',
              type=click.Path(),
              default=None,
              help='Path for output classification file (default: auto-generated as {base}_class.txt)')
def classify_trajectory(input_file_path, contacts, summary_csv, csv, cutoff_distance, max_frames, window_size, cutoff, output):
    """
    Classify trajectory by determining cluster formation order.
    
    Works backwards from the first time all clusters are formed to determine
    the order in which clusters were last formed (formation order).
    
    INPUT_FILE_PATH: Path to the PDB trajectory file (used as reference)
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]TRAJECTORY CLASSIFICATION[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    config_table.add_row("📁 Input PDB file:", f"[cyan]{input_file_path}[/cyan]")
    if summary_csv:
        config_table.add_row("📄 Summary CSV:", f"[cyan]{summary_csv}[/cyan]")
    else:
        if csv:
            config_table.add_row("📄 Trajectory CSV:", f"[cyan]{csv}[/cyan]")
        else:
            config_table.add_row("⚙️  Cutoff distance:", f"[yellow]{cutoff_distance}[/yellow]")
            if max_frames:
                config_table.add_row("🔢 Max frames:", f"[yellow]{max_frames}[/yellow] [dim](debug mode)[/dim]")
            else:
                config_table.add_row("🔢 Max frames:", "[green]All[/green]")
            config_table.add_row("📊 Window size:", f"[yellow]{window_size}[/yellow]")
        config_table.add_row("🔢 Binary cutoff:", f"[yellow]{cutoff}[/yellow]")
    config_table.add_row("📁 Contacts file:", f"[cyan]{contacts}[/cyan]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Validate contacts CSV file
        if not os.path.exists(contacts):
            raise FileNotFoundError(f"Contacts CSV file not found: {contacts}")
        
        # Check if it's a CSV file
        if not contacts.lower().endswith('.csv'):
            console.print(f"[yellow]Warning:[/yellow] File '{contacts}' does not have .csv extension. "
                        f"Trajectory class expects a CSV file with columns: i, j, r, cluster")
        
        # Create Trajectory object
        console.print("📊 [bold]Loading native contacts from CSV...[/bold]")
        trajectory = Trajectory(input_file_path, contacts)
        console.print(f"✅ [green]Loaded {trajectory.total_contacts} native contacts[/green]")
        
        # Load or generate binary summary
        summary_df = None
        if summary_csv:
            # Load existing binary summary CSV
            console.print(f"\n📄 [bold]Loading binary summary from CSV...[/bold]")
            if not os.path.exists(summary_csv):
                raise FileNotFoundError(f"Summary CSV file not found: {summary_csv}")
            summary_df = pd.read_csv(summary_csv, index_col=0)
            summary_df.index = summary_df.index.astype(int)
            console.print(f"✅ [green]Loaded summary with {len(summary_df)} windows[/green]")
        else:
            # Generate binary summary from trajectory data
            # Read trajectory if needed
            if not csv:
                console.print(f"\n🔍 [bold]Reading trajectory...[/bold]")
                trajectory.read_trajectory(
                    cutoff_distance=cutoff_distance,
                    max_frames=max_frames
                )
                console.print(f"✅ [green]Trajectory read and saved[/green]")
            
            # Generate binary summary
            console.print(f"\n📊 [bold]Generating binary summary...[/bold]")
            trajectory.summarize_trajectory(
                window_size=window_size,
                cutoff=cutoff
            )
            
            # Load summary for classification
            base_name = os.path.splitext(os.path.basename(input_file_path))[0]
            output_dir = os.path.dirname(input_file_path) or '.'
            summary_csv_path = os.path.join(output_dir, f"{base_name}_summary_binary_{cutoff}.csv")
            console.print(f"✅ [green]Generated binary summary[/green]")
        
        # Classify trajectory
        console.print(f"\n🔍 [bold]Classifying trajectory...[/bold]")
        formation_order = trajectory.classify(
            summary_df_or_path=summary_csv_path,
            output_path=output
        )
        
        # Determine output path
        if output:
            output_path = output
        else:
            base_name = os.path.splitext(os.path.basename(input_file_path))[0]
            output_dir = os.path.dirname(input_file_path) or '.'
            output_path = os.path.join(output_dir, f"{base_name}_class.txt")
        
        console.print(f"\n✅ [bold green]Classification complete![/bold green]")
        console.print(f"📝 Formation order: [cyan]{formation_order}[/cyan]")
        console.print(f"💾 Classification saved to: [cyan]{output_path}[/cyan]")
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@trajectory.command('animate')
@click.argument('input_file_path', type=click.Path(exists=True, readable=True))
@click.option('--contacts', '-c',
              type=click.Path(exists=True, readable=True),
              default='data/BDD_clustered.csv',
              help='Path to native contacts CSV file (must contain columns: i, j, r, cluster). Default: data/BDD_clustered.csv')
@click.option('--csv',
              type=click.Path(exists=True, readable=True),
              default=None,
              help='Path to CSV file with parsed trajectory data (if not provided, will read PDB)')
@click.option('--cutoff-distance',
              type=float,
              default=1.2,
              help='Multiplier for native distance cutoff (default: 1.2, only used if reading PDB)')
@click.option('--max-frames',
              type=int,
              default=None,
              help='Maximum number of frames to process (only used if reading PDB, default: all frames)')
@click.option('--interval',
              type=int,
              default=50,
              help='Delay between frames in milliseconds (default: 50)')
@click.option('--output-format',
              type=click.Choice(['gif', 'mp4'], case_sensitive=False),
              default='gif',
              help='Output format: gif or mp4 (default: gif)')
@click.option('--fps',
              type=int,
              default=20,
              help='Frames per second for saved animation (default: 20, only used for mp4)')
def animate_trajectory(input_file_path, contacts, csv, cutoff_distance, max_frames, interval, output_format, fps):
    """
    Create an animated contact map showing which contacts are formed in each frame.
    
    Creates an animation where contacts are shown as colored dots on a contact map,
    with different colors for each cluster. Contacts appear/disappear as they form/break.
    
    INPUT_FILE_PATH: Path to the PDB trajectory file (or reference path if using --csv)
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]TRAJECTORY ANIMATION[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    if csv:
        config_table.add_row("📄 CSV file:", f"[cyan]{csv}[/cyan]")
        config_table.add_row("📁 PDB file (reference):", f"[cyan]{input_file_path}[/cyan]")
    else:
        config_table.add_row("📁 Input PDB file:", f"[cyan]{input_file_path}[/cyan]")
        config_table.add_row("⚙️  Cutoff distance:", f"[yellow]{cutoff_distance}[/yellow]")
        if max_frames:
            config_table.add_row("🔢 Max frames:", f"[yellow]{max_frames}[/yellow] [dim](debug mode)[/dim]")
        else:
            config_table.add_row("🔢 Max frames:", "[green]All[/green]")
    config_table.add_row("📁 Contacts file:", f"[cyan]{contacts}[/cyan]")
    config_table.add_row("⏱️  Interval:", f"[yellow]{interval} ms[/yellow]")
    config_table.add_row("🖼️  Output format:", f"[cyan]{output_format}[/cyan]")
    if output_format == 'mp4':
        config_table.add_row("🎬 FPS:", f"[yellow]{fps}[/yellow]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Validate contacts CSV file
        if not os.path.exists(contacts):
            raise FileNotFoundError(f"Contacts CSV file not found: {contacts}")
        
        # Check if it's a CSV file
        if not contacts.lower().endswith('.csv'):
            console.print(f"[yellow]Warning:[/yellow] File '{contacts}' does not have .csv extension. "
                        f"Trajectory class expects a CSV file with columns: i, j, r, cluster")
        
        # Create Trajectory object
        console.print("📊 [bold]Loading native contacts from CSV...[/bold]")
        trajectory = Trajectory(input_file_path, contacts)
        console.print(f"✅ [green]Loaded {trajectory.total_contacts} native contacts[/green]")
        
        # Load or read data
        # Read trajectory if needed
        if not csv:
            console.print(f"\n🔍 [bold]Reading trajectory...[/bold]")
            trajectory.read_trajectory(
                cutoff_distance=cutoff_distance,
                max_frames=max_frames
            )
            console.print(f"✅ [green]Trajectory read and saved[/green]")
        
        # Determine CSV path for animation
        if csv:
            csv_path = csv
        else:
            base_name = os.path.splitext(os.path.basename(input_file_path))[0]
            output_dir = os.path.dirname(input_file_path) or '.'
            csv_path = os.path.join(output_dir, f"{base_name}_parsed.csv")
        
        # Create animation
        console.print(f"\n🎬 [bold]Creating animation...[/bold]")
        trajectory.animate(
            csv_path=csv_path,
            interval=interval,
            save_animation=True,
            output_format=output_format,
            fps=fps
        )
        
        # Determine output path
        base_name = os.path.splitext(os.path.basename(input_file_path))[0]
        output_path = f'results/plots/{base_name}_contact_animation.{output_format}'
        
        console.print(f"\n✅ [bold green]Animation complete![/bold green]")
        console.print(f"💾 Animation saved to: [cyan]{output_path}[/cyan]")
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@cli.group()
def contacts():
    """Contact map analysis commands."""
    pass


@contacts.command('parse')
@click.argument('input_file', type=click.Path(exists=True, readable=True))
@click.option('--clusters',
              type=int,
              default=10,
              help='Number of clusters to create (default: 10)')
@click.option('--cutoff-size',
              type=int,
              default=5,
              help='Minimum cluster size to retain (default: 5)')
@click.option('--output-name',
              type=str,
              default=None,
              help='Base name for output files (default: auto-generated as contacts_{n}_clusters)')
@click.option('--no-draw',
              is_flag=True,
              default=False,
              help='Do not draw the contact map (default: False, map is drawn by default)')
@click.option('--output-format',
              type=click.Choice(['png', 'svg'], case_sensitive=False),
              default='png',
              help='Output format for the plot (default: png, only used if drawing)')
@click.option('--no-show',
              is_flag=True,
              default=False,
              help='Do not display the plot (only save it, only used if drawing)')
def parse_contacts(input_file, clusters, cutoff_size, output_name, no_draw, output_format, no_show):
    """
    Parse contact file and perform clustering.
    
    Creates a CSV file with contact indices (i, j), r6, r12, distance (r), and cluster number.
    By default, also draws and saves the contact map visualization.
    
    INPUT_FILE: Path to the contact file (whitespace-separated table with columns: i, j, r6, r12)
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]CONTACT MAP PARSING[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    config_table.add_row("📁 Input file:", f"[cyan]{input_file}[/cyan]")
    config_table.add_row("🔢 Number of clusters:", f"[yellow]{clusters}[/yellow]")
    config_table.add_row("✂️  Cutoff size:", f"[yellow]{cutoff_size}[/yellow]")
    if output_name:
        config_table.add_row("📝 Output base name:", f"[cyan]{output_name}[/cyan]")
    config_table.add_row("🎨 Draw map:", f"[cyan]{'Yes' if not no_draw else 'No'}[/cyan]")
    if not no_draw:
        config_table.add_row("🖼️  Output format:", f"[cyan]{output_format}[/cyan]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Load contact map and perform clustering (drawing is integrated)
        console.print("📊 [bold]Loading contact map and performing clustering...[/bold]")
        contact_map = ContactMap(
            input_file=input_file,
            cluster_number=clusters,
            cutoff_size=cutoff_size,
            draw_map=not no_draw,
            output_base_name=output_name,
            output_format=output_format,
            show_plot=not no_show
        )
        console.print(f"✅ [green]Loaded {len(contact_map)} contacts[/green]")
        console.print(f"✅ [green]Clustering complete![/green]")
        
        # Get the dataframe with cluster assignments
        df = contact_map.df
        
        # Ensure we have the required columns
        required_columns = ['i', 'j', 'r6', 'r12', 'r', 'cluster']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"Missing required columns. Expected: {required_columns}, got: {list(df.columns)}")
        
        # Select only the required columns in the specified order
        output_df = df[required_columns].copy()
        
        # Display summary
        console.print(f"\n✅ [bold green]Parsing complete![/bold green]")
        console.print(f"📊 Total contacts: [cyan]{len(output_df)}[/cyan]")
        console.print(f"🔢 Unique clusters: [cyan]{output_df['cluster'].nunique()}[/cyan]")
        
        # Cluster distribution
        cluster_dist = output_df['cluster'].value_counts().sort_index()
        stats_table = Table(title="Cluster Distribution", box=box.ROUNDED, show_header=True)
        stats_table.add_column("Cluster", style="cyan", no_wrap=True)
        stats_table.add_column("Count", style="yellow")
        
        for cluster_num, count in cluster_dist.items():
            cluster_label = "Unassigned" if cluster_num == 0 else f"Cluster {cluster_num}"
            stats_table.add_row(cluster_label, str(count))
        
        console.print(f"\n📋 [bold]Cluster distribution:[/bold]")
        console.print(stats_table)
        
        # Display output paths
        csv_path = f'results/{contact_map.output_base_name}.csv'
        console.print(f"\n💾 [bold]Output files:[/bold]")
        console.print(f"📄 CSV: [cyan]{csv_path}[/cyan]")
        if not no_draw:
            plot_path = f'results/plots/{contact_map.output_base_name}.{output_format}'
            console.print(f"🖼️  Plot: [cyan]{plot_path}[/cyan]")
        
        # Display first few rows
        console.print(f"\n📄 [bold]First few rows:[/bold]")
        console.print(output_df.head().to_string())
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@contacts.command('draw')
@click.argument('csv_file', type=click.Path(exists=True, readable=True))
@click.option('--output-name',
              type=str,
              default=None,
              help='Base name for output plot file (default: auto-generated from CSV name)')
@click.option('--output-format',
              type=click.Choice(['png', 'svg'], case_sensitive=False),
              default='png',
              help='Output format for the plot (default: png)')
@click.option('--no-show',
              is_flag=True,
              default=False,
              help='Do not display the plot (only save it)')
def draw_contacts(csv_file, output_name, output_format, no_show):
    """
    Draw a contact map visualization from an existing clustered CSV file.
    
    Creates a contact map showing contacts as colored rectangles below the diagonal,
    with different colors for each cluster. Contacts with cluster=0 are shown in black.
    
    CSV_FILE: Path to CSV file with clustered contacts (must contain columns: i, j, r6, r12, r, cluster)
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]CONTACT MAP VISUALIZATION[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    config_table.add_row("📄 CSV file:", f"[cyan]{csv_file}[/cyan]")
    if output_name:
        config_table.add_row("📝 Output base name:", f"[cyan]{output_name}[/cyan]")
    config_table.add_row("🖼️  Output format:", f"[cyan]{output_format}[/cyan]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Load from CSV
        console.print("📄 [bold]Loading clustered contacts from CSV...[/bold]")
        df = pd.read_csv(csv_file)
        
        # Check if required columns exist
        required_columns = ['i', 'j', 'r6', 'r12', 'r', 'cluster']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"CSV file must contain columns: {required_columns}. Got: {list(df.columns)}")
        
        # Create a ContactMap instance for drawing from existing CSV
        # We need to create a dummy input file temporarily since ContactMap.__init__ requires it
        # The dataframe will be replaced with the one from CSV
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cont', delete=False) as tmp_file:
            # Write minimal contact file (required by ContactMap.__init__)
            tmp_file.write("1 2 1.0 1.0\n")
            tmp_path = tmp_file.name
        
        try:
            # Create ContactMap with dummy file (drawing disabled, we'll do it manually)
            contact_map_obj = ContactMap(
                input_file=tmp_path,
                cluster_number=10,
                cutoff_size=5,
                draw_map=False
            )
            # Replace the dataframe with the one from CSV
            contact_map_obj.df = df
            # Set output base name
            if output_name:
                contact_map_obj.output_base_name = output_name
            else:
                base_name = os.path.splitext(os.path.basename(csv_file))[0]
                contact_map_obj.output_base_name = base_name
            
            console.print(f"✅ [green]Loaded {len(df)} contacts from CSV[/green]")
            
            # Display cluster information
            cluster_dist = df['cluster'].value_counts().sort_index()
            console.print(f"\n📊 Total contacts: [cyan]{len(df)}[/cyan]")
            console.print(f"🔢 Unique clusters: [cyan]{df['cluster'].nunique()}[/cyan]")
            
            # Create visualization
            console.print(f"\n🎨 [bold]Creating contact map visualization...[/bold]")
            
            contact_map_obj.draw(
                save_plot=True,
                output_format=output_format,
                show_plot=not no_show,
                output_base_name=contact_map_obj.output_base_name
            )
            
            # Determine output path
            plot_path = f'results/plots/{contact_map_obj.output_base_name}.{output_format}'
            
            console.print(f"\n✅ [bold green]Visualization complete![/bold green]")
            console.print(f"💾 Plot saved to: [cyan]{plot_path}[/cyan]")
        finally:
            # Clean up temp file
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


if __name__ == '__main__':
    cli()

