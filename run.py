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


@trajectory.command('analyze')
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
@click.option('--save-csv/--no-save-csv',
              default=True,
              help='Save results to CSV file (default: True)')
@click.option('--output-csv',
              type=click.Path(),
              default=None,
              help='Path for output CSV file (default: auto-generated from input path)')
def analyze(input_file_path, contacts, cutoff_distance, max_frames, save_csv, output_csv):
    """
    Analyze trajectory PDB file and calculate native contact formation.
    
    INPUT_FILE_PATH: Path to the PDB trajectory file to analyze
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]TRAJECTORY ANALYSIS[/bold cyan]",
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
    config_table.add_row("💾 Save CSV:", f"[green]{'Yes' if save_csv else 'No'}[/green]")
    if output_csv:
        config_table.add_row("📄 Output CSV:", f"[cyan]{output_csv}[/cyan]")
    
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
        
        # Determine output CSV path (same logic as in Trajectory.parse)
        actual_output_csv = output_csv
        if save_csv and actual_output_csv is None:
            base_name = os.path.splitext(os.path.basename(input_file_path))[0]
            output_dir = os.path.dirname(input_file_path) or '.'
            actual_output_csv = os.path.join(output_dir, f"{base_name}_contacts.csv")
        
        # Create Trajectory object (it will load the CSV internally)
        console.print("📊 [bold]Loading native contacts from CSV...[/bold]")
        trajectory = Trajectory(input_file_path, contacts)
        console.print(f"✅ [green]Loaded {trajectory.total_contacts} native contacts[/green]")
        
        # Parse trajectory
        console.print(f"\n🔍 [bold]Analyzing trajectory...[/bold]")
        results_df = trajectory.parse(
            cutoff_distance=cutoff_distance,
            max_frames=max_frames,
            save_csv=save_csv,
            output_csv_path=output_csv
        )
        
        # Display summary
        console.print(f"\n✅ [bold green]Analysis complete![/bold green]")
        console.print(f"📊 Processed [cyan]{len(results_df)}[/cyan] frames")
        console.print(f"📈 DataFrame shape: [cyan]{results_df.shape}[/cyan]")
        
        if not results_df.empty:
            # Summary statistics table
            stats_table = Table(title="Summary Statistics", box=box.ROUNDED, show_header=True)
            stats_table.add_column("Metric", style="cyan", no_wrap=True)
            stats_table.add_column("Value", style="yellow")
            
            stats_table.add_row("Q range", f"{results_df['q'].min():.4f} - {results_df['q'].max():.4f}")
            stats_table.add_row("Average Q", f"{results_df['q'].mean():.4f}")
            stats_table.add_row("Contacts range", f"{results_df['contacts'].min()} - {results_df['contacts'].max()}")
            stats_table.add_row("Average contacts", f"{results_df['contacts'].mean():.1f}")
            
            console.print(f"\n📋 [bold]Summary statistics:[/bold]")
            console.print(stats_table)
            
            # Display first few rows
            console.print(f"\n📄 [bold]First few rows:[/bold]")
            console.print(results_df.head().to_string())
            
            if save_csv:
                console.print(f"\n✅ [green]Results saved to:[/green] [cyan]{actual_output_csv}[/cyan]")
        
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
              help='Path to CSV file with parsed trajectory data (if not provided, will parse PDB file)')
@click.option('--cutoff-distance', 
              type=float,
              default=1.2,
              help='Multiplier for native distance cutoff (default: 1.2, only used if parsing PDB)')
@click.option('--max-frames',
              type=int,
              default=None,
              help='Maximum number of frames to process (for debugging, default: all frames, only used if parsing PDB)')
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
def draw(input_file_path, contacts, csv, cutoff_distance, max_frames, window_size, output_format, no_show):
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
        
        # Load or parse data
        df = None
        if csv:
            console.print(f"\n📄 [bold]Loading trajectory data from CSV...[/bold]")
            df = pd.read_csv(csv)
            # Convert string representations back to dicts/lists if needed
            if 'clusters_filling' in df.columns and df['clusters_filling'].dtype == 'object':
                # Check if it's string representation
                if isinstance(df['clusters_filling'].iloc[0], str):
                    df['clusters_filling'] = df['clusters_filling'].apply(ast.literal_eval)
            if 'contact_list' in df.columns and df['contact_list'].dtype == 'object':
                if isinstance(df['contact_list'].iloc[0], str):
                    df['contact_list'] = df['contact_list'].apply(ast.literal_eval)
            console.print(f"✅ [green]Loaded {len(df)} frames from CSV[/green]")
        else:
            console.print(f"\n🔍 [bold]Parsing trajectory...[/bold]")
            df = trajectory.parse(
                cutoff_distance=cutoff_distance,
                max_frames=max_frames,
                save_csv=False
            )
            console.print(f"✅ [green]Parsed {len(df)} frames[/green]")
        
        # Create visualization
        console.print(f"\n🎨 [bold]Creating visualization...[/bold]")
        
        trajectory.draw(
            df=df,
            window_size=window_size,
            save_plot=True,
            output_format=output_format,
            show_plot=not no_show
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
@click.option('--output-csv',
              type=click.Path(),
              default=None,
              help='Path for output CSV file (default: auto-generated from input path)')
def parse_contacts(input_file, clusters, cutoff_size, output_csv):
    """
    Parse contact file and perform clustering.
    
    Creates a CSV file with contact indices (i, j), r6, r12, distance (r), and cluster number.
    
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
    
    # Determine output CSV path
    actual_output_csv = output_csv
    if actual_output_csv is None:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_dir = os.path.dirname(input_file) or '.'
        actual_output_csv = os.path.join(output_dir, f"{base_name}_clustered.csv")
    
    config_table.add_row("📄 Output CSV:", f"[cyan]{actual_output_csv}[/cyan]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Load contact map
        console.print("📊 [bold]Loading contact map...[/bold]")
        contact_map = ContactMap(input_file)
        console.print(f"✅ [green]Loaded {len(contact_map)} contacts[/green]")
        
        # Perform clustering
        console.print(f"\n🔍 [bold]Performing clustering...[/bold]")
        contact_map.cluster(cluster_number=clusters, cutoff_size=cutoff_size)
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
        
        # Save to CSV
        console.print(f"\n💾 [bold]Saving to CSV...[/bold]")
        # Ensure output directory exists
        output_dir = os.path.dirname(actual_output_csv)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        
        output_df.to_csv(actual_output_csv, index=False)
        console.print(f"✅ [green]Results saved to:[/green] [cyan]{actual_output_csv}[/cyan]")
        
        # Display first few rows
        console.print(f"\n📄 [bold]First few rows:[/bold]")
        console.print(output_df.head().to_string())
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


@contacts.command('draw')
@click.argument('input_file', type=click.Path(exists=True, readable=True))
@click.option('--csv',
              type=click.Path(exists=True, readable=True),
              default=None,
              help='Path to CSV file with clustered contacts (if not provided, will cluster from input file)')
@click.option('--clusters',
              type=int,
              default=10,
              help='Number of clusters to create (default: 10, only used if not using --csv)')
@click.option('--cutoff-size',
              type=int,
              default=5,
              help='Minimum cluster size to retain (default: 5, only used if not using --csv)')
@click.option('--output-format',
              type=click.Choice(['png', 'svg'], case_sensitive=False),
              default='png',
              help='Output format for the plot (default: png)')
@click.option('--no-show',
              is_flag=True,
              default=False,
              help='Do not display the plot (only save it)')
def draw_contacts(input_file, csv, clusters, cutoff_size, output_format, no_show):
    """
    Draw a contact map visualization with cluster assignments.
    
    Creates a contact map showing contacts as colored rectangles below the diagonal,
    with different colors for each cluster. Contacts with cluster=0 are shown in black.
    
    INPUT_FILE: Path to the contact file (or reference path if using --csv)
    """
    console = Console()
    
    # Display header
    console.print(Panel.fit(
        "[bold cyan]CONTACT MAP VISUALIZATION[/bold cyan]",
        border_style="cyan"
    ))
    
    # Display configuration
    config_table = Table(show_header=False, box=None, padding=(0, 2))
    if csv:
        config_table.add_row("📄 CSV file:", f"[cyan]{csv}[/cyan]")
        config_table.add_row("📁 Contact file (reference):", f"[cyan]{input_file}[/cyan]")
    else:
        config_table.add_row("📁 Input file:", f"[cyan]{input_file}[/cyan]")
        config_table.add_row("🔢 Number of clusters:", f"[yellow]{clusters}[/yellow]")
        config_table.add_row("✂️  Cutoff size:", f"[yellow]{cutoff_size}[/yellow]")
    config_table.add_row("🖼️  Output format:", f"[cyan]{output_format}[/cyan]")
    
    console.print(config_table)
    console.print()
    
    try:
        # Load contact map
        if csv:
            console.print("📄 [bold]Loading clustered contacts from CSV...[/bold]")
            # Load from CSV
            df = pd.read_csv(csv)
            # Check if cluster column exists
            if 'cluster' not in df.columns:
                raise ValueError("CSV file must contain a 'cluster' column. Use 'contacts parse' first or provide a clustered CSV.")
            
            # Create ContactMap from original file to get the structure
            contact_map = ContactMap(input_file)
            # Replace the dataframe with the one from CSV
            contact_map.df = df
            console.print(f"✅ [green]Loaded {len(contact_map)} contacts from CSV[/green]")
        else:
            console.print("📊 [bold]Loading contact map...[/bold]")
            contact_map = ContactMap(input_file)
            console.print(f"✅ [green]Loaded {len(contact_map)} contacts[/green]")
            
            # Check if clustering has been done
            if 'cluster' not in contact_map.df.columns:
                console.print(f"\n🔍 [bold]Performing clustering...[/bold]")
                contact_map.cluster(cluster_number=clusters, cutoff_size=cutoff_size)
                console.print(f"✅ [green]Clustering complete![/green]")
            else:
                console.print(f"✅ [green]Using existing cluster assignments[/green]")
        
        # Display cluster information
        cluster_dist = contact_map.df['cluster'].value_counts().sort_index()
        console.print(f"\n📊 Total contacts: [cyan]{len(contact_map.df)}[/cyan]")
        console.print(f"🔢 Unique clusters: [cyan]{contact_map.df['cluster'].nunique()}[/cyan]")
        
        # Create visualization
        console.print(f"\n🎨 [bold]Creating contact map visualization...[/bold]")
        
        contact_map.draw(
            save_plot=True,
            output_format=output_format,
            show_plot=not no_show
        )
        
        # Determine output path
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_path = f'results/plots/{base_name}_contact_map.{output_format}'
        
        console.print(f"\n✅ [bold green]Visualization complete![/bold green]")
        console.print(f"💾 Plot saved to: [cyan]{output_path}[/cyan]")
        
    except Exception as e:
        console.print(f"\n❌ [bold red]Error:[/bold red] {e}")
        import traceback
        console.print_exception()
        sys.exit(1)


if __name__ == '__main__':
    cli()

