import pandas as pd
import random
import numpy as np
import ast
from typing import Dict, List, Tuple, Optional, Union
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from utils import get_file_lines


class Trajectory:
    """
    Class for analyzing trajectory PDB files and calculating native contact formation.
    """
    
    def __init__(self, trajectory_file_path: str, native_contacts_csv: str):
        """
        Initialize Trajectory object.
        
        Args:
            trajectory_file_path (str): Path to PDB file (can be multi-model)
            native_contacts_csv (str): Path to CSV file with native contacts.
                Must contain columns: 'i', 'j', 'r', 'cluster'
        """
        self.trajectory_file_path = trajectory_file_path
        
        # Load native contacts from CSV
        if not os.path.exists(native_contacts_csv):
            raise FileNotFoundError(f"Native contacts CSV file not found: {native_contacts_csv}")
        
        self.native_contacts = pd.read_csv(native_contacts_csv)
        
        # Validate required columns
        required_columns = ['i', 'j', 'r', 'cluster']
        missing_columns = [col for col in required_columns if col not in self.native_contacts.columns]
        if missing_columns:
            raise ValueError(f"Native contacts CSV is missing required columns: {missing_columns}. "
                           f"Found columns: {list(self.native_contacts.columns)}")
        
        self.total_contacts = len(self.native_contacts)
        
        # Pre-compute cluster information
        self.cluster_contacts = {}
        for _, row in self.native_contacts.iterrows():
            cluster_num = int(row['cluster'])
            if cluster_num not in self.cluster_contacts:
                self.cluster_contacts[cluster_num] = []
            self.cluster_contacts[cluster_num].append((int(row['i']), int(row['j'])))
        
        self.cluster_sizes = {cluster: len(contacts) for cluster, contacts in self.cluster_contacts.items()}
    
    def _read_pdb_manual(self, max_frames: int = None) -> Dict[int, Dict[int, np.ndarray]]:
        """
        Read PDB file manually and extract CA atom coordinates.
        
        Args:
            max_frames (int): Maximum number of frames to read (None for all frames)
        
        Returns:
            Dict mapping frame number to dict mapping residue number to CA coordinates
        """
        # First, count total lines for progress bar
        total_lines = get_file_lines(self.trajectory_file_path)
        frames_data = {}
        current_model = None
        current_residue_coords = {}
        model_found = False
        model_saved = False  # Track if current model has been saved
        
        # Read file with progress bar
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TextColumn("•"),
            TextColumn("Lines: {task.completed}/{task.total}"),
            TextColumn("•"),
            TextColumn("Models: {task.fields[models]}"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            transient=False
        ) as progress:
            task = progress.add_task(
                f"[cyan]Reading PDB file: {os.path.basename(self.trajectory_file_path)}",
                total=total_lines,
                models=0
            )
            
            with open(self.trajectory_file_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    progress.update(task, advance=1)
                    
                    if line.startswith('MODEL'):
                        # Start of new model
                        model_found = True
                        # Only save previous model if it hasn't been saved yet (no ENDMDL encountered)
                        if current_model is not None and not model_saved:
                            frames_data[current_model] = current_residue_coords
                            progress.update(task, models=len(frames_data))
                        
                        # Check if we've reached max_frames limit
                        if max_frames is not None and len(frames_data) >= max_frames:
                            break
                        
                        current_model = int(line.split()[1])
                        current_residue_coords = {}
                        model_saved = False  # New model, not saved yet
                        progress.update(task, models=len(frames_data) + 1)
                    
                    elif line.startswith('ATOM') and ' CA ' in line:
                        # Extract residue number and coordinates
                        residue_num = int(line[22:26].strip())
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        
                        current_residue_coords[residue_num] = np.array([x, y, z])
                    
                    elif line.startswith('ENDMDL'):
                        # End of model - save the data
                        if current_model is not None and not model_saved:
                            frames_data[current_model] = current_residue_coords
                            model_saved = True
                            progress.update(task, models=len(frames_data))
                            # Reset for next model (if any)
                            current_residue_coords = {}
                        
                        # Check if we've reached max_frames limit
                        if max_frames is not None and len(frames_data) >= max_frames:
                            break
            
            # Handle last model if file doesn't end with ENDMDL (only if we haven't reached limit)
            if max_frames is None or len(frames_data) < max_frames:
                if current_model is not None and not model_saved and current_residue_coords:
                    frames_data[current_model] = current_residue_coords
                    progress.update(task, models=len(frames_data))
                # Handle single-model files without MODEL/ENDMDL markers
                elif not model_found and current_residue_coords:
                    frames_data[1] = current_residue_coords
                    progress.update(task, models=1)
        
        return frames_data

    def read_trajectory(self,
                        cutoff_distance: float = 1.2,
                        max_frames: Optional[int] = None,
                        save_csv: bool = True,
                        output_csv_path: Optional[str] = None,
                        return_df: bool = False) -> Optional[pd.DataFrame]:
        """
        Read trajectory file and calculate native contact formation.
    
        For each model/frame in the PDB file, calculates distances between contact pairs
        and determines which contacts exist based on cutoff distance.
        
        Args:
            cutoff_distance (float): Multiplier for native distance 'r' (default: 1.2)
                max_frames (int, optional): Maximum number of frames to process (None for all frames, useful for debugging)
                save_csv (bool): Whether to save results to CSV file (default: True)
                output_csv_path (str, optional): Path for output CSV file. If None, auto-generates from input path
                return_df (bool): Whether to return the DataFrame (default: False)
            
        Returns:
                Optional[pd.DataFrame]: DataFrame with columns:
                - 'frame': frame number (int)
                - 'contacts': number of existing contacts (int)
                - 'q': fraction of native contacts existing (float)
                - 'contact_list': list of existing contact pairs (List[Tuple[int, int]])
                - 'clusters_filling': dict mapping cluster number to fraction (Dict[int, float])
                    Returns None if return_df=False, otherwise returns the DataFrame
        """
        # Read PDB file using raw file reading (much faster than BioPython)
        # Pass max_frames to stop reading early if limit is reached
        frames_data = self._read_pdb_manual(max_frames=max_frames)

        if not frames_data:
            return None
    
        results = []
        sorted_frame_nums = sorted(frames_data.keys())
        total_frames = len(sorted_frame_nums)
        
        # Process each frame with progress bar
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TextColumn("•"),
            TextColumn("Frame: {task.completed}/{task.total}"),
            TextColumn("•"),
            TextColumn("Contacts: {task.fields[contacts]}"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            transient=False
        ) as progress:
            task = progress.add_task(
                "[cyan]Processing frames and calculating contacts",
                total=total_frames,
                contacts=0
            )
            
            for frame_num in sorted_frame_nums:
                residue_coords = frames_data[frame_num]
                
                existing_contacts = []
                cluster_counts = {cluster: 0 for cluster in self.cluster_contacts.keys()}
                
                # Check each native contact
                for _, row in self.native_contacts.iterrows():
                    i = int(row['i'])
                    j = int(row['j'])
                    r_native = row['r']
                    cluster_num = int(row['cluster'])
                    
                    # Check if both residues exist in the structure
                    if i in residue_coords and j in residue_coords:
                        # Calculate distance between CA atoms
                        coord_i = residue_coords[i]
                        coord_j = residue_coords[j]
                        distance = np.linalg.norm(coord_i - coord_j)
                        
                        # Check if contact exists (distance < r_native * cutoff_distance)
                        if distance < r_native * cutoff_distance:
                            existing_contacts.append((i, j))
                            cluster_counts[cluster_num] += 1
                
                # Calculate q (fraction of native contacts existing)
                q = len(existing_contacts) / self.total_contacts if self.total_contacts > 0 else 0.0
                
                # Calculate clusters_filling (fraction of contacts in each cluster)
                clusters_filling = {}
                for cluster_num, count in cluster_counts.items():
                    cluster_size = self.cluster_sizes[cluster_num]
                    if cluster_size > 0:
                        clusters_filling[cluster_num] = count / cluster_size
                    else:
                        clusters_filling[cluster_num] = 0.0
                
                results.append({
                    'frame': frame_num,
                    'contacts': len(existing_contacts),
                    'q': q,
                    'contact_list': existing_contacts,
                    'clusters_filling': clusters_filling
                })
            
                # Update progress bar
                progress.update(task, advance=1, contacts=len(existing_contacts))
        
        # Convert results to DataFrame
        df = pd.DataFrame(results)
        
        # Sort by frame number
        if not df.empty:
            df = df.sort_values('frame').reset_index(drop=True)
        
        # Save to CSV if requested
        if save_csv and not df.empty:
            if output_csv_path is None:
                # Auto-generate output path from input path
                # Remove .pdb extension and add _parsed.csv suffix
                base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
                output_dir = os.path.dirname(self.trajectory_file_path) or '.'
                output_csv_path = os.path.join(output_dir, f"{base_name}_parsed.csv")
            
            # Ensure output directory exists
            output_dir = os.path.dirname(output_csv_path)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
            
            # Save DataFrame to CSV
            # Note: contact_list and clusters_filling are complex types, so we'll convert them to strings
            df_to_save = df.copy()
            df_to_save['contact_list'] = df_to_save['contact_list'].apply(lambda x: str(x))
            df_to_save['clusters_filling'] = df_to_save['clusters_filling'].apply(lambda x: str(x))
            
            df_to_save.to_csv(output_csv_path, index=False)
        
        if return_df:
            return df
        else:
            return None
    
    def parse(self,
              cutoff_distance: float = 1.2,
              max_frames: Optional[int] = None,
              window_size: int = 10000,
              cutoff: float = 0.5,
              output_format: str = "png",
              show_plot: bool = False,
              animate: bool = False,
              return_summary: bool = False) -> Optional[pd.DataFrame]:
        """
        Main method for parsing trajectory: reads, draws, summarizes, and classifies.
        
        This method orchestrates the complete trajectory analysis pipeline:
        1. Reads trajectory and saves as _parsed.csv
        2. Draws trajectory visualization
        3. Summarizes trajectory (float mode) and saves as _summarized.csv
        4. Plots summary and saves as _summarized.{format}
        5. Summarizes trajectory (binary mode) and saves as _summarized_binary.csv
        6. Plots binary summary and saves as _summarized_binary.{format}
        7. Classifies trajectory and saves as _class.txt
        8. Optionally animates trajectory and saves as .gif
        
        All output files are saved in the same directory as the trajectory file.
        
        Args:
            cutoff_distance (float): Multiplier for native distance cutoff (default: 1.2)
            max_frames (int, optional): Maximum number of frames to process (default: None)
            window_size (int): Window size for summarization (default: 10000)
            cutoff (float): Cutoff for binary classification (default: 0.5)
            output_format (str): Output format for plots - 'png' or 'svg' (default: 'png')
            show_plot (bool): Whether to display plots (default: False)
            animate (bool): Whether to create animation (default: False)
            return_summary (bool): Whether to return the binary summary DataFrame (default: False)
            
        Returns:
            Optional[pd.DataFrame]: Binary summary DataFrame if return_summary=True, else None
        """
        import shutil
        
        # Get trajectory directory for all output files
        traj_dir = os.path.dirname(self.trajectory_file_path) or '.'
        base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
        
        # Step 1: Read trajectory
        df = self.read_trajectory(
            cutoff_distance=cutoff_distance,
            max_frames=max_frames,
            save_csv=True,
            return_df=True
        )
        
        if df is None or df.empty:
            raise ValueError("Failed to read trajectory or trajectory is empty")
        
        # Step 2: Draw trajectory
        self.draw(
            df=df,
            save_plot=True,
            output_format=output_format,
            show_plot=show_plot
        )
        # Move plot to trajectory directory
        old_plot_path = f'results/plots/{base_name}_trajectory_analysis.{output_format}'
        new_plot_path = os.path.join(traj_dir, f"{base_name}_trajectory_analysis.{output_format}")
        if os.path.exists(old_plot_path):
            os.makedirs(traj_dir, exist_ok=True)
            shutil.move(old_plot_path, new_plot_path)
        
        # Step 3: Summarize trajectory (float mode)
        summary_df = self.summarize_trajectory(
            df=df,
            window_size=window_size,
            save_csv=True,
            output_csv_path=os.path.join(traj_dir, f"{base_name}_summarized.csv"),
            cutoff=None
        )
        
        # Step 4: Plot summary (float mode)
        self.plot_summary(
            summary_df=summary_df,
            save_plot=True,
            output_format=output_format,
            show_plot=show_plot,
            cutoff=None
        )
        # Move plot to trajectory directory
        old_summary_plot = f'results/plots/{base_name}_summary.{output_format}'
        new_summary_plot = os.path.join(traj_dir, f"{base_name}_summarized.{output_format}")
        if os.path.exists(old_summary_plot):
            os.makedirs(traj_dir, exist_ok=True)
            shutil.move(old_summary_plot, new_summary_plot)
        
        # Step 5: Summarize trajectory (binary mode)
        binary_summary_df = self.summarize_trajectory(
            df=df,
            window_size=window_size,
            save_csv=True,
            output_csv_path=os.path.join(traj_dir, f"{base_name}_summarized_binary.csv"),
            cutoff=cutoff
        )
        
        # Step 6: Plot binary summary
        self.plot_summary(
            summary_df=binary_summary_df,
            save_plot=True,
            output_format=output_format,
            show_plot=show_plot,
            cutoff=cutoff
        )
        # Move plot to trajectory directory
        old_binary_plot = f'results/plots/{base_name}_summary.{output_format}'
        new_binary_plot = os.path.join(traj_dir, f"{base_name}_summarized_binary.{output_format}")
        if os.path.exists(old_binary_plot):
            os.makedirs(traj_dir, exist_ok=True)
            shutil.move(old_binary_plot, new_binary_plot)
        
        # Step 7: Classify trajectory
        formation_order = self.classify(
            summary_df_or_path=binary_summary_df,
            output_path=os.path.join(traj_dir, f"{base_name}_class.txt")
        )
        
        # Step 8: Optionally animate
        if animate:
            self.animate(
                df=df,
                save_animation=True,
                output_format="gif"
            )
            # Move animation to trajectory directory
            old_anim = f'results/plots/{base_name}_contact_animation.gif'
            new_anim = os.path.join(traj_dir, f"{base_name}_animation.gif")
            if os.path.exists(old_anim):
                os.makedirs(traj_dir, exist_ok=True)
                shutil.move(old_anim, new_anim)
        
        if return_summary:
            return binary_summary_df
        else:
            return None
    
    def draw(self, 
             df: Optional[pd.DataFrame] = None,
             window_size: int = 100,
             figsize: Tuple[int, int] = (16, 6),
             save_plot: bool = True,
             output_format: str = "png",
             show_plot: bool = True,
             raw_values: bool = False) -> None:
        """
        Create a two-panel plot showing q and cluster filling over time.
        
        Args:
            df (pd.DataFrame): DataFrame from read_trajectory() method. If None, will call read_trajectory() automatically.
            window_size (int): Window size for running average (default: 100)
            figsize (Tuple[int, int]): Figure size (default: (16, 6))
            save_plot (bool): Whether to save the plot (default: True)
            output_format (str): Output format - 'png' or 'svg' (default: 'png')
            show_plot (bool): Whether to display the plot (default: True)
            raw_values (bool): Whether to plot the raw values (default: False)
        """
        # If no DataFrame provided, read the trajectory
        if df is None:
            df = self.read_trajectory(save_csv=False, return_df=True)
            if df is None:
                df = pd.DataFrame()
        
        if df.empty:
            raise ValueError("DataFrame is empty. Cannot create plot.")
        
        # Ensure DataFrame is sorted by frame
        df = df.sort_values('frame').reset_index(drop=True)
        
        # Create figure with two subplots side by side
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Left panel: q vs frame
        frames = df['frame'].values
        q_values = df['q'].values
        
        # Calculate running average for q
        q_series = pd.Series(q_values)
        q_running_avg = q_series.rolling(window=window_size, center=True, min_periods=1).mean()
        
        # Plot raw data (light) and running average (bold)
        if raw_values:
            ax1.plot(frames, q_values, alpha=0.3, color='blue', label='Raw', linewidth=0.5)
        ax1.plot(frames, q_running_avg, color='blue', label=f'Running avg (window={window_size})', linewidth=2)
        ax1.set_xlabel('Frame Number', fontsize=12)
        ax1.set_ylabel('Fraction of Native Contacts (q)', fontsize=12)
        ax1.set_title('Native Contact Formation Over Time', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Right panel: clusters_filling vs frame
        # Extract cluster data from clusters_filling column
        # clusters_filling is a dict, so we need to extract each cluster's values
        cluster_numbers = sorted(self.cluster_contacts.keys())
        
        # Create a color map for clusters
        colors = plt.cm.tab10(np.linspace(0, 1, len(cluster_numbers)))
        
        for idx, cluster_num in enumerate(cluster_numbers):
            # Extract cluster filling values for this cluster
            cluster_values = []
            for clusters_dict in df['clusters_filling']:
                if isinstance(clusters_dict, dict):
                    cluster_values.append(clusters_dict.get(cluster_num, 0.0))
                else:
                    # Handle string representation of dict
                    if isinstance(clusters_dict, str):
                        clusters_dict = ast.literal_eval(clusters_dict)
                    cluster_values.append(clusters_dict.get(cluster_num, 0.0))
            
            cluster_series = pd.Series(cluster_values)
            cluster_running_avg = cluster_series.rolling(window=window_size, center=True, min_periods=1).mean()
            
            # Use black for cluster_0, otherwise use color from colormap
            plot_color = 'black' if cluster_num == 0 else colors[idx]
            
            # Plot raw data (light) and running average (bold)
            if raw_values:
                ax2.plot(frames, cluster_values, alpha=0.3, color=plot_color, linewidth=0.5)
            ax2.plot(frames, cluster_running_avg, color=plot_color, 
                    label=f'Cluster {cluster_num} (avg)', linewidth=1)
        
        ax2.set_xlabel('Frame Number', fontsize=12)
        ax2.set_ylabel('Cluster Filling Fraction', fontsize=12)
        ax2.set_title('Cluster Formation Over Time', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='best', fontsize=9)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot if requested
        if save_plot:
            os.makedirs('results/plots', exist_ok=True)
            base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
            filename = f'results/plots/{base_name}_trajectory_analysis.{output_format}'
            plt.savefig(filename, format=output_format, dpi=300, bbox_inches='tight')
        
        if show_plot:
            plt.show()
        else:
            plt.close()
    
    def summarize_trajectory(self,
                            df: Optional[pd.DataFrame] = None,
                            window_size: int = 10000,
                            save_csv: bool = True,
                            output_csv_path: str = None,
                            cutoff: Optional[float] = None) -> pd.DataFrame:
        """
        Summarize trajectory by calculating mean cluster filling fractions in windows.
        
        Groups frames into windows and calculates the mean fraction of contacts
        established for each cluster in each window.
        
        Args:
            df (pd.DataFrame): DataFrame from parse() method. If None, will attempt to load
                from auto-generated CSV path based on trajectory file name.
            window_size (int): Number of frames per window (default: 10000)
            save_csv (bool): Whether to save summary to CSV file (default: True)
            output_csv_path (str): Path for output CSV file. If None, auto-generates from input path
            cutoff (float, optional): If provided, convert probabilities to binary (0 or 1).
                Values >= cutoff become 1, values < cutoff become 0. If None, returns float probabilities.
            
        Returns:
            pd.DataFrame: DataFrame with:
                - Index: frame numbers (integers, starting frame of each window)
                - Columns: cluster numbers (e.g., 'cluster_0', 'cluster_1', etc.)
                - Values: mean fraction of contacts established in that window for that cluster (float),
                  or binary values (0 or 1) if cutoff is provided
        """
        # Load DataFrame if not provided
        if df is None:
            # Try to load from auto-generated CSV path
            base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
            output_dir = os.path.dirname(self.trajectory_file_path) or '.'
            csv_path = os.path.join(output_dir, f"{base_name}_contacts.csv")
            
            if not os.path.exists(csv_path):
                raise FileNotFoundError(
                    f"Trajectory data not found. Please run parse() first or provide a DataFrame. "
                    f"Expected file: {csv_path}"
                )
            
            df = pd.read_csv(csv_path)
            # Convert string representations back to dicts
            if 'clusters_filling' in df.columns:
                df['clusters_filling'] = df['clusters_filling'].apply(
                    lambda x: ast.literal_eval(x) if isinstance(x, str) else x
                )
        
        if df.empty:
            raise ValueError("DataFrame is empty. Cannot create summary.")
        
        # Ensure DataFrame is sorted by frame
        df = df.sort_values('frame').reset_index(drop=True)
        
        # Get all unique cluster numbers
        cluster_numbers = sorted(self.cluster_contacts.keys())
        
        # Group frames into windows
        total_frames = len(df)
        num_windows = (total_frames + window_size - 1) // window_size  # Ceiling division
        
        # Initialize summary data
        summary_data = []
        
        for window_idx in range(num_windows):
            start_idx = window_idx * window_size
            end_idx = min(start_idx + window_size, total_frames)
            
            window_df = df.iloc[start_idx:end_idx]
            
            # Get the frame number for this window (use the first frame in the window)
            frame_number = int(window_df.iloc[0]['frame'])
            
            # Calculate mean fraction for each cluster in this window
            window_summary = {'frame': frame_number}
            
            for cluster_num in cluster_numbers:
                # Extract cluster filling values for this cluster from all frames in window
                cluster_values = []
                for clusters_dict in window_df['clusters_filling']:
                    if isinstance(clusters_dict, dict):
                        cluster_values.append(clusters_dict.get(cluster_num, 0.0))
                    elif isinstance(clusters_dict, str):
                        # Handle string representation
                        clusters_dict = ast.literal_eval(clusters_dict)
                        cluster_values.append(clusters_dict.get(cluster_num, 0.0))
                    else:
                        cluster_values.append(0.0)
                
                # Calculate mean for this cluster in this window
                mean_fraction = np.mean(cluster_values) if cluster_values else 0.0
                window_summary[f'cluster_{cluster_num}'] = mean_fraction
            
            summary_data.append(window_summary)
        
        # Create summary DataFrame
        summary_df = pd.DataFrame(summary_data)
        
        # Set frame as index
        summary_df = summary_df.set_index('frame')
        
        # Apply binary conversion if cutoff is provided
        if cutoff is not None:
            # Convert all cluster columns to binary (0 or 1)
            cluster_columns = [col for col in summary_df.columns if col.startswith('cluster_')]
            for col in cluster_columns:
                summary_df[col] = (summary_df[col] >= cutoff).astype(int)
        
        # Save to CSV if requested
        if save_csv:
            if output_csv_path is None:
                # Auto-generate output path from input path
                base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
                output_dir = os.path.dirname(self.trajectory_file_path) or '.'
                output_csv_path = os.path.join(output_dir, f"{base_name}_summary.csv")
            
            # Ensure output directory exists
            output_dir = os.path.dirname(output_csv_path)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
            
            summary_df.to_csv(output_csv_path)
        
        return summary_df
    
    def plot_summary(self,
                    summary_df: Optional[pd.DataFrame] = None,
                    df: Optional[pd.DataFrame] = None,
                    window_size: int = 10000,
                    figsize: Tuple[int, int] = (14, 8),
                    save_plot: bool = True,
                    output_format: str = "png",
                    show_plot: bool = True,
                    cutoff: Optional[float] = None) -> None:
        """
        Plot each cluster column from the summary DataFrame as a time series.
        
        Args:
            summary_df (pd.DataFrame): Summary DataFrame from summarize_trajectory().
                If None, will call summarize_trajectory() automatically.
            df (pd.DataFrame): DataFrame from parse() method. Required if summary_df is None.
            window_size (int): Number of frames per window (default: 10000, only used if summary_df is None)
            figsize (Tuple[int, int]): Figure size (default: (14, 8))
            save_plot (bool): Whether to save the plot (default: True)
            output_format (str): Output format - 'png' or 'svg' (default: 'png')
            show_plot (bool): Whether to display the plot (default: True)
            cutoff (float, optional): If provided, convert probabilities to binary (0 or 1) before plotting.
                Values >= cutoff become 1, values < cutoff become 0. If None, plots float probabilities.
        """
        # Generate summary if not provided
        if summary_df is None:
            if df is None:
                raise ValueError("Either summary_df or df must be provided")
            summary_df = self.summarize_trajectory(df=df, window_size=window_size, save_csv=True, cutoff=cutoff)
        
        if summary_df.empty:
            raise ValueError("Summary DataFrame is empty. Cannot create plot.")
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Get frame numbers (index) and cluster columns
        frame_numbers = summary_df.index.values
        cluster_columns = [col for col in summary_df.columns if col.startswith('cluster_')]
        
        # Create color map for clusters
        colors = plt.cm.tab10(np.linspace(0, 1, len(cluster_columns)))
        
        # Plot each cluster as a time series
        for idx, cluster_col in enumerate(cluster_columns):
            cluster_values = summary_df[cluster_col].values
            cluster_num = cluster_col.replace('cluster_', '')
            
            # Use black for cluster 0, otherwise use color map
            if cluster_num == '0':
                plot_color = 'black'
            else:
                plot_color = colors[idx]
            
            ax.plot(frame_numbers, cluster_values, 
                   color=plot_color, 
                   label=f'Cluster {cluster_num}',
                   linewidth=2,
                   marker='o',
                   markersize=4,
                   alpha=0.8)
        
        ax.set_xlabel('Frame Number', fontsize=12)
        if cutoff is not None:
            ax.set_ylabel('Cluster State (Binary)', fontsize=12)
            ax.set_title(f'Cluster Formation Over Time (Binary, cutoff={cutoff})', fontsize=14, fontweight='bold')
            # Set y-axis limits for binary data
            ax.set_ylim(-0.1, 1.1)
        else:
            ax.set_ylabel('Mean Cluster Filling Fraction', fontsize=12)
            ax.set_title('Cluster Formation Over Time (Windowed Summary)', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=9, ncol=2)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot if requested
        if save_plot:
            os.makedirs('results/plots', exist_ok=True)
            base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
            filename = f'results/plots/{base_name}_summary.{output_format}'
            plt.savefig(filename, format=output_format, dpi=300, bbox_inches='tight')
        
        if show_plot:
            plt.show()
        else:
            plt.close()
  
    def classify(self,
                summary_df_or_path: Union[pd.DataFrame, str],
                output_path: Optional[str] = None) -> List[int]:
        """
        Classify trajectory by determining cluster formation order.
        
        Works backwards from the first time all clusters are formed to determine
        the order in which clusters were last formed (formation order).
        
        Args:
            summary_df_or_path (Union[pd.DataFrame, str]): DataFrame from summarize_trajectory() method
                or path to CSV file (binary representation with frame index and cluster columns)
            output_path (str, optional): Path to save the formation order. If None, auto-generates
                from trajectory file path as {base_name}_class.txt
            
        Returns:
            List[int]: List of cluster numbers in formation order (last formed first)
        """
        # Load summary DataFrame
        if isinstance(summary_df_or_path, pd.DataFrame):
            summary_df = summary_df_or_path.copy()
            # Ensure index is integer (frame numbers)
            if not isinstance(summary_df.index, pd.RangeIndex):
                summary_df.index = summary_df.index.astype(int)
            # Use trajectory file path for auto-generating output path
            summary_csv_path = self.trajectory_file_path
        else:
            summary_csv_path = summary_df_or_path
            if not os.path.exists(summary_csv_path):
                raise FileNotFoundError(f"Summary CSV file not found: {summary_csv_path}")
            summary_df = pd.read_csv(summary_csv_path, index_col=0)
            # Ensure index is integer (frame numbers)
            summary_df.index = summary_df.index.astype(int)
        
        # Get all cluster columns (excluding cluster_0)
        cluster_columns = [col for col in summary_df.columns if col.startswith('cluster_') and col != 'cluster_0']
        
        if not cluster_columns:
            raise ValueError("No cluster columns found in summary DataFrame (excluding cluster_0)")
        
        # Sort by frame number (ascending)
        summary_df = summary_df.sort_index()
        
        # Step 2: Find first time all clusters (except cluster_0) are formed (all 1)
        all_formed_mask = (summary_df[cluster_columns] == 1).all(axis=1)
        all_formed_frames = summary_df.index[all_formed_mask]
        
        if len(all_formed_frames) == 0:
            # No frame where all clusters are formed - find frame with maximum clusters formed
            cluster_counts = (summary_df[cluster_columns] == 1).sum(axis=1)
            max_clusters = cluster_counts.max()
            
            if max_clusters == 0:
                raise ValueError("No clusters are formed in the trajectory. Cannot determine formation order.")
            
            # Use first frame where maximum number of clusters are formed
            max_cluster_frames = cluster_counts[cluster_counts == max_clusters].index
            start_frame = max_cluster_frames[0]
        else:
            # Use first frame where all clusters are formed
            start_frame = all_formed_frames[0]
        
        # Step 3-5: Go backwards in time and track cluster breaks
        # Find the index of the start_frame
        start_idx = summary_df.index.get_loc(start_frame)
        
        # Track which clusters are currently formed (start with all formed)
        currently_formed = set()
        for col in cluster_columns:
            if summary_df.loc[start_frame, col] == 1:
                cluster_num = int(col.replace('cluster_', ''))
                currently_formed.add(cluster_num)
        
        # List to store order of cluster breaks (going backwards)
        cluster_breaks_order = []
        
        # Go backwards from start_frame
        for idx in range(start_idx, -1, -1):
            current_frame = summary_df.index[idx]
            current_row = summary_df.loc[current_frame]
            
            # Check which clusters are formed at this frame
            formed_at_frame = set()
            for col in cluster_columns:
                if current_row[col] == 1:
                    cluster_num = int(col.replace('cluster_', ''))
                    formed_at_frame.add(cluster_num)
            
            # Find clusters that were formed before but are now broken (1 -> 0)
            broken_clusters = currently_formed - formed_at_frame
            
            # Add broken clusters to the list (in sorted order for consistency)
            if broken_clusters:
                for cluster_num in sorted(broken_clusters):
                    cluster_breaks_order.append(cluster_num)
            
            # Update currently_formed to reflect current state
            currently_formed = formed_at_frame.copy()
            
            # If all clusters (except cluster_0) are broken, we're done
            if not currently_formed:
                break
        
        # Step 6: If we reached the beginning and some clusters are still formed,
        # append them to the end (these were formed before the trajectory started)
        if currently_formed:
            for cluster_num in sorted(currently_formed):
                cluster_breaks_order.append(cluster_num)
        
        # Reverse the order to get formation order (last formed first in breaks = first formed last in formation)
        formation_order = list(reversed(cluster_breaks_order))
        
        # Save to file
        if output_path is None:
            # Auto-generate output path from trajectory file path
            base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
            output_dir = os.path.dirname(self.trajectory_file_path) or '.'
            output_path = os.path.join(output_dir, f"{base_name}_class.txt")
        
        # Ensure output directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        
        # Save as one-line file (comma-separated cluster numbers)
        with open(output_path, 'w') as f:
            f.write(','.join(map(str, formation_order)) + '\n')
        
        return formation_order
        
    def animate(self,
                df: Optional[pd.DataFrame] = None,
                interval: int = 50,
                figsize: Tuple[int, int] = (12, 10),
                save_animation: bool = True,
                output_format: str = "gif",
                fps: int = 20):
        """
        Create an animated contact map showing which contacts are formed in each frame.
        
        Args:
            df (pd.DataFrame): DataFrame from parse() method. If None, will call parse() automatically.
            interval (int): Delay between frames in milliseconds (default: 50)
            figsize (Tuple[int, int]): Figure size (default: (12, 10))
            save_animation (bool): Whether to save the animation (default: True)
            output_format (str): Output format - 'gif' or 'mp4' (default: 'gif')
            fps (int): Frames per second for saved animation (default: 20, only used for mp4)
        """
        # If no DataFrame provided, read the trajectory
        if df is None:
            df = self.read_trajectory(save_csv=False, return_df=True)
            if df is None:
                df = pd.DataFrame()
        
        if df.empty:
            raise ValueError("DataFrame is empty. Cannot create animation.")
        
        # Ensure DataFrame is sorted by frame
        df = df.sort_values('frame').reset_index(drop=True)
        
        # Create mapping from contact (i, j) to cluster number
        contact_to_cluster = {}
        for _, row in self.native_contacts.iterrows():
            i = int(row['i'])
            j = int(row['j'])
            cluster_num = int(row['cluster'])
            contact_to_cluster[(i, j)] = cluster_num
        
        # Get all unique contacts and their clusters
        all_contacts = list(contact_to_cluster.keys())
        cluster_numbers = sorted(self.cluster_contacts.keys())
        
        # Create color map for clusters
        colors = plt.cm.tab10(np.linspace(0, 1, len(cluster_numbers)))
        cluster_color_map = {cluster_num: colors[idx] for idx, cluster_num in enumerate(cluster_numbers)}
        
        # Find max residue number for axis limits
        max_residue = max(max(contact) for contact in all_contacts) if all_contacts else 100
        
        # Create figure and axis
        fig, ax = plt.subplots(figsize=figsize)
        
        # Initialize scatter plots for each cluster (empty initially)
        scatter_plots = {}
        for cluster_num in cluster_numbers:
            scatter_plots[cluster_num] = ax.scatter([], [], 
                                                    c=[cluster_color_map[cluster_num]], 
                                                    s=50, 
                                                    alpha=0.7,
                                                    edgecolors='black',
                                                    linewidths=0.5,
                                                    label=f'Cluster {cluster_num}')
        
        # Set axis properties
        ax.set_xlim(0.5, max_residue + 0.5)
        ax.set_ylim(0.5, max_residue + 0.5)
        ax.set_xlabel('Residue Index (j)', fontsize=12)
        ax.set_ylabel('Residue Index (i)', fontsize=12)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Add diagonal line
        ax.plot([0.5, max_residue + 0.5], [0.5, max_residue + 0.5], 
               'k--', linewidth=1, alpha=0.5, label='Diagonal')
        
        # Add title and frame counter
        title_text = ax.text(0.5, 0.95, '', transform=ax.transAxes, 
                           fontsize=14, fontweight='bold',
                           ha='center', va='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Add legend
        ax.legend(loc='upper right', fontsize=8)
        
        # Prepare data for animation - normalize all contacts to (min, max) format
        def normalize_contact(contact):
            """Normalize contact to (min, max) format for consistent matching."""
            if isinstance(contact, (list, tuple)) and len(contact) == 2:
                return (min(contact[0], contact[1]), max(contact[0], contact[1]))
            return contact
        
        frames_data = []
        for _, row in df.iterrows():
            frame_num = row['frame']
            # Get contact list for this frame
            contact_list = row['contact_list']
            if isinstance(contact_list, str):
                contact_list = ast.literal_eval(contact_list)
            
            # Normalize and convert to set for faster lookup
            existing_contacts = set(normalize_contact(c) for c in contact_list)
            frames_data.append((frame_num, existing_contacts))
        
        def animate_frame(frame_idx):
            """Update the plot for each frame."""
            frame_num, existing_contacts = frames_data[frame_idx]
            
            # Get q value for this frame
            frame_row = df.iloc[frame_idx]
            q_value = frame_row['q'] if 'q' in df.columns else len(existing_contacts) / len(all_contacts) if all_contacts else 0.0
            
            # Clear previous scatter plots
            for cluster_num in cluster_numbers:
                scatter_plots[cluster_num].remove()
            
            # Plot contacts for each cluster
            for cluster_num in cluster_numbers:
                x_coords = []
                y_coords = []
                
                # Get all contacts in this cluster
                cluster_contact_list = self.cluster_contacts[cluster_num]
                
                # Filter to only show contacts that exist in this frame
                for contact in cluster_contact_list:
                    # Normalize contact for matching
                    contact_normalized = normalize_contact(contact)
                    
                    # Check if contact exists in this frame
                    if contact_normalized in existing_contacts:
                        # Get i, j for plotting (ensure i < j for below diagonal)
                        i, j = contact_normalized
                        x_coords.append(j)
                        y_coords.append(i)
                
                # Create new scatter plot for this cluster
                if x_coords:
                    scatter_plots[cluster_num] = ax.scatter(x_coords, y_coords,
                                                           c=[cluster_color_map[cluster_num]],
                                                           s=50,
                                                           alpha=0.7,
                                                           edgecolors='black',
                                                           linewidths=0.5,
                                                           label=f'Cluster {cluster_num}')
                else:
                    # Empty scatter plot
                    scatter_plots[cluster_num] = ax.scatter([], [],
                                                           c=[cluster_color_map[cluster_num]],
                                                           s=50,
                                                           alpha=0.7,
                                                           edgecolors='black',
                                                           linewidths=0.5,
                                                           label=f'Cluster {cluster_num}')
            
            # Update title with frame number and contact count
            num_contacts = len(existing_contacts)
            title_text.set_text(f'Frame {frame_num} | Contacts: {num_contacts} | Q: {q_value:.3f}')
        
        # Create animation
        # Note: blit=False because we're removing and recreating scatter plots
        anim = animation.FuncAnimation(fig, animate_frame, frames=len(frames_data),
                                     interval=interval, blit=False, repeat=True)
        
        # Save animation if requested
        if save_animation:
            os.makedirs('results/plots', exist_ok=True)
            base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
            
            if output_format.lower() == 'gif':
                filename = f'results/plots/{base_name}_contact_animation.gif'
                print(f"Saving animation to {filename}...")
                anim.save(filename, writer='pillow', fps=1000/interval)
            elif output_format.lower() == 'mp4':
                filename = f'results/plots/{base_name}_contact_animation.mp4'
                print(f"Saving animation to {filename}...")
                anim.save(filename, writer='ffmpeg', fps=fps)
        else:
            raise ValueError(f"Unsupported output format: {output_format}. Use 'gif' or 'mp4'.")
        print(f"Animation saved successfully!")
        
        plt.tight_layout()
        plt.show()
        
        return anim