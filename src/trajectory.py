import pandas as pd
import random
import numpy as np
import ast
from typing import Dict, List, Tuple, Optional
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

    def parse(self,
              cutoff_distance: float = 1.2,
              max_frames: int = None,
              save_csv: bool = True,
              output_csv_path: str = None) -> pd.DataFrame:
        """
        Parse trajectory and calculate native contact formation.
        
        For each model/frame in the PDB file, calculates distances between contact pairs
        and determines which contacts exist based on cutoff distance.
        
        Args:
            cutoff_distance (float): Multiplier for native distance 'r' (default: 1.2)
            max_frames (int): Maximum number of frames to process (None for all frames, useful for debugging)
            save_csv (bool): Whether to save results to CSV file (default: True)
            output_csv_path (str): Path for output CSV file. If None, auto-generates from input path
            
        Returns:
            pd.DataFrame: DataFrame with columns:
                - 'frame': frame number (int)
                - 'contacts': number of existing contacts (int)
                - 'q': fraction of native contacts existing (float)
                - 'contact_list': list of existing contact pairs (List[Tuple[int, int]])
                - 'clusters_filling': dict mapping cluster number to fraction (Dict[int, float])
        """
        # Read PDB file using raw file reading (much faster than BioPython)
        # Pass max_frames to stop reading early if limit is reached
        frames_data = self._read_pdb_manual(max_frames=max_frames)
        
        if not frames_data:
            return pd.DataFrame(columns=['frame', 'contacts', 'q', 'contact_list', 'clusters_filling'])
        
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
                base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
                output_dir = os.path.dirname(self.trajectory_file_path) or '.'
                output_csv_path = os.path.join(output_dir, f"{base_name}_contacts.csv")
            
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
        
        return df
    
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
            df (pd.DataFrame): DataFrame from parse() method. If None, will call parse() automatically.
            window_size (int): Window size for running average (default: 100)
            figsize (Tuple[int, int]): Figure size (default: (16, 6))
            save_plot (bool): Whether to save the plot (default: True)
            output_format (str): Output format - 'png' or 'svg' (default: 'png')
            show_plot (bool): Whether to display the plot (default: True)
            raw_values (bool): Whether to plot the raw values (default: False)
        """
        # If no DataFrame provided, parse the trajectory
        if df is None:
            df = self.parse(save_csv=False)
        
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
            
            # Plot raw data (light) and running average (bold)
            if raw_values:
                ax2.plot(frames, cluster_values, alpha=0.3, color=colors[idx], linewidth=0.5)
            ax2.plot(frames, cluster_running_avg, color=colors[idx], 
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
        # If no DataFrame provided, parse the trajectory
        if df is None:
            df = self.parse(save_csv=False)
        
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