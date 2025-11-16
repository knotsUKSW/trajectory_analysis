import pandas as pd
import numpy as np
import ast
from typing import List, Tuple, Optional, Union
import os
import shutil
import matplotlib.pyplot as plt
import matplotlib.animation as animation

try:
    import folding_analysis_rs
except ImportError:
    raise ImportError(
        "Rust bindings (folding_analysis_rs) are required. "
        "Please install them by running: cd rust && maturin develop --features python"
    )


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
        self.native_contacts_csv_path = native_contacts_csv  # Store path for Rust bindings
        
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
    
    def _get_base_paths(self) -> Tuple[str, str]:
        """
        Get base name and output directory for trajectory file.
    
    Returns:
            Tuple[str, str]: (base_name, output_dir)
        """
        base_name = os.path.splitext(os.path.basename(self.trajectory_file_path))[0]
        output_dir = os.path.dirname(self.trajectory_file_path) or '.'
        return base_name, output_dir
    
    def _parse_clusters_filling(self, clusters_filling_str: Union[str, dict]) -> dict:
        """
        Parse clusters_filling from string representation or return dict as-is.
        
        Args:
            clusters_filling_str: String representation of dict or dict itself
            
        Returns:
            dict: Parsed clusters_filling dictionary
        """
        if isinstance(clusters_filling_str, dict):
            return clusters_filling_str
        elif isinstance(clusters_filling_str, str):
            return ast.literal_eval(clusters_filling_str)
        else:
            return {}

    def read_trajectory(self,
                        cutoff_distance: float = 1.2,
                        max_frames: Optional[int] = None) -> None:
        """
        Read trajectory file and calculate native contact formation.
        
        For each model/frame in the PDB file, calculates distances between contact pairs
        and determines which contacts exist based on cutoff distance.
        Always saves results to {base_name}_parsed.csv in the same directory as the trajectory file.
        
        Args:
            cutoff_distance (float): Multiplier for native distance 'r' (default: 1.2)
            max_frames (int, optional): Maximum number of frames to process (None for all frames, useful for debugging)
        """
        # Determine output path
        base_name, output_dir = self._get_base_paths()
        output_csv_path = os.path.join(output_dir, f"{base_name}_parsed.csv")
        
        # Call Rust function
        folding_analysis_rs.read_trajectory(
            trajectory_file=self.trajectory_file_path,
            contacts_file=self.native_contacts_csv_path,
            cutoff_distance=cutoff_distance,
            max_frames=max_frames,
            output_csv=output_csv_path
        )
    
    def parse(self,
              cutoff_distance: float = 1.2,
              max_frames: Optional[int] = None,
              window_size: int = 10000,
              cutoff: float = 0.5,
              output_format: str = "png",
              animate: bool = False):
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
            animate (bool): Whether to create animation (default: False)
        
        """
        # Get trajectory directory for all output files
        base_name, traj_dir = self._get_base_paths()
        
        # Step 1: Read trajectory
        self.read_trajectory(
            cutoff_distance=cutoff_distance,
            max_frames=max_frames
        )
        
        # Step 2: Draw trajectory (uses default window_size=100, which will call smooth if needed)
        self.draw(output_format=output_format)
        # Move plot to trajectory directory
        old_plot_path = f'results/plots/{base_name}_trajectory_analysis.{output_format}'
        new_plot_path = os.path.join(traj_dir, f"{base_name}_trajectory_analysis.{output_format}")
        if os.path.exists(old_plot_path):
            os.makedirs(traj_dir, exist_ok=True)
            shutil.move(old_plot_path, new_plot_path)
        
        # Step 3: Summarize trajectory (float mode)
        self.summarize_trajectory(
            window_size=window_size,
            cutoff=None
        )
        
        # Step 4: Plot summary (float mode)
        self.plot_summary(
            cutoff=None,
            output_format=output_format
        )
        # Move plot to trajectory directory
        old_summary_plot = f'results/plots/{base_name}_summary.{output_format}'
        new_summary_plot = os.path.join(traj_dir, f"{base_name}_summarized.{output_format}")
        if os.path.exists(old_summary_plot):
            os.makedirs(traj_dir, exist_ok=True)
            shutil.move(old_summary_plot, new_summary_plot)
        
        # Step 5: Summarize trajectory (binary mode)
        self.summarize_trajectory(
            window_size=window_size,
            cutoff=cutoff
        )
        
        # Step 6: Plot binary summary
        self.plot_summary(
            cutoff=cutoff,
            output_format=output_format
        )
        # Move plot to trajectory directory
        old_binary_plot = f'results/plots/{base_name}_summary_binary_{cutoff}.{output_format}'
        new_binary_plot = os.path.join(traj_dir, f"{base_name}_summarized_binary.{output_format}")
        if os.path.exists(old_binary_plot):
            os.makedirs(traj_dir, exist_ok=True)
            shutil.move(old_binary_plot, new_binary_plot)
        
        # Step 7: Classify trajectory
        binary_summary_path = os.path.join(traj_dir, f"{base_name}_summary_binary_{cutoff}.csv")
        formation_order = self.classify(
            summary_df_or_path=binary_summary_path,
            output_path=os.path.join(traj_dir, f"{base_name}_class.txt")
        )
        
        # Step 8: Optionally animate
        if animate:
            parsed_csv_path = os.path.join(traj_dir, f"{base_name}_parsed.csv")
            self.animate(
                csv_path=parsed_csv_path,
                save_animation=True,
                output_format="gif"
            )
            # Move animation to trajectory directory
            old_anim = f'results/plots/{base_name}_contact_animation.gif'
            new_anim = os.path.join(traj_dir, f"{base_name}_animation.gif")
            if os.path.exists(old_anim):
                os.makedirs(traj_dir, exist_ok=True)
                shutil.move(old_anim, new_anim)
    
    def smooth(self,
               window_size: int = 100) -> None:
        """
        Smooth trajectory data by calculating running averages.
        
        Calculates running average of q (fraction of all contacts) and per-cluster
        fraction of contacts formed using a centered window.
        Always saves results to {base_name}_{window_size}_smoothed.csv in the same directory as the trajectory file.
    
    Args:
            window_size (int): Window size for running average (default: 100)
        """
        # Determine output path
        base_name, output_dir = self._get_base_paths()
        output_csv_path = os.path.join(output_dir, f"{base_name}_{window_size}_smoothed.csv")
        
        # Call Rust function
        folding_analysis_rs.smooth(
            trajectory_file=self.trajectory_file_path,
            window_size=window_size,
            output_csv=output_csv_path
        )
    
    def draw(self, 
             window_size: int = 100,
             figsize: Tuple[int, int] = (16, 6),
             output_format: str = "png",
             raw_values: bool = False) -> None:
        """
        Create a two-panel plot showing q and cluster filling over time.
        
        Args:
            window_size (int): Window size for running average (default: 100).
                If window_size > 1, uses smoothed data from {base_name}_{window_size}_smoothed.csv.
                If window_size == 1, uses raw data from {base_name}_parsed.csv.
            figsize (Tuple[int, int]): Figure size (default: (16, 6))
            output_format (str): Output format - 'png' or 'svg' (default: 'png')
            raw_values (bool): Whether to plot the raw values (default: False)
        """
        base_name, output_dir = self._get_base_paths()
        
        # Determine input CSV path based on window_size
        if window_size > 1:
            # First, ensure smoothed data exists
            smoothed_csv_path = os.path.join(output_dir, f"{base_name}_{window_size}_smoothed.csv")
            if not os.path.exists(smoothed_csv_path):
                # Generate smoothed data
                self.smooth(window_size=window_size)
            csv_path = smoothed_csv_path
        else:
            csv_path = os.path.join(output_dir, f"{base_name}_parsed.csv")
        
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Trajectory CSV file not found: {csv_path}")
        
        # Load DataFrame from CSV
        df_plot = pd.read_csv(csv_path)
        
        if df_plot.empty:
            raise ValueError("DataFrame is empty. Cannot create plot.")
        
        # Ensure DataFrame is sorted by frame
        df_plot = df_plot.sort_values('frame').reset_index(drop=True)
        
        # Load raw data if needed for raw_values overlay
        if raw_values and window_size > 1:
            raw_csv_path = os.path.join(output_dir, f"{base_name}_parsed.csv")
            if os.path.exists(raw_csv_path):
                df = pd.read_csv(raw_csv_path)
                if 'clusters_filling' in df.columns:
                    df['clusters_filling'] = df['clusters_filling'].apply(self._parse_clusters_filling)
                df = df.sort_values('frame').reset_index(drop=True)
            else:
                df = None
        else:
            df = None
        
        # Create figure with two subplots side by side
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Left panel: q vs frame
        frames = df_plot['frame'].values
        
        if window_size > 1:
            q_values = df_plot['q_smooth'].values
            q_label = f'Running avg (window={window_size})'
        else:
            q_values = df_plot['q'].values
            q_label = 'q'
        
        # Plot raw data (light) and running average (bold)
        if raw_values and window_size > 1:
            ax1.plot(df['frame'].values, df['q'].values, alpha=0.3, color='blue', label='Raw', linewidth=0.5)
        ax1.plot(frames, q_values, color='blue', label=q_label, linewidth=2)
        ax1.set_xlabel('Frame Number', fontsize=12)
        ax1.set_ylabel('Fraction of Native Contacts (q)', fontsize=12)
        ax1.set_title('Native Contact Formation Over Time', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Right panel: clusters_filling vs frame
        cluster_numbers = sorted(self.cluster_contacts.keys())
        
        # Create a color map for clusters
        colors = plt.cm.tab10(np.linspace(0, 1, len(cluster_numbers)))
        
        for idx, cluster_num in enumerate(cluster_numbers):
            if window_size > 1:
                cluster_col = f'cluster_{cluster_num}_smooth'
                if cluster_col in df_plot.columns:
                    cluster_values = df_plot[cluster_col].values
                else:
                    continue
            else:
                # Extract cluster filling values from raw data
                cluster_values = []
                for clusters_dict in df_plot['clusters_filling']:
                    clusters_dict = self._parse_clusters_filling(clusters_dict)
                    cluster_values.append(clusters_dict.get(cluster_num, 0.0))
                cluster_values = np.array(cluster_values)
            
            # Use black for cluster_0, otherwise use color from colormap
            plot_color = 'black' if cluster_num == 0 else colors[idx]
            
            # Plot raw data (light) if requested and using smoothed data
            if raw_values and window_size > 1 and df is not None:
                df_raw = df.sort_values('frame').reset_index(drop=True)
                cluster_raw_values = [
                    self._parse_clusters_filling(clusters_dict).get(cluster_num, 0.0)
                    for clusters_dict in df_raw['clusters_filling']
                ]
                ax2.plot(df_raw['frame'].values, cluster_raw_values, alpha=0.3, color=plot_color, linewidth=0.5)
            
            label_suffix = ' (avg)' if window_size > 1 else ''
            ax2.plot(frames, cluster_values, color=plot_color, 
                    label=f'Cluster {cluster_num}{label_suffix}', linewidth=1)
        
        ax2.set_xlabel('Frame Number', fontsize=12)
        ax2.set_ylabel('Cluster Filling Fraction', fontsize=12)
        ax2.set_title('Cluster Formation Over Time', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='best', fontsize=9)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        os.makedirs('results/plots', exist_ok=True)
        filename = f'results/plots/{base_name}_trajectory_analysis.{output_format}'
        plt.savefig(filename, format=output_format, dpi=300, bbox_inches='tight')
        plt.close()
    
    def summarize_trajectory(self,
                            window_size: int = 10000,
                            cutoff: Optional[float] = None) -> None:
        """
        Summarize trajectory by calculating mean cluster filling fractions in windows.
        
        Groups frames into windows and calculates the mean fraction of contacts
        established for each cluster in each window.
        Always saves results to {base_name}_summary.csv or {base_name}_summary_binary_{cutoff}.csv
        in the same directory as the trajectory file.
        
        Args:
            window_size (int): Number of frames per window (default: 10000)
            cutoff (float, optional): If provided, convert probabilities to binary (0 or 1).
                Values >= cutoff become 1, values < cutoff become 0. If None, returns float probabilities.
        """
        # Determine output path
        base_name, output_dir = self._get_base_paths()
        if cutoff is not None:
            output_csv_path = os.path.join(output_dir, f"{base_name}_summary_binary_{cutoff}.csv")
        else:
            output_csv_path = os.path.join(output_dir, f"{base_name}_summary.csv")
        
        # Call Rust function
        folding_analysis_rs.summarize_trajectory(
            trajectory_file=self.trajectory_file_path,
            window_size=window_size,
            cutoff=cutoff,
            output_csv=output_csv_path
        )
    
    def plot_summary(self,
                    cutoff: Optional[float] = None,
                    figsize: Tuple[int, int] = (14, 8),
                    output_format: str = "png") -> None:
        """
        Plot each cluster column from the summary CSV file as a time series.
        
        Loads summary data from {base_name}_summary.csv or {base_name}_summary_binary_{cutoff}.csv
        based on the cutoff parameter.
        
        Args:
            cutoff (float, optional): Cutoff value used for binary summary. If provided, loads
                {base_name}_summary_binary_{cutoff}.csv. If None, loads {base_name}_summary.csv.
            figsize (Tuple[int, int]): Figure size (default: (14, 8))
            output_format (str): Output format - 'png' or 'svg' (default: 'png')
        """
        # Determine input CSV path
        base_name, output_dir = self._get_base_paths()
        if cutoff is not None:
            summary_csv_path = os.path.join(output_dir, f"{base_name}_summary_binary_{cutoff}.csv")
        else:
            summary_csv_path = os.path.join(output_dir, f"{base_name}_summary.csv")
        
        if not os.path.exists(summary_csv_path):
            raise FileNotFoundError(f"Summary CSV file not found: {summary_csv_path}")
        
        # Load summary DataFrame
        summary_df = pd.read_csv(summary_csv_path, index_col=0)
        
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
        
        # Save plot
        os.makedirs('results/plots', exist_ok=True)
        base_name, _ = self._get_base_paths()
        if cutoff is not None:
            filename = f'results/plots/{base_name}_summary_binary_{cutoff}.{output_format}'
        else:
            filename = f'results/plots/{base_name}_summary.{output_format}'
        plt.savefig(filename, format=output_format, dpi=300, bbox_inches='tight')
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
        # Determine summary CSV path
        if isinstance(summary_df_or_path, pd.DataFrame):
            # If DataFrame provided, use auto-load from expected path
            summary_csv_path = None
        else:
            summary_csv_path = summary_df_or_path
            if not os.path.exists(summary_csv_path):
                raise FileNotFoundError(f"Summary CSV file not found: {summary_csv_path}")
        
        # Determine output path
        if output_path is None:
            base_name, output_dir = self._get_base_paths()
            output_path = os.path.join(output_dir, f"{base_name}_class.txt")
        
        # Call Rust function
        formation_order = folding_analysis_rs.classify(
            trajectory_file=self.trajectory_file_path,
            summary_csv=summary_csv_path,
            output_path=output_path
        )
        
        return formation_order
        
    def animate(self,
                csv_path: Optional[str] = None,
                interval: int = 50,
                figsize: Tuple[int, int] = (12, 10),
                save_animation: bool = True,
                output_format: str = "gif",
                fps: int = 20):
        """
        Create an animated contact map showing which contacts are formed in each frame.
        
        Args:
            csv_path (str, optional): Path to parsed CSV file from read_trajectory().
                If None, auto-generates from trajectory file path.
            interval (int): Delay between frames in milliseconds (default: 50)
            figsize (Tuple[int, int]): Figure size (default: (12, 10))
            save_animation (bool): Whether to save the animation (default: True)
            output_format (str): Output format - 'gif' or 'mp4' (default: 'gif')
            fps (int): Frames per second for saved animation (default: 20, only used for mp4)
        """
        # Determine CSV path
        if csv_path is None:
            base_name, output_dir = self._get_base_paths()
            csv_path = os.path.join(output_dir, f"{base_name}_parsed.csv")
        
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Trajectory CSV file not found: {csv_path}")
        
        # Load DataFrame from CSV
        df = pd.read_csv(csv_path)
        
        if df.empty:
            raise ValueError("DataFrame is empty. Cannot create animation.")
        
        # Convert string representations back to lists if needed
        if 'contact_list' in df.columns:
            df['contact_list'] = df['contact_list'].apply(
                lambda x: ast.literal_eval(x) if isinstance(x, str) else x
            )
        
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
        def normalize_contact(contact: Union[tuple, list]) -> tuple:
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
            base_name, _ = self._get_base_paths()
            
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
        plt.close()
        
        return anim