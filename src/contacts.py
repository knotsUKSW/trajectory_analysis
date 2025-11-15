import pandas as pd
import numpy as np
from typing import Optional, List, Tuple, Dict
from collections import defaultdict
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os


class ContactMap:
    """
    Class for loading and managing native contact maps.
    """
    
    def __init__(self, input_file: str):
        """
        Initialize ContactMap by loading contacts from a file.
        
        Args:
            input_file (str): Path to contacts file (whitespace-separated table with columns: i, j, r6, r12)
        """
        self.input_file = input_file
        self.df = self._load_contacts()
    
    def _load_contacts(self) -> pd.DataFrame:
        """
        Load contacts from file and calculate native distances.
        
        Returns:
            pd.DataFrame: DataFrame with columns 'i', 'j', 'r6', 'r12', 'r'
        """
        # Read whitespace-separated file
        df = pd.read_csv(self.input_file, sep=r'\s+', header=None, names=['i', 'j', 'r6', 'r12'])
        
        # Calculate native distance 'r'
        # Formula: r = 10 * (6/5 * (r12/r6))^(1/2)
        df['r'] = 10 * (6/5 * (df['r12'] / df['r6'])) ** (1/2)
        
        return df
    
    @property
    def contacts(self) -> pd.DataFrame:
        """
        Get the contacts DataFrame.
        
        Returns:
            pd.DataFrame: Contacts DataFrame
        """
        return self.df
    
    def __len__(self) -> int:
        """Return the number of contacts."""
        return len(self.df)
    
    def __repr__(self) -> str:
        """String representation of ContactMap."""
        has_clusters = 'cluster' in self.df.columns
        cluster_info = f", clusters={self.df['cluster'].nunique()}" if has_clusters else ""
        return f"ContactMap(file='{self.input_file}', contacts={len(self.df)}{cluster_info})"
    
    def cluster(self, 
                cluster_number: int = 10, 
                cutoff_size: int = 5) -> None:
        """
        Cluster contacts using adjacency-based clustering and add 'cluster' column to dataframe.
        
        Contacts are clustered based on adjacency rules (Manhattan distance 1 or diagonal adjacency).
        Contacts not assigned to any cluster are assigned cluster number 0.
        
        Args:
            cluster_number (int): Desired number of clusters (default: 10)
            cutoff_size (int): Minimum cluster size to retain (default: 5)
        """
        # Extract contact pairs from dataframe
        contacts = [(int(row['i']), int(row['j'])) for _, row in self.df.iterrows()]
        
        # Step 1: Create raw clusters using adjacency rules
        raw_clusters = self._create_raw_clusters(contacts)
        
        # Step 2: Filter clusters by cutoff size
        filtered_clusters = [cluster for cluster in raw_clusters if len(cluster) >= cutoff_size]
        
        if not filtered_clusters:
            # If no clusters meet cutoff, assign all contacts to cluster 0
            self.df['cluster'] = 0
            return
        
        # Step 3: Split clusters to obtain desired number
        final_clusters = self._split_clusters_to_number(filtered_clusters, cluster_number)
        
        # Step 4: Create mapping from contact to cluster number
        # Initialize all contacts with cluster 0 (not assigned)
        contact_to_cluster = {contact: 0 for contact in contacts}
        
        # Assign cluster numbers (1, 2, 3, ...) to clustered contacts
        clustered_contacts = set()
        for cluster_idx, cluster in enumerate(final_clusters, start=1):
            for contact in cluster:
                contact_to_cluster[contact] = cluster_idx
                clustered_contacts.add(contact)
        
        # Step 5: Map cluster assignments to dataframe
        # Create a function to get cluster number for each row
        def get_cluster(row):
            contact = (int(row['i']), int(row['j']))
            return contact_to_cluster.get(contact, 0)
        
        # Add cluster column
        self.df['cluster'] = self.df.apply(get_cluster, axis=1)
    
    def _are_adjacent_contacts(self, contact1: Tuple[int, int], contact2: Tuple[int, int]) -> bool:
        """
        Check if two contacts are adjacent according to the clustering rules.
        
        Two contacts (i1, j1) and (i2, j2) are adjacent if:
        1. Manhattan distance = 1: |i1-i2| + |j1-j2| = 1
        2. Diagonal adjacency: |i1-i2| = 1 and |j1-j2| = 1
        
        Args:
            contact1 (Tuple[int, int]): First contact pair
            contact2 (Tuple[int, int]): Second contact pair
            
        Returns:
            bool: True if contacts are adjacent
        """
        i1, j1 = contact1
        i2, j2 = contact2
        
        # Check Manhattan distance = 1
        manhattan_dist = abs(i1 - i2) + abs(j1 - j2)
        if manhattan_dist == 1:
            return True
        
        # Check diagonal adjacency: |i1-i2| = 1 and |j1-j2| = 1
        diagonal_adjacent = (abs(i1 - i2) == 1) and (abs(j1 - j2) == 1)
        
        return diagonal_adjacent
    
    def _create_raw_clusters(self, contacts: List[Tuple[int, int]]) -> List[List[Tuple[int, int]]]:
        """
        Create raw clusters using adjacency-based clustering.
        
        Args:
            contacts (List[Tuple[int, int]]): List of contact pairs
            
        Returns:
            List[List[Tuple[int, int]]]: List of raw clusters
        """
        # Create adjacency graph
        adjacency_graph = defaultdict(set)
        
        # Build adjacency relationships
        for i, contact1 in enumerate(contacts):
            for j, contact2 in enumerate(contacts):
                if i >= j:  # Avoid duplicate pairs
                    continue
                    
                if self._are_adjacent_contacts(contact1, contact2):
                    adjacency_graph[contact1].add(contact2)
                    adjacency_graph[contact2].add(contact1)
        
        # Find connected components (clusters) using DFS
        visited = set()
        clusters = []
        
        for contact in contacts:
            if contact not in visited:
                cluster = []
                stack = [contact]
                
                while stack:
                    current = stack.pop()
                    if current not in visited:
                        visited.add(current)
                        cluster.append(current)
                        
                        # Add all adjacent contacts to stack
                        for neighbor in adjacency_graph[current]:
                            if neighbor not in visited:
                                stack.append(neighbor)
                
                clusters.append(cluster)
        
        return clusters
    
    def _calculate_cluster_centroid(self, cluster: List[Tuple[int, int]]) -> Tuple[float, float]:
        """
        Calculate the centroid of a cluster.
        
        Args:
            cluster (List[Tuple[int, int]]): List of contact pairs in the cluster
            
        Returns:
            Tuple[float, float]: Centroid coordinates
        """
        if not cluster:
            return (0.0, 0.0)
        
        all_i = [contact[0] for contact in cluster]
        all_j = [contact[1] for contact in cluster]
        
        centroid_i = np.mean(all_i)
        centroid_j = np.mean(all_j)
        
        return (centroid_i, centroid_j)
    
    def _split_cluster(self, cluster: List[Tuple[int, int]]) -> List[List[Tuple[int, int]]]:
        """
        Split a cluster into two smaller clusters using k-means.
        
        Args:
            cluster (List[Tuple[int, int]]): Cluster to split
            
        Returns:
            List[List[Tuple[int, int]]]: List of split clusters
        """
        if len(cluster) < 4:
            return [cluster]
        
        # Convert contacts to feature vectors for k-means
        features = np.array(cluster)
        
        # Use k-means to split into 2 clusters
        kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
        labels = kmeans.fit_predict(features)
        
        # Split based on labels
        cluster1 = [cluster[i] for i in range(len(cluster)) if labels[i] == 0]
        cluster2 = [cluster[i] for i in range(len(cluster)) if labels[i] == 1]
        
        # Ensure both clusters have at least 2 contacts
        if len(cluster1) < 2 or len(cluster2) < 2:
            return [cluster]
        
        return [cluster1, cluster2]
    
    def _split_clusters_to_number(self, 
                                  clusters: List[List[Tuple[int, int]]], 
                                  target_number: int) -> List[List[Tuple[int, int]]]:
        """
        Split clusters to obtain the desired number of clusters.
        
        Args:
            clusters (List[List[Tuple[int, int]]]): List of clusters to split
            target_number (int): Desired number of clusters
            
        Returns:
            List[List[Tuple[int, int]]]: List of final clusters
        """
        if len(clusters) >= target_number:
            # If we have enough clusters, return the largest ones
            clusters.sort(key=len, reverse=True)
            return clusters[:target_number]
        
        # If we need more clusters, split the largest ones
        result_clusters = clusters.copy()
        
        while len(result_clusters) < target_number:
            # Find the largest cluster to split
            largest_cluster_idx = max(range(len(result_clusters)), 
                                     key=lambda i: len(result_clusters[i]))
            largest_cluster = result_clusters[largest_cluster_idx]
            
            if len(largest_cluster) < 4:  # Can't split clusters smaller than 4
                break
            
            # Split the largest cluster
            split_clusters = self._split_cluster(largest_cluster)
            
            # Replace the largest cluster with its splits
            result_clusters.pop(largest_cluster_idx)
            result_clusters.extend(split_clusters)
        
        return result_clusters
    
    def draw(self,
             figsize: Tuple[int, int] = (12, 10),
             title: str = "Contact Map with Clusters",
             output_format: str = "png",
             save_plot: bool = True,
             show_plot: bool = True) -> None:
        """
        Draw a contact map showing contacts as rectangles below the diagonal, 
        with different colors for each cluster.
        
        Requires that clustering has been performed (cluster column must exist).
        Contacts with cluster=0 are shown in black as unassigned/removed contacts.
        
        Args:
            figsize (Tuple[int, int]): Figure size for the plot (default: (12, 10))
            title (str): Title for the plot (default: "Contact Map with Clusters")
            output_format (str): Output format ('png' or 'svg', default: 'png')
            save_plot (bool): Whether to save the plot to results/plots (default: True)
            show_plot (bool): Whether to display the plot (default: True)
        """
        # Check if clustering has been performed
        if 'cluster' not in self.df.columns:
            raise ValueError("Clustering must be performed before drawing. Call cluster() method first.")
        
        # Get all contacts
        all_contacts = [(int(row['i']), int(row['j'])) for _, row in self.df.iterrows()]
        
        if not all_contacts:
            raise ValueError("No contacts to plot")
        
        # Find the maximum residue number
        max_residue = max(max(contact) for contact in all_contacts)
        
        # Group contacts by cluster
        clusters_dict = {}  # Maps cluster number to list of contacts
        cluster_numbers = []  # Keep track of cluster numbers in order
        removed_contacts = []
        
        for cluster_num in sorted(self.df['cluster'].unique()):
            cluster_df = self.df[self.df['cluster'] == cluster_num]
            cluster_contacts = [(int(row['i']), int(row['j'])) for _, row in cluster_df.iterrows()]
            
            if cluster_num == 0:
                removed_contacts = cluster_contacts
            else:
                # Store cluster number and contacts
                clusters_dict[cluster_num] = cluster_contacts
                cluster_numbers.append(cluster_num)
        
        if not clusters_dict and not removed_contacts:
            raise ValueError("No clusters or contacts to plot")
        
        # Create figure and axis
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set up colors for clusters
        num_clusters = len(clusters_dict)
        if num_clusters > 0:
            colors = plt.cm.tab20(np.linspace(0, 1, num_clusters))
        else:
            colors = []
        
        cluster_colors = {}
        
        # Plot contacts for each cluster (only below diagonal)
        for i, cluster_num in enumerate(cluster_numbers):
            cluster_contacts = clusters_dict[cluster_num]
            color = colors[i]
            cluster_colors[f'Cluster {cluster_num}'] = color
            
            # Plot contacts as rectangles (only below diagonal)
            for contact in cluster_contacts:
                i_res, j_res = contact
                
                # Only plot contacts below the diagonal (i < j)
                if i_res < j_res:
                    # Swap axes: plot (j_res, i_res) so that y < x
                    rect = patches.Rectangle(
                        (j_res - 0.5, i_res - 0.5),  # Bottom-left corner (swapped)
                        1, 1,  # Width and height
                        linewidth=0.5,
                        edgecolor='black',
                        facecolor=color,
                        alpha=0.7
                    )
                    ax.add_patch(rect)
        
        # Plot removed contacts (cluster=0) in black (only below diagonal)
        if removed_contacts:
            cluster_colors['Unassigned (cluster 0)'] = 'black'
            for contact in removed_contacts:
                i_res, j_res = contact
                
                # Only plot contacts below the diagonal (i < j)
                if i_res < j_res:
                    # Swap axes: plot (j_res, i_res) so that y < x
                    rect = patches.Rectangle(
                        (j_res - 0.5, i_res - 0.5),  # Bottom-left corner (swapped)
                        1, 1,
                        linewidth=0.5,
                        edgecolor='black',
                        facecolor='black',
                        alpha=0.7
                    )
                    ax.add_patch(rect)
        
        # Set axis properties
        ax.set_xlim(0.5, max_residue + 0.5)
        ax.set_ylim(0.5, max_residue + 0.5)
        ax.set_xlabel('Residue Index (j)', fontsize=12)
        ax.set_ylabel('Residue Index (i)', fontsize=12)
        ax.set_title(title, fontsize=14)
        
        # Set aspect ratio to be equal
        ax.set_aspect('equal')
        
        # Add diagonal line
        ax.plot([0.5, max_residue + 0.5], [0.5, max_residue + 0.5], 
                'k--', alpha=0.5, linewidth=1, label='Diagonal')
        
        # Create legend
        legend_elements = []
        for cluster_name, color in cluster_colors.items():
            legend_elements.append(patches.Patch(color=color, label=cluster_name))
        
        # Add diagonal to legend
        legend_elements.append(plt.Line2D([0], [0], color='black', linestyle='--', 
                                         alpha=0.5, label='Diagonal'))
        
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Add statistics text
        clustered_contacts = sum(len(cluster) for cluster in clusters_dict.values())
        removed_count = len(removed_contacts)
        total_contacts = clustered_contacts + removed_count
        
        contacts_below_diagonal = sum(1 for contact in all_contacts if contact[0] < contact[1])
        
        stats_text = f'Total contacts: {total_contacts}\nClustered: {clustered_contacts}\nUnassigned: {removed_count}\nBelow diagonal: {contacts_below_diagonal}\nClusters: {num_clusters}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # Save plot if requested
        if save_plot:
            os.makedirs('results/plots', exist_ok=True)
            base_name = os.path.splitext(os.path.basename(self.input_file))[0]
            filename = f'results/plots/{base_name}_contact_map.{output_format}'
            plt.savefig(filename, format=output_format, dpi=300, bbox_inches='tight')
        
        if show_plot:
            plt.show()
        else:
            plt.close()

