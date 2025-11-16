use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::contacts::Contact;
use crate::structure::{Coordinate, FrameData};
use indicatif::{ProgressBar, ProgressStyle};

/// Trajectory data: maps frame number to frame data
pub type TrajectoryData = HashMap<i32, FrameData>;

/// Result of trajectory analysis for a single frame
#[derive(Debug, Clone)]
pub struct FrameResult {
    pub frame: i32,
    pub contacts: usize,
    pub q: f64,
    pub contact_list: Vec<(i32, i32)>,
    pub clusters_filling: HashMap<i32, f64>,
}

/// Summary of trajectory analysis for a window of frames
#[derive(Debug, Clone)]
pub struct WindowSummary {
    pub frame: i32,  // Starting frame of the window
    pub cluster_means: HashMap<i32, f64>,  // Mean fraction for each cluster
}

/// Smoothed result for a single frame (running average)
#[derive(Debug, Clone)]
pub struct SmoothedResult {
    pub frame: i32,
    pub q_smooth: f64,  // Smoothed q value
    pub cluster_smooth: HashMap<i32, f64>,  // Smoothed fraction for each cluster
}

/// Trait for reading and processing trajectory files
pub trait Trajectory {
    /// Read PDB file and extract CA atom coordinates for each frame
    /// 
    /// # Arguments
    /// * `max_frames` - Maximum number of frames to read (None for all frames)
    /// 
    /// # Returns
    /// HashMap mapping frame number to HashMap mapping residue number to CA coordinates
    fn read_pdb(&self, max_frames: Option<usize>) -> Result<TrajectoryData, String>;
    
    /// Read trajectory and calculate native contact formation
    /// 
    /// For each frame, calculates which contacts are formed based on distance cutoff.
    /// 
    /// # Arguments
    /// * `contacts` - Vector of Contact structs with native contact information
    /// * `cutoff_distance` - Multiplier for native distance cutoff (default: 1.2)
    /// * `max_frames` - Maximum number of frames to process (None for all frames)
    /// * `output_csv_path` - Optional path to save CSV output. If None, auto-generates from input path
    /// 
    /// # Returns
    /// Vector of FrameResult for each frame
    fn read_trajectory(
        &self,
        contacts: &[Contact],
        cutoff_distance: f64,
        max_frames: Option<usize>,
        output_csv_path: Option<&Path>,
    ) -> Result<Vec<FrameResult>, String>;
    
    /// Summarize trajectory by calculating mean cluster filling fractions in windows.
    /// 
    /// Groups frames into windows and calculates the mean fraction of contacts
    /// established for each cluster in each window.
    /// 
    /// # Arguments
    /// * `results` - Vector of FrameResult from read_trajectory. If None, will attempt to load
    ///   from auto-generated CSV path based on trajectory file name.
    /// * `window_size` - Number of frames per window (default: 10000)
    /// * `cutoff` - Optional cutoff for binary conversion. If Some(value), convert probabilities
    ///   to binary (0 or 1). Values >= cutoff become 1, values < cutoff become 0.
    /// * `output_csv_path` - Optional path to save summary CSV. If None, auto-generates from input path
    /// 
    /// # Returns
    /// Vector of WindowSummary for each window
    fn summarize_trajectory(
        &self,
        results: Option<&[FrameResult]>,
        window_size: usize,
        cutoff: Option<f64>,
        output_csv_path: Option<&Path>,
    ) -> Result<Vec<WindowSummary>, String>;
    
    /// Smooth trajectory data by calculating running averages.
    /// 
    /// Calculates running average of q (fraction of all contacts) and per-cluster
    /// fraction of contacts formed using a centered window.
    /// 
    /// # Arguments
    /// * `results` - Vector of FrameResult from read_trajectory. If None, will attempt to load
    ///   from auto-generated CSV path based on trajectory file name.
    /// * `window_size` - Window size for running average (default: 100)
    /// * `output_csv_path` - Optional path to save smoothed CSV. If None, auto-generates from input path
    /// 
    /// # Returns
    /// Vector of SmoothedResult for each frame
    fn smooth(
        &self,
        results: Option<&[FrameResult]>,
        window_size: usize,
        output_csv_path: Option<&Path>,
    ) -> Result<Vec<SmoothedResult>, String>;
    
    /// Classify trajectory by determining cluster formation order.
    /// 
    /// Works backwards from the first time all clusters are formed to determine
    /// the order in which clusters were last formed (formation order).
    /// 
    /// # Arguments
    /// * `summary_csv_path` - Path to summary CSV file (binary representation from summarize_trajectory with cutoff)
    /// * `output_path` - Optional path to save the formation order. If None, auto-generates from input path
    /// 
    /// # Returns
    /// Vector of cluster numbers in formation order
    fn classify(
        &self,
        summary_csv_path: Option<&Path>,
        output_path: Option<&Path>,
    ) -> Result<Vec<i32>, String>;
}

/// Implementation of Trajectory trait for PDB files
pub struct PdbTrajectory {
    file_path: String,
}

impl PdbTrajectory {
    pub fn new(file_path: impl AsRef<Path>) -> Self {
        Self {
            file_path: file_path.as_ref().to_string_lossy().to_string(),
        }
    }
}

impl Trajectory for PdbTrajectory {
    fn read_pdb(&self, max_frames: Option<usize>) -> Result<TrajectoryData, String> {
        let file = File::open(&self.file_path)
            .map_err(|e| format!("Failed to open PDB file {}: {}", self.file_path, e))?;
        
        let reader = BufReader::new(file);
        let mut frames_data = TrajectoryData::new();
        
        let mut current_model: Option<i32> = None;
        let mut current_residue_coords = FrameData::new();
        let mut model_found = false;
        let mut model_saved = false;
        
        for line_result in reader.lines() {
            let line = line_result.map_err(|e| format!("Error reading line: {}", e))?;
            
            if line.starts_with("MODEL") {
                // Start of new model
                model_found = true;
                
                // Save previous model if it hasn't been saved yet
                if let Some(model_num) = current_model.take() {
                    if !model_saved {
                        frames_data.insert(model_num, current_residue_coords.clone());
                    }
                }
                
                // Check if we've reached max_frames limit
                if let Some(max) = max_frames {
                    if frames_data.len() >= max {
                        break;
                    }
                }
                
                // Parse model number
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() < 2 {
                    return Err(format!("Invalid MODEL line: {}", line));
                }
                
                current_model = Some(
                    parts[1]
                        .parse::<i32>()
                        .map_err(|e| format!("Failed to parse model number: {}", e))?,
                );
                current_residue_coords = FrameData::new();
                model_saved = false;
                
            } else if line.starts_with("ATOM") && line.contains(" CA ") {
                // Extract CA atom coordinates
                // PDB format: columns 22-26 = residue number, 30-38 = x, 38-46 = y, 46-54 = z
                if line.len() < 54 {
                    continue; // Skip malformed lines
                }
                
                // Extract residue number (columns 22-26, 1-indexed)
                let residue_str = line.get(22..26).unwrap_or("").trim();
                let residue_num = match residue_str.parse::<i32>() {
                    Ok(n) => n,
                    Err(_) => continue, // Skip if can't parse residue number
                };
                
                // Extract coordinates
                let x_str = line.get(30..38).unwrap_or("").trim();
                let y_str = line.get(38..46).unwrap_or("").trim();
                let z_str = line.get(46..54).unwrap_or("").trim();
                
                let x = match x_str.parse::<f64>() {
                    Ok(v) => v,
                    Err(_) => continue,
                };
                let y = match y_str.parse::<f64>() {
                    Ok(v) => v,
                    Err(_) => continue,
                };
                let z = match z_str.parse::<f64>() {
                    Ok(v) => v,
                    Err(_) => continue,
                };
                
                current_residue_coords.insert(residue_num, Coordinate::new(x, y, z));
                
            } else if line.starts_with("ENDMDL") {
                // End of model - save the data
                if let Some(model_num) = current_model {
                    if !model_saved {
                        frames_data.insert(model_num, current_residue_coords.clone());
                        model_saved = true;
                        current_residue_coords = FrameData::new();
                    }
                }
                
                // Check if we've reached max_frames limit
                if let Some(max) = max_frames {
                    if frames_data.len() >= max {
                        break;
                    }
                }
            }
        }
        
        // Handle last model if file doesn't end with ENDMDL
        if max_frames.is_none() || frames_data.len() < max_frames.unwrap_or(0) {
            if let Some(model_num) = current_model {
                if !model_saved && !current_residue_coords.is_empty() {
                    frames_data.insert(model_num, current_residue_coords.clone());
                }
            }
            
            // Handle single-model files without MODEL/ENDMDL markers
            if !model_found && !current_residue_coords.is_empty() {
                frames_data.insert(1, current_residue_coords);
            }
        }
        
        Ok(frames_data)
    }
    
    fn read_trajectory(
        &self,
        contacts: &[Contact],
        cutoff_distance: f64,
        max_frames: Option<usize>,
        output_csv_path: Option<&Path>,
    ) -> Result<Vec<FrameResult>, String> {
        // Read PDB file
        let frames_data = self.read_pdb(max_frames)?;
        
        if frames_data.is_empty() {
            return Ok(Vec::new());
        }
        
        // Pre-compute cluster information
        let mut cluster_contacts: HashMap<i32, Vec<usize>> = HashMap::new();
        for (idx, contact) in contacts.iter().enumerate() {
            cluster_contacts
                .entry(contact.cluster)
                .or_insert_with(Vec::new)
                .push(idx);
        }
        
        let cluster_sizes: HashMap<i32, usize> = cluster_contacts
            .iter()
            .map(|(cluster, indices)| (*cluster, indices.len()))
            .collect();
        
        let total_contacts = contacts.len();
        
        // Process each frame with progress bar
        let mut results = Vec::new();
        let mut sorted_frames: Vec<i32> = frames_data.keys().copied().collect();
        sorted_frames.sort();
        
        let total_frames = sorted_frames.len();
        let pb = ProgressBar::new(total_frames as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} frames ({percent}%) | ETA: {eta}")
                .unwrap()
                .progress_chars("#>-")
        );
        pb.set_message("Processing trajectory frames");
        
        for frame_num in sorted_frames {
            let residue_coords = &frames_data[&frame_num];
            
            let mut existing_contacts = Vec::new();
            let mut cluster_counts: HashMap<i32, usize> = HashMap::new();
            
            // Check each native contact
            for contact in contacts {
                // Check if both residues exist in the structure
                if let (Some(coord_i), Some(coord_j)) = 
                    (residue_coords.get(&contact.i), residue_coords.get(&contact.j)) {
                    // Calculate distance between CA atoms
                    let distance = coord_i.distance_to(coord_j);
                    
                    // Check if contact exists (distance < r_native * cutoff_distance)
                    if distance < contact.r * cutoff_distance {
                        existing_contacts.push((contact.i, contact.j));
                        *cluster_counts.entry(contact.cluster).or_insert(0) += 1;
                    }
                }
            }
            
            // Calculate q (fraction of native contacts existing)
            let q = if total_contacts > 0 {
                existing_contacts.len() as f64 / total_contacts as f64
            } else {
                0.0
            };
            
            // Calculate clusters_filling (fraction of contacts in each cluster)
            let mut clusters_filling = HashMap::new();
            for (cluster_num, count) in cluster_counts.iter() {
                let cluster_size = cluster_sizes.get(cluster_num).copied().unwrap_or(0);
                if cluster_size > 0 {
                    clusters_filling.insert(*cluster_num, *count as f64 / cluster_size as f64);
                } else {
                    clusters_filling.insert(*cluster_num, 0.0);
                }
            }
            
            // Ensure all clusters are represented (even if no contacts formed)
            for cluster_num in cluster_sizes.keys() {
                clusters_filling.entry(*cluster_num).or_insert(0.0);
            }
            
            results.push(FrameResult {
                frame: frame_num,
                contacts: existing_contacts.len(),
                q,
                contact_list: existing_contacts,
                clusters_filling,
            });
            
            pb.inc(1);
        }
        
        pb.finish_with_message("Processing complete");
        
        // Save to CSV if requested
        if let Some(output_path) = output_csv_path {
            save_results_to_csv(&results, output_path)?;
        }
        
        Ok(results)
    }
    
    fn summarize_trajectory(
        &self,
        results: Option<&[FrameResult]>,
        window_size: usize,
        cutoff: Option<f64>,
        output_csv_path: Option<&Path>,
    ) -> Result<Vec<WindowSummary>, String> {
        // Load results if not provided
        let frame_results = if let Some(res) = results {
            res.to_vec()
        } else {
            // Try to load from auto-generated CSV path
            let base_name = Path::new(&self.file_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("trajectory");
            let traj_dir = Path::new(&self.file_path)
                .parent()
                .unwrap_or(Path::new("."));
            let csv_path = traj_dir.join(format!("{}_parsed.csv", base_name));
            
            if !csv_path.exists() {
                return Err(format!(
                    "Trajectory data not found. Please run read_trajectory() first or provide results. \
                     Expected file: {}",
                    csv_path.display()
                ));
            }
            
            load_results_from_csv(&csv_path)?
        };
        
        if frame_results.is_empty() {
            return Err("Frame results are empty. Cannot create summary.".to_string());
        }
        
        // Sort by frame number
        let mut sorted_results = frame_results;
        sorted_results.sort_by_key(|r| r.frame);
        
        // Get all unique cluster numbers
        let mut cluster_numbers: Vec<i32> = sorted_results
            .iter()
            .flat_map(|r| r.clusters_filling.keys().copied())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        cluster_numbers.sort();
        
        // Group frames into windows
        let total_frames = sorted_results.len();
        let num_windows = (total_frames + window_size - 1) / window_size;  // Ceiling division
        
        let mut window_summaries = Vec::new();
        
        for window_idx in 0..num_windows {
            let start_idx = window_idx * window_size;
            let end_idx = (start_idx + window_size).min(total_frames);
            
            let window_results = &sorted_results[start_idx..end_idx];
            
            // Get the frame number for this window (use the first frame in the window)
            let frame_number = window_results[0].frame;
            
            // Calculate mean fraction for each cluster in this window
            let mut cluster_means = HashMap::new();
            
            for cluster_num in &cluster_numbers {
                let mut cluster_values = Vec::new();
                
                for result in window_results {
                    let fraction = result.clusters_filling.get(cluster_num).copied().unwrap_or(0.0);
                    cluster_values.push(fraction);
                }
                
                // Calculate mean for this cluster in this window
                let mean_fraction = if !cluster_values.is_empty() {
                    cluster_values.iter().sum::<f64>() / cluster_values.len() as f64
                } else {
                    0.0
                };
                
                // Apply binary conversion if cutoff is provided
                let final_value = if let Some(cutoff_val) = cutoff {
                    if mean_fraction >= cutoff_val {
                        1.0
                    } else {
                        0.0
                    }
                } else {
                    mean_fraction
                };
                
                cluster_means.insert(*cluster_num, final_value);
            }
            
            window_summaries.push(WindowSummary {
                frame: frame_number,
                cluster_means,
            });
        }
        
        // Auto-generate output path if not provided
        let output_path_buf: Option<std::path::PathBuf> = if let Some(path) = output_csv_path {
            Some(path.to_path_buf())
        } else {
            // Auto-generate from trajectory file path
            let base_name = Path::new(&self.file_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("trajectory");
            let traj_dir = Path::new(&self.file_path)
                .parent()
                .unwrap_or(Path::new("."));
            
            // Use _summary_binary.csv if cutoff is provided, otherwise _summary.csv
            let suffix = if cutoff.is_some() {
                "_summary_binary.csv"
            } else {
                "_summary.csv"
            };
            Some(traj_dir.join(format!("{}{}", base_name, suffix)))
        };
        
        // Save to CSV
        if let Some(output_path) = output_path_buf.as_deref() {
            save_summary_to_csv(&window_summaries, &cluster_numbers, output_path)?;
        }
        
        Ok(window_summaries)
    }
    
    fn smooth(
        &self,
        results: Option<&[FrameResult]>,
        window_size: usize,
        output_csv_path: Option<&Path>,
    ) -> Result<Vec<SmoothedResult>, String> {
        // Load results if not provided
        let frame_results = if let Some(res) = results {
            res.to_vec()
        } else {
            // Try to load from auto-generated CSV path
            let base_name = Path::new(&self.file_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("trajectory");
            let traj_dir = Path::new(&self.file_path)
                .parent()
                .unwrap_or(Path::new("."));
            let csv_path = traj_dir.join(format!("{}_parsed.csv", base_name));
            
            if !csv_path.exists() {
                return Err(format!(
                    "Trajectory data not found. Please run read_trajectory() first or provide results. \
                     Expected file: {}",
                    csv_path.display()
                ));
            }
            
            load_results_from_csv(&csv_path)?
        };
        
        if frame_results.is_empty() {
            return Err("Frame results are empty. Cannot smooth data.".to_string());
        }
        
        // Sort by frame number
        let mut sorted_results = frame_results;
        sorted_results.sort_by_key(|r| r.frame);
        
        // Get all unique cluster numbers
        let mut cluster_numbers: Vec<i32> = sorted_results
            .iter()
            .flat_map(|r| r.clusters_filling.keys().copied())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        cluster_numbers.sort();
        
        let total_frames = sorted_results.len();
        let half_window = window_size / 2;
        
        // Calculate running averages
        let mut smoothed_results = Vec::new();
        
        for (idx, result) in sorted_results.iter().enumerate() {
            // Calculate running average for q
            let start_idx = idx.saturating_sub(half_window);
            let end_idx = (idx + half_window + 1).min(total_frames);
            
            let q_window: Vec<f64> = sorted_results[start_idx..end_idx]
                .iter()
                .map(|r| r.q)
                .collect();
            let q_smooth = q_window.iter().sum::<f64>() / q_window.len() as f64;
            
            // Calculate running average for each cluster
            let mut cluster_smooth = HashMap::new();
            
            for cluster_num in &cluster_numbers {
                let cluster_window: Vec<f64> = sorted_results[start_idx..end_idx]
                    .iter()
                    .map(|r| r.clusters_filling.get(cluster_num).copied().unwrap_or(0.0))
                    .collect();
                let cluster_mean = cluster_window.iter().sum::<f64>() / cluster_window.len() as f64;
                cluster_smooth.insert(*cluster_num, cluster_mean);
            }
            
            smoothed_results.push(SmoothedResult {
                frame: result.frame,
                q_smooth,
                cluster_smooth,
            });
        }
        
        // Auto-generate output path if not provided
        let output_path_buf: Option<std::path::PathBuf> = if let Some(path) = output_csv_path {
            Some(path.to_path_buf())
        } else {
            // Auto-generate from trajectory file path
            let base_name = Path::new(&self.file_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("trajectory");
            let traj_dir = Path::new(&self.file_path)
                .parent()
                .unwrap_or(Path::new("."));
            Some(traj_dir.join(format!("{}_smoothed.csv", base_name)))
        };
        
        // Save to CSV
        if let Some(output_path) = output_path_buf.as_deref() {
            save_smoothed_to_csv(&smoothed_results, &cluster_numbers, output_path)?;
        }
        
        Ok(smoothed_results)
    }
    
    fn classify(
        &self,
        summary_csv_path: Option<&Path>,
        output_path: Option<&Path>,
    ) -> Result<Vec<i32>, String> {
        // Determine summary CSV path
        let csv_path_buf: std::path::PathBuf = if let Some(path) = summary_csv_path {
            path.to_path_buf()
        } else {
            // Auto-generate from trajectory file path
            let base_name = Path::new(&self.file_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("trajectory");
            let traj_dir = Path::new(&self.file_path)
                .parent()
                .unwrap_or(Path::new("."));
            traj_dir.join(format!("{}_summary_binary.csv", base_name))
        };
        
        if !csv_path_buf.exists() {
            return Err(format!(
                "Summary CSV file not found. Please run summarize_trajectory() with cutoff first. \
                 Expected file: {}",
                csv_path_buf.display()
            ));
        }
        
        // Load summary CSV
        let summary_data = load_summary_csv(&csv_path_buf)?;
        
        if summary_data.is_empty() {
            return Err("Summary data is empty. Cannot determine formation order.".to_string());
        }
        
        // Get all cluster numbers (excluding cluster_0)
        let mut cluster_numbers: Vec<i32> = summary_data[0]
            .cluster_means
            .keys()
            .filter(|&&k| k != 0)
            .copied()
            .collect();
        cluster_numbers.sort();
        
        if cluster_numbers.is_empty() {
            return Err("No clusters found in summary data (excluding cluster_0). Cannot determine formation order.".to_string());
        }
        
        // Sort by frame number
        let mut sorted_data = summary_data;
        sorted_data.sort_by_key(|w| w.frame);
        
        // Step 2: Find first time all clusters (except cluster_0) are formed (all 1)
        let mut start_idx = None;
        let mut max_clusters_formed = 0;
        let mut max_clusters_frame_idx = 0;
        
        for (idx, window) in sorted_data.iter().enumerate() {
            let clusters_formed: usize = cluster_numbers
                .iter()
                .map(|&cn| {
                    if window.cluster_means.get(&cn).copied().unwrap_or(0.0) >= 1.0 {
                        1
                    } else {
                        0
                    }
                })
                .sum();
            
            if clusters_formed > max_clusters_formed {
                max_clusters_formed = clusters_formed;
                max_clusters_frame_idx = idx;
            }
            
            // Check if all clusters are formed
            if clusters_formed == cluster_numbers.len() {
                start_idx = Some(idx);
                break;
            }
        }
        
        let start_idx = start_idx.unwrap_or(max_clusters_frame_idx);
        
        if max_clusters_formed == 0 {
            return Err("No clusters are formed in the trajectory. Cannot determine formation order.".to_string());
        }
        
        // Step 3-5: Go backwards in time and track cluster breaks
        // Track which clusters are currently formed (start with all formed at start_frame)
        let mut currently_formed: std::collections::HashSet<i32> = cluster_numbers
            .iter()
            .filter(|&&cn| {
                sorted_data[start_idx].cluster_means.get(&cn).copied().unwrap_or(0.0) >= 1.0
            })
            .copied()
            .collect();
        
        // List to store order of cluster breaks (going backwards)
        // Track which clusters have already been recorded (each cluster should appear only once)
        let mut recorded_clusters = std::collections::HashSet::new();
        let mut cluster_breaks_order = Vec::new();
        
        // Go backwards from start_idx
        for idx in (0..=start_idx).rev() {
            let window = &sorted_data[idx];
            
            // Check which clusters are formed at this frame
            let mut formed_at_frame = std::collections::HashSet::new();
            for &cluster_num in &cluster_numbers {
                if window.cluster_means.get(&cluster_num).copied().unwrap_or(0.0) >= 1.0 {
                    formed_at_frame.insert(cluster_num);
                }
            }
            
            // Find clusters that were formed before but are now broken (1 -> 0)
            let broken_clusters: Vec<i32> = currently_formed
                .difference(&formed_at_frame)
                .copied()
                .collect::<Vec<_>>();
            
            // Add broken clusters to the list ONLY if they haven't been recorded yet
            // (in sorted order for consistency)
            if !broken_clusters.is_empty() {
                let mut sorted_broken: Vec<i32> = broken_clusters
                    .into_iter()
                    .filter(|&cn| !recorded_clusters.contains(&cn))
                    .collect();
                sorted_broken.sort();
                
                // Mark these clusters as recorded
                for &cn in &sorted_broken {
                    recorded_clusters.insert(cn);
                }
                
                cluster_breaks_order.extend(sorted_broken);
            }
            
            // Update currently_formed to reflect current state
            currently_formed = formed_at_frame;
            
            // If all clusters (except cluster_0) have been recorded, we're done
            if recorded_clusters.len() == cluster_numbers.len() {
                break;
            }
        }
        
        // Step 6: If we reached the beginning and some clusters are still formed,
        // append them to the end (these were formed before the trajectory started)
        // Only add clusters that haven't been recorded yet
        if !currently_formed.is_empty() {
            let mut still_formed: Vec<i32> = currently_formed
                .into_iter()
                .filter(|&cn| !recorded_clusters.contains(&cn))
                .collect();
            still_formed.sort();
            cluster_breaks_order.extend(still_formed);
        }
        
        // Reverse the order to get formation order (last formed first in breaks = first formed last in formation)
        let formation_order: Vec<i32> = cluster_breaks_order.into_iter().rev().collect();
        
        // Auto-generate output path if not provided
        let output_path_buf: Option<std::path::PathBuf> = if let Some(path) = output_path {
            Some(path.to_path_buf())
        } else {
            // Auto-generate from trajectory file path
            let base_name = Path::new(&self.file_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("trajectory");
            let traj_dir = Path::new(&self.file_path)
                .parent()
                .unwrap_or(Path::new("."));
            Some(traj_dir.join(format!("{}_class.txt", base_name)))
        };
        
        // Save to file
        if let Some(output_path) = output_path_buf.as_deref() {
            save_classification_to_file(&formation_order, output_path)?;
        }
        
        Ok(formation_order)
    }
}

/// Save frame results to CSV file
fn save_results_to_csv(results: &[FrameResult], output_path: &Path) -> Result<(), String> {
    let mut writer = csv::Writer::from_path(output_path)
        .map_err(|e| format!("Failed to create CSV file {}: {}", output_path.display(), e))?;
    
    // Write header
    writer.write_record(&["frame", "contacts", "q", "contact_list", "clusters_filling"])
        .map_err(|e| format!("Failed to write CSV header: {}", e))?;
    
    // Write each row
    for result in results {
        // Format contact_list as string: "[(i1, j1), (i2, j2), ...]"
        let contact_list_str = format!(
            "[{}]",
            result.contact_list
                .iter()
                .map(|(i, j)| format!("({}, {})", i, j))
                .collect::<Vec<_>>()
                .join(", ")
        );
        
        // Format clusters_filling as string: "{cluster1: fraction1, cluster2: fraction2, ...}"
        let mut cluster_entries: Vec<(i32, f64)> = result.clusters_filling.iter()
            .map(|(k, v)| (*k, *v))
            .collect();
        cluster_entries.sort_by_key(|(k, _)| *k);
        let clusters_filling_str = format!(
            "{{{}}}",
            cluster_entries
                .iter()
                .map(|(k, v)| format!("{}: {}", k, v))
                .collect::<Vec<_>>()
                .join(", ")
        );
        
        writer.write_record(&[
            result.frame.to_string(),
            result.contacts.to_string(),
            result.q.to_string(),
            contact_list_str,
            clusters_filling_str,
        ])
        .map_err(|e| format!("Failed to write CSV row: {}", e))?;
    }
    
    writer.flush()
        .map_err(|e| format!("Failed to flush CSV file: {}", e))?;
    
    Ok(())
}

/// Save smoothed results to CSV file
fn save_smoothed_to_csv(
    smoothed: &[SmoothedResult],
    cluster_numbers: &[i32],
    output_path: &Path,
) -> Result<(), String> {
    let mut writer = csv::Writer::from_path(output_path)
        .map_err(|e| format!("Failed to create CSV file {}: {}", output_path.display(), e))?;
    
    // Write header: frame, q_smooth, cluster_0_smooth, cluster_1_smooth, ...
    let mut header = vec!["frame".to_string(), "q_smooth".to_string()];
    for cluster_num in cluster_numbers {
        header.push(format!("cluster_{}_smooth", cluster_num));
    }
    writer.write_record(&header)
        .map_err(|e| format!("Failed to write CSV header: {}", e))?;
    
    // Write each row
    for result in smoothed {
        let mut row = vec![result.frame.to_string(), result.q_smooth.to_string()];
        for cluster_num in cluster_numbers {
            let value = result.cluster_smooth.get(cluster_num).copied().unwrap_or(0.0);
            row.push(value.to_string());
        }
        writer.write_record(&row)
            .map_err(|e| format!("Failed to write CSV row: {}", e))?;
    }
    
    writer.flush()
        .map_err(|e| format!("Failed to flush CSV file: {}", e))?;
    
    Ok(())
}

/// Load summary CSV file (from summarize_trajectory)
fn load_summary_csv(csv_path: &Path) -> Result<Vec<WindowSummary>, String> {
    let mut reader = csv::Reader::from_path(csv_path)
        .map_err(|e| format!("Failed to open summary CSV file {}: {}", csv_path.display(), e))?;
    
    let mut summaries = Vec::new();
    let headers = reader.headers()
        .map_err(|e| format!("Failed to read CSV headers: {}", e))?;
    
    // Find cluster columns
    let cluster_columns: Vec<(usize, i32)> = headers
        .iter()
        .enumerate()
        .filter_map(|(idx, col)| {
            if col.starts_with("cluster_") {
                let cluster_num_str = &col[8..]; // Skip "cluster_"
                cluster_num_str.parse::<i32>().ok().map(|num| (idx, num))
            } else {
                None
            }
        })
        .collect();
    
    let frame_col_idx = headers.iter()
        .position(|h| h == "frame")
        .ok_or("Missing 'frame' column in summary CSV")?;
    
    for record_result in reader.records() {
        let record = record_result
            .map_err(|e| format!("Failed to read CSV record: {}", e))?;
        
        let frame: i32 = record.get(frame_col_idx)
            .ok_or("Missing frame value")?
            .parse()
            .map_err(|e| format!("Failed to parse frame: {}", e))?;
        
        let mut cluster_means = HashMap::new();
        for (col_idx, cluster_num) in &cluster_columns {
            let value_str = record.get(*col_idx)
                .ok_or(format!("Missing cluster_{} value", cluster_num))?;
            let value: f64 = value_str.parse()
                .map_err(|e| format!("Failed to parse cluster_{} value '{}': {}", cluster_num, value_str, e))?;
            cluster_means.insert(*cluster_num, value);
        }
        
        summaries.push(WindowSummary {
            frame,
            cluster_means,
        });
    }
    
    Ok(summaries)
}

/// Save classification result to text file
fn save_classification_to_file(formation_order: &[i32], output_path: &Path) -> Result<(), String> {
    use std::fs::File;
    use std::io::Write;
    
    // Ensure output directory exists
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create output directory: {}", e))?;
    }
    
    let mut file = File::create(output_path)
        .map_err(|e| format!("Failed to create output file {}: {}", output_path.display(), e))?;
    
    // Write as comma-separated values on one line
    let line = formation_order
        .iter()
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(",");
    
    writeln!(file, "{}", line)
        .map_err(|e| format!("Failed to write to output file: {}", e))?;
    
    Ok(())
}

/// Load frame results from CSV file
fn load_results_from_csv(csv_path: &Path) -> Result<Vec<FrameResult>, String> {
    let mut reader = csv::Reader::from_path(csv_path)
        .map_err(|e| format!("Failed to open CSV file {}: {}", csv_path.display(), e))?;
    
    let mut results = Vec::new();
    
    for record_result in reader.records() {
        let record = record_result
            .map_err(|e| format!("Failed to read CSV record: {}", e))?;
        
        let frame: i32 = record.get(0)
            .ok_or("Missing frame column")?
            .parse()
            .map_err(|e| format!("Failed to parse frame: {}", e))?;
        
        let contacts: usize = record.get(1)
            .ok_or("Missing contacts column")?
            .parse()
            .map_err(|e| format!("Failed to parse contacts: {}", e))?;
        
        let q: f64 = record.get(2)
            .ok_or("Missing q column")?
            .parse()
            .map_err(|e| format!("Failed to parse q: {}", e))?;
        
        // Parse contact_list (skip for now, not needed for summary)
        // Parse clusters_filling
        let clusters_filling_str = record.get(4)
            .ok_or("Missing clusters_filling column")?;
        
        // Parse clusters_filling string like "{0: 0.5, 1: 0.3, ...}"
        let clusters_filling = parse_clusters_filling(clusters_filling_str)?;
        
        results.push(FrameResult {
            frame,
            contacts,
            q,
            contact_list: Vec::new(),  // Not needed for summary
            clusters_filling,
        });
    }
    
    Ok(results)
}

/// Parse clusters_filling string like "{0: 0.5, 1: 0.3, ...}"
fn parse_clusters_filling(s: &str) -> Result<HashMap<i32, f64>, String> {
    let mut clusters = HashMap::new();
    
    // Remove outer braces
    let s = s.trim();
    if !s.starts_with('{') || !s.ends_with('}') {
        return Err(format!("Invalid clusters_filling format: {}", s));
    }
    
    let s = &s[1..s.len()-1];  // Remove { and }
    
    if s.is_empty() {
        return Ok(clusters);
    }
    
    // Split by comma and parse each key:value pair
    for pair in s.split(',') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        
        let parts: Vec<&str> = pair.split(':').collect();
        if parts.len() != 2 {
            return Err(format!("Invalid cluster pair format: {}", pair));
        }
        
        let cluster_num: i32 = parts[0].trim().parse()
            .map_err(|e| format!("Failed to parse cluster number '{}': {}", parts[0], e))?;
        
        let fraction: f64 = parts[1].trim().parse()
            .map_err(|e| format!("Failed to parse fraction '{}': {}", parts[1], e))?;
        
        clusters.insert(cluster_num, fraction);
    }
    
    Ok(clusters)
}

/// Save window summaries to CSV file
fn save_summary_to_csv(
    summaries: &[WindowSummary],
    cluster_numbers: &[i32],
    output_path: &Path,
) -> Result<(), String> {
    let mut writer = csv::Writer::from_path(output_path)
        .map_err(|e| format!("Failed to create CSV file {}: {}", output_path.display(), e))?;
    
    // Write header: frame, cluster_0, cluster_1, ...
    let mut header = vec!["frame".to_string()];
    for cluster_num in cluster_numbers {
        header.push(format!("cluster_{}", cluster_num));
    }
    writer.write_record(&header)
        .map_err(|e| format!("Failed to write CSV header: {}", e))?;
    
    // Write each row
    for summary in summaries {
        let mut row = vec![summary.frame.to_string()];
        for cluster_num in cluster_numbers {
            let value = summary.cluster_means.get(cluster_num).copied().unwrap_or(0.0);
            row.push(value.to_string());
        }
        writer.write_record(&row)
            .map_err(|e| format!("Failed to write CSV row: {}", e))?;
    }
    
    writer.flush()
        .map_err(|e| format!("Failed to flush CSV file: {}", e))?;
    
    Ok(())
}

