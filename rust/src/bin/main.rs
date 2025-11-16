use clap::{Parser, Subcommand};
use folding_analysis_rs::{load_contacts_from_csv, load_results_from_csv, PdbTrajectory, Trajectory};
use std::path::{Path, PathBuf};

/// Command-line tool for analyzing protein folding trajectories
#[derive(Parser)]
#[command(name = "folding-analysis")]
#[command(about = "Analyze protein folding trajectories from PDB files", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Read trajectory and calculate native contact formation
    Read {
        /// Path to the PDB trajectory file
        #[arg(short, long)]
        trajectory: PathBuf,
        
        /// Path to the contacts CSV file (must contain columns: i, j, r, cluster)
        #[arg(short, long)]
        contacts: PathBuf,
        
        /// Multiplier for native distance cutoff (default: 1.2)
        #[arg(long, default_value_t = 1.2)]
        cutoff_distance: f64,
        
        /// Maximum number of frames to process (default: all frames)
        #[arg(long)]
        max_frames: Option<usize>,
        
        /// Output CSV path (default: auto-generated from trajectory path)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    
    /// Summarize trajectory by calculating mean cluster filling fractions in windows
    Summarize {
        /// Path to the parsed trajectory CSV (from read command)
        #[arg(short, long)]
        input: PathBuf,
        
        /// Window size for summarization (default: 10000)
        #[arg(short, long, default_value_t = 10000)]
        window_size: usize,
        
        /// Optional cutoff for binary conversion (values >= cutoff become 1, values < cutoff become 0)
        #[arg(long)]
        cutoff: Option<f64>,
        
        /// Output CSV path (default: auto-generated from input path)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    
    /// Smooth trajectory data by calculating running averages
    Smooth {
        /// Path to the parsed trajectory CSV (from read command)
        #[arg(short, long)]
        input: PathBuf,
        
        /// Window size for running average (default: 100)
        #[arg(short, long, default_value_t = 100)]
        window_size: usize,
        
        /// Output CSV path (default: auto-generated from input path)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    
    /// Classify trajectory by determining cluster formation order
    Classify {
        /// Path to the summary CSV file (binary representation from summarize with cutoff)
        #[arg(short, long)]
        input: PathBuf,
        
        /// Output text file path (default: auto-generated from input path)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

fn main() {
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Read {
            trajectory,
            contacts,
            cutoff_distance,
            max_frames,
            output,
        } => {
            println!("Reading trajectory: {:?}", trajectory);
            println!("Using contacts: {:?}", contacts);
            println!("Cutoff distance: {}", cutoff_distance);
            
            // Load contacts
            let contacts_vec = match load_contacts_from_csv(contacts.to_str().unwrap()) {
                Ok(c) => {
                    println!("✅ Loaded {} contacts", c.len());
                    c
                }
                Err(e) => {
                    eprintln!("❌ Error loading contacts: {}", e);
                    std::process::exit(1);
                }
            };
            
            // Create trajectory reader
            let traj = PdbTrajectory::new(&trajectory);
            
            // Read trajectory
            match traj.read_trajectory(
                &contacts_vec,
                cutoff_distance,
                max_frames,
                output.as_deref(),
            ) {
                Ok(results) => {
                    println!("✅ Successfully processed {} frames", results.len());
                    if let Some(output_path) = output {
                        println!("📄 Results saved to: {:?}", output_path);
                    } else {
                        // Auto-generated path
                        let base = trajectory.file_stem().unwrap().to_str().unwrap();
                        let dir = trajectory.parent().unwrap_or(std::path::Path::new("."));
                        let auto_path = dir.join(format!("{}_parsed.csv", base));
                        println!("📄 Results saved to: {:?}", auto_path);
                    }
                }
                Err(e) => {
                    eprintln!("❌ Error reading trajectory: {}", e);
                    std::process::exit(1);
                }
            }
        }
        
        Commands::Summarize {
            input,
            window_size,
            cutoff,
            output,
        } => {
            println!("Summarizing trajectory: {:?}", input);
            println!("Window size: {}", window_size);
            if let Some(c) = cutoff {
                println!("Binary cutoff: {}", c);
            }
            
            // Load results from CSV
            let results = match load_results_from_csv(&input) {
                Ok(r) => {
                    println!("✅ Loaded {} frames from CSV", r.len());
                    r
                }
                Err(e) => {
                    eprintln!("❌ Error loading CSV: {}", e);
                    std::process::exit(1);
                }
            };
            
            // Create trajectory with input path (for output path generation)
            let traj = PdbTrajectory::new(&input);
            
            match traj.summarize_trajectory(Some(&results), window_size, cutoff, output.as_deref()) {
                Ok(summaries) => {
                    println!("✅ Successfully summarized {} windows", summaries.len());
                    if let Some(output_path) = output {
                        println!("📄 Summary saved to: {:?}", output_path);
                    } else {
                        // Auto-generated path
                        let base = input.file_stem().unwrap().to_str().unwrap();
                        let dir = input.parent().unwrap_or(std::path::Path::new("."));
                        let suffix = if cutoff.is_some() {
                            format!("_summary_binary_{}.csv", cutoff.unwrap())
                        } else {
                            "_summary.csv".to_string()
                        };
                        let auto_path = dir.join(format!("{}{}", base, suffix));
                        println!("📄 Summary saved to: {:?}", auto_path);
                    }
                }
                Err(e) => {
                    eprintln!("❌ Error summarizing trajectory: {}", e);
                    std::process::exit(1);
                }
            }
        }
        
        Commands::Smooth {
            input,
            window_size,
            output,
        } => {
            println!("Smoothing trajectory: {:?}", input);
            println!("Window size: {}", window_size);
            
            // Load results from CSV
            let results = match load_results_from_csv(&input) {
                Ok(r) => {
                    println!("✅ Loaded {} frames from CSV", r.len());
                    r
                }
                Err(e) => {
                    eprintln!("❌ Error loading CSV: {}", e);
                    std::process::exit(1);
                }
            };
            
            // Create trajectory with input path (for output path generation)
            let traj = PdbTrajectory::new(&input);
            
            match traj.smooth(Some(&results), window_size, output.as_deref()) {
                Ok(smoothed) => {
                    println!("✅ Successfully smoothed {} frames", smoothed.len());
                    if let Some(output_path) = output {
                        println!("📄 Smoothed data saved to: {:?}", output_path);
                    } else {
                        // Auto-generated path
                        let base = input.file_stem().unwrap().to_str().unwrap();
                        let dir = input.parent().unwrap_or(std::path::Path::new("."));
                        let auto_path = dir.join(format!("{}_{}_smoothed.csv", base, window_size));
                        println!("📄 Smoothed data saved to: {:?}", auto_path);
                    }
                }
                Err(e) => {
                    eprintln!("❌ Error smoothing trajectory: {}", e);
                    std::process::exit(1);
                }
            }
        }
        
        Commands::Classify {
            input,
            output,
        } => {
            println!("Classifying trajectory: {:?}", input);
            
            let traj = PdbTrajectory::new(&input);
            
            match traj.classify(Some(&input), output.as_deref()) {
                Ok(order) => {
                    println!("✅ Cluster formation order: {:?}", order);
                    if let Some(output_path) = output {
                        println!("📄 Classification saved to: {:?}", output_path);
                    } else {
                        // Auto-generated path
                        let base = input.file_stem().unwrap().to_str().unwrap();
                        let dir = input.parent().unwrap_or(std::path::Path::new("."));
                        let auto_path = dir.join(format!("{}_class.txt", base));
                        println!("📄 Classification saved to: {:?}", auto_path);
                    }
                }
                Err(e) => {
                    eprintln!("❌ Error classifying trajectory: {}", e);
                    std::process::exit(1);
                }
            }
        }
    }
}

