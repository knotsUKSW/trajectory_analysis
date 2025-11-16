use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use std::path::{Path, PathBuf};

use crate::contacts::load_contacts_from_csv;
use crate::trajectory::{PdbTrajectory, Trajectory};

/// Python binding for read_trajectory function
#[pyfunction]
#[pyo3(signature = (trajectory_file, contacts_file, cutoff_distance=1.2, max_frames=None, output_csv=None))]
fn read_trajectory(
    py: Python<'_>,
    trajectory_file: &str,
    contacts_file: &str,
    cutoff_distance: f64,
    max_frames: Option<usize>,
    output_csv: Option<&str>,
) -> PyResult<PyObject> {
    // Load contacts
    let contacts = load_contacts_from_csv(contacts_file)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to load contacts: {}", e)))?;
    
    // Create trajectory reader
    let trajectory = PdbTrajectory::new(trajectory_file);
    
    // Determine output path - auto-generate if None
    let output_path = if let Some(csv_path) = output_csv {
        Some(PathBuf::from(csv_path))
    } else {
        // Auto-generate from trajectory file path
        let traj_path = Path::new(trajectory_file);
        let traj_dir = traj_path.parent().unwrap_or(Path::new("."));
        let base_name = traj_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("trajectory");
        let output_file = format!("{}_parsed.csv", base_name);
        Some(traj_dir.join(output_file))
    };
    
    // Read trajectory
    let results = trajectory.read_trajectory(
        &contacts,
        cutoff_distance,
        max_frames,
        output_path.as_deref(),
    )
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to read trajectory: {}", e)))?;
    
    // Convert results to Python list of dicts
    let py_results = PyList::empty_bound(py);
    
    for result in results {
        let py_dict = PyDict::new_bound(py);
        
        // Convert contact_list to Python list of tuples
        let py_contact_list = PyList::empty_bound(py);
        for (i, j) in &result.contact_list {
            let py_tuple = PyTuple::new_bound(py, &[i.into_py(py), j.into_py(py)]);
            py_contact_list.append(py_tuple)?;
        }
        
        // Convert clusters_filling to Python dict
        let py_clusters = PyDict::new_bound(py);
        for (cluster, fraction) in &result.clusters_filling {
            py_clusters.set_item(cluster, fraction)?;
        }
        
        // Build the dictionary
        py_dict.set_item("frame", result.frame)?;
        py_dict.set_item("contacts", result.contacts)?;
        py_dict.set_item("q", result.q)?;
        py_dict.set_item("contact_list", py_contact_list)?;
        py_dict.set_item("clusters_filling", py_clusters)?;
        
        py_results.append(py_dict)?;
    }
    
    Ok(py_results.into())
}

/// Python binding for summarize_trajectory function
#[pyfunction]
#[pyo3(signature = (trajectory_file, window_size=10000, cutoff=None, output_csv=None))]
fn summarize_trajectory(
    py: Python<'_>,
    trajectory_file: &str,
    window_size: usize,
    cutoff: Option<f64>,
    output_csv: Option<&str>,
) -> PyResult<PyObject> {
    // Create trajectory reader
    let trajectory = PdbTrajectory::new(trajectory_file);
    
    // Determine output path - auto-generate if None
    let output_path = if let Some(csv_path) = output_csv {
        Some(PathBuf::from(csv_path))
    } else {
        None  // Will be auto-generated in summarize_trajectory
    };
    
    // Summarize trajectory (loads from CSV if results not provided)
    let summaries = trajectory.summarize_trajectory(
        None,  // Load from CSV
        window_size,
        cutoff,
        output_path.as_deref(),
    )
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to summarize trajectory: {}", e)))?;
    
    // Convert summaries to Python list of dicts
    let py_summaries = PyList::empty_bound(py);
    
    for summary in summaries {
        let py_dict = PyDict::new_bound(py);
        
        // Convert cluster_means to Python dict
        let py_clusters = PyDict::new_bound(py);
        for (cluster, mean) in &summary.cluster_means {
            py_clusters.set_item(format!("cluster_{}", cluster), mean)?;
        }
        
        // Build the dictionary
        py_dict.set_item("frame", summary.frame)?;
        py_dict.set_item("cluster_means", py_clusters)?;
        
        py_summaries.append(py_dict)?;
    }
    
    Ok(py_summaries.into())
}

/// Python binding for smooth function
#[pyfunction]
#[pyo3(signature = (trajectory_file, window_size=100, output_csv=None))]
fn smooth(
    py: Python<'_>,
    trajectory_file: &str,
    window_size: usize,
    output_csv: Option<&str>,
) -> PyResult<PyObject> {
    // Create trajectory reader
    let trajectory = PdbTrajectory::new(trajectory_file);
    
    // Determine output path - auto-generate if None
    let output_path = if let Some(csv_path) = output_csv {
        Some(PathBuf::from(csv_path))
    } else {
        None  // Will be auto-generated in smooth
    };
    
    // Smooth trajectory (loads from CSV if results not provided)
    let smoothed = trajectory.smooth(
        None,  // Load from CSV
        window_size,
        output_path.as_deref(),
    )
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to smooth trajectory: {}", e)))?;
    
    // Convert smoothed results to Python list of dicts
    let py_smoothed = PyList::empty_bound(py);
    
    for result in smoothed {
        let py_dict = PyDict::new_bound(py);
        
        // Convert cluster_smooth to Python dict
        let py_clusters = PyDict::new_bound(py);
        for (cluster, smooth_val) in &result.cluster_smooth {
            py_clusters.set_item(format!("cluster_{}_smooth", cluster), smooth_val)?;
        }
        
        // Build the dictionary
        py_dict.set_item("frame", result.frame)?;
        py_dict.set_item("q_smooth", result.q_smooth)?;
        py_dict.set_item("cluster_smooth", py_clusters)?;
        
        py_smoothed.append(py_dict)?;
    }
    
    Ok(py_smoothed.into())
}

/// Python binding for classify function
#[pyfunction]
#[pyo3(signature = (trajectory_file, summary_csv=None, output_path=None))]
fn classify(
    py: Python<'_>,
    trajectory_file: &str,
    summary_csv: Option<&str>,
    output_path: Option<&str>,
) -> PyResult<PyObject> {
    // Create trajectory reader
    let trajectory = PdbTrajectory::new(trajectory_file);
    
    // Determine paths
    let summary_path = summary_csv.map(|p| PathBuf::from(p));
    let output_path_buf = output_path.map(|p| PathBuf::from(p));
    
    // Classify trajectory
    let formation_order = trajectory.classify(
        summary_path.as_deref(),
        output_path_buf.as_deref(),
    )
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to classify trajectory: {}", e)))?;
    
    // Convert to Python list
    let py_order = PyList::empty_bound(py);
    for cluster_num in formation_order {
        py_order.append(cluster_num)?;
    }
    
    Ok(py_order.into())
}

/// Python module definition
#[pymodule]
fn folding_analysis_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_trajectory, m)?)?;
    m.add_function(wrap_pyfunction!(summarize_trajectory, m)?)?;
    m.add_function(wrap_pyfunction!(smooth, m)?)?;
    m.add_function(wrap_pyfunction!(classify, m)?)?;
    m.add("__doc__", "Folding analysis Rust library with Python bindings")?;
    Ok(())
}

