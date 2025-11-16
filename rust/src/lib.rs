pub mod contacts;
pub mod structure;
pub mod trajectory;

#[cfg(feature = "python")]
pub mod python_bindings;

// Re-export commonly used types and traits
pub use contacts::{Contact, load_contacts_from_csv};
pub use structure::{Coordinate, FrameData};
pub use trajectory::{FrameResult, PdbTrajectory, SmoothedResult, Trajectory, TrajectoryData, WindowSummary};
