use serde::Deserialize;

/// Contact information from CSV file
#[derive(Debug, Clone, Deserialize)]
pub struct Contact {
    pub i: i32,
    pub j: i32,
    pub r: f64,
    pub cluster: i32,
}

/// Load contacts from CSV file
pub fn load_contacts_from_csv(csv_path: &str) -> Result<Vec<Contact>, String> {
    let mut reader = csv::Reader::from_path(csv_path)
        .map_err(|e| format!("Failed to open contacts CSV file {}: {}", csv_path, e))?;
    
    let mut contacts = Vec::new();
    
    for result in reader.deserialize() {
        let contact: Contact = result
            .map_err(|e| format!("Failed to parse contact from CSV: {}", e))?;
        contacts.push(contact);
    }
    
    Ok(contacts)
}

