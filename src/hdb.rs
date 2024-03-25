pub struct HDBItem {
    pub h_num: i32,
    pub h_type: i32,
    pub h_atom: String,
    pub heavy_atoms: Vec<String>
}

impl HDBItem {
    pub fn new(h_num: i32, h_type: i32, h_atom: String, heavy_atoms: Vec<String>) -> HDBItem {
        HDBItem{h_num, h_type, h_atom, heavy_atoms}
    }
}
