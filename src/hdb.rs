use crate::mol2::MOL2;

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

pub fn get_htype_from_heavy_atom(mol2: &MOL2, ref_id: usize) -> i32 {
    // 相连原子
    let adj = get_adj_atoms_id(mol2, ref_id);
    // 相连H
    let adj_h = filter_h(mol2, &adj);
    match (adj.len(), adj_h.len()) {
        (3, 1) => 1,
        (2, 1) => 2,
        (3, 2) => 3,
        (4, 3) => 4,
        (4, 1) => 5,
        (4, 2) => 6,
        _ => 0,
    }
}

pub fn get_adj_atoms_id(mol2: &MOL2, ref_id: usize) -> Vec<usize> {
    mol2.bonds.iter()
        .filter(|&b| (b.a1 as i32 - ref_id as i32) * (b.a2 as i32 - ref_id as i32) == 0)
        .map(|b| b.a1 + b.a2 - ref_id)
        .collect()
}

pub fn get_adj_h_id(mol2: &MOL2, ref_id: usize) -> Vec<usize> {
    let adj = get_adj_atoms_id(mol2, ref_id);
    filter_h(mol2, &adj)
}

pub fn get_adj_heavy_id(mol2: &MOL2, ref_id: usize) -> Vec<usize> {
    let adj = get_adj_atoms_id(mol2, ref_id);
    filter_heavy(mol2, &adj)
}

fn filter_h(mol2: &MOL2, atoms_id: &Vec<usize>) -> Vec<usize> {
    atoms_id.iter()
        .filter(|&&a| mol2.atoms[a as usize - 1].element.eq("H"))
        .cloned()
        .collect()
}

fn filter_heavy(mol2: &MOL2, atoms_id: &Vec<usize>) -> Vec<usize> {
    atoms_id.iter()
        .filter(|&&a| mol2.atoms[a as usize - 1].element.ne("H"))
        .cloned()
        .collect()
}