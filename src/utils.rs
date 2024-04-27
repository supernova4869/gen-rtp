use std::io;
use std::path::{Path, PathBuf};
use regex::Regex;
use std::str::FromStr;
use std::fmt::Debug;

pub fn read_file() -> String {
    let inp = get_input("".to_string());
    let inp: String = match inp.starts_with("\"") {
        true => inp[1..inp.len() - 1].to_string(),
        false => inp.to_string()
    };
    let re = Regex::new(r"\\").unwrap();
    re.replace_all(&inp, "/").to_string()
}

pub fn get_stemname(fname: &str) -> String {
    let file = Path::new(fname);
    file.file_stem().unwrap().to_str().unwrap().to_string()
}

pub fn get_parent_path(fname: &str) -> PathBuf {
    let file = Path::new(fname);
    Path::new(file.parent().unwrap()).to_owned()
}

pub fn atrange2atlist(atom_selection_str: &str) -> Vec<usize> {
    let mut selection_range: Vec<usize> = vec![];
    if atom_selection_str.trim().is_empty() {
        return selection_range;
    }
    let atom_selection_str = atom_selection_str.replace(" ", "");
    let sub_selections: Vec<&str> = atom_selection_str.split(',').collect();
    for s in sub_selections {
        if s.contains('-') {
            let r: Vec<&str> = s.split('-').collect();
            let l: usize = r[0].parse().expect("Range lower bound not int");
            let u: usize = r[1].parse().expect("Range upper bound not int");
            selection_range.append(&mut (l..=u).collect());
        } else {
            let s: usize = s.parse().expect("Index not int");
            selection_range.push(s);
        }
    }
    return selection_range
}

pub fn get_input<T: FromStr>(default: T) -> T where <T as FromStr>::Err: Debug {
    let mut inp: String = String::new();
    io::stdin().read_line(&mut inp).expect("Failed to read line");
    let inp = inp.trim();
    match inp.is_empty() {
        true => default,
        false => inp.parse().expect("Failed to parse input")
    }
}

pub fn get_exclude_atoms() -> (Vec<usize>, Vec<usize>, 
                               Option<usize>, Option<usize>, 
                               Option<String>, Option<String>,
                               Option<usize>, Option<usize>, 
                               Option<String>, Option<String>) {
    println!("Input atoms id of the previous residue, e.g., 1-3, 5 (leave blank if it is the first residue): ");
    let prev_atoms = get_input("".to_string());
    let prev_atoms = atrange2atlist(prev_atoms.as_str());
    let (prev_con_atom, prev_con_atom_name, prev_adj_atom, prev_adj_atom_name) = match prev_atoms.is_empty() {
        false => {
            println!("Connection atom id of the previous residue (default: {}): ", prev_atoms[0]);
            let prev_con_atom = get_input(prev_atoms[0]);
            println!("Rename connection atom name to (default: -C): ");
            let prev_atom_name = get_input("-C".to_string());
            println!("Connection atom id of the current residue to previous (default: {})): ", prev_atoms[prev_atoms.len() - 1] + 1);
            let prev_adj_atom = get_input(prev_atoms[prev_atoms.len() - 1] + 1);
            println!("Rename connection atom name to (default: N): ");
            let prev_adj_atom_name = get_input("N".to_string());
            (Some(prev_con_atom), Some(prev_atom_name), Some(prev_adj_atom), Some(prev_adj_atom_name))
        },
        true => {
            (None, None, None, None)
        }
    };
    println!("Input atoms id of the next residue, e.g., 1-3, 5 (leave blank if it is the last residue): ");
    let next_atoms = get_input("".to_string());
    let next_atoms = atrange2atlist(next_atoms.as_str());
    let (next_con_atom, next_atom_name, next_adj_atom, next_adj_atom_name) = match next_atoms.is_empty() {
        false => {
            println!("Connection atom id of the next residue (default: {}): ", next_atoms[0]);
            let next_con_atom = get_input(next_atoms[0]);
            println!("Rename connection atom to (default: +N): ");
            let next_con_atom_name = get_input("+N".to_string());
            println!("Connection atom id of the current residue to next (default: {})): ", next_atoms[next_atoms.len() - 1] + 1);
            let next_adj_atom = get_input(next_atoms[0] - 1);
            println!("Rename connection atom to (default: C): ");
            let next_adj_atom_name = get_input("C".to_string());
            (Some(next_con_atom), Some(next_con_atom_name), Some(next_adj_atom), Some(next_adj_atom_name))
        },
        true => {
            (None, None, None, None)
        }
    };
    (prev_atoms, next_atoms, 
        prev_con_atom, next_con_atom, 
        prev_con_atom_name, next_atom_name, 
        prev_adj_atom, next_adj_atom,
        prev_adj_atom_name, next_adj_atom_name)
}