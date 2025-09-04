use std::process::exit;
use std::{fs, io::Write};
use std::path::Path;
use crate::hdb::HDBItem;
use crate::hdb::{get_adj_h_id, get_adj_heavy_id, get_htype_from_heavy_atom};
use std::fmt::{self, Debug, Display};

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct MOL2 {
    pub mol: Molecule,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

impl MOL2 {
    pub fn new(mol: Molecule, atoms: Vec<Atom>, bonds: Vec<Bond>) -> MOL2 {
        MOL2 { mol, atoms, bonds }
    }
    pub fn from(file: &str) -> MOL2 {
        // 读取文件
        let mol2_content = fs::read_to_string(file).unwrap();
        let mut mol2_content: Vec<&str> = mol2_content.split("\n").collect();
        mol2_content.iter_mut().for_each(|s| *s = s.trim());

        // Molecule定位
        let mol_ln = mol2_content.iter().enumerate().find_map(|(index, &s)| {
            if s.eq("@<TRIPOS>MOLECULE") {
                Some(index)
            } else {
                None
            }
        }).unwrap();

        // Molecule字段
        let sys_name = Path::new(file).file_stem().unwrap().to_str().unwrap();
        let num: Vec<i32> = mol2_content[mol_ln + 2].trim().split_whitespace().map(|s| s.parse().unwrap()).collect();
        let at_num = num.get(0);
        if at_num.is_none() {
            println!("Error: atom number is 0.");
            exit(0);
        }
        let bond_num = num.get(1);
        if bond_num.is_none() {
            println!("Error: bond number is 0.");
            exit(0);
        }
        let sub_struct_num = num.get(2);
        let prop_num = num.get(3);
        let set_num = num.get(4);
        let sys_type = Some(mol2_content[mol_ln + 3].trim().to_string());
        let at_charge = Some(mol2_content[mol_ln + 4].trim().to_string());
        let mol = Molecule::new(sys_name, at_num, bond_num, sub_struct_num, prop_num, set_num, sys_type, at_charge);

        // Atom字段
        let atom_ln = mol2_content.iter().enumerate().find_map(|(index, &s)| {
            if s.eq("@<TRIPOS>ATOM") {
                Some(index)
            } else {
                None
            }
        }).unwrap();
        let mut atoms: Vec<Atom> = Vec::new();
        for &at_line in &mol2_content[atom_ln + 1 .. atom_ln + 1 + *at_num.unwrap() as usize] {
            let mut atom = Atom::from(at_line);
            atom.sub_struct_name = sys_name.to_string();
            atoms.push(atom);
        }
        
        // Bond字段
        let bond_ln = mol2_content.iter().enumerate().find_map(|(index, &s)| {
            if s.eq("@<TRIPOS>BOND") {
                Some(index)
            } else {
                None
            }
        }).unwrap();
        let mut bonds: Vec<Bond> = Vec::new();
        for &bond_line in &mol2_content[bond_ln + 1 .. bond_ln + 1 + *bond_num.unwrap() as usize] {
            bonds.push(Bond::from(bond_line));
        }

        // total mol2
        MOL2::new(mol, atoms, bonds)
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct Molecule {
    pub sys_name: String,
    at_num: i32,
    bond_num: i32,
    sub_struct_num: i32,
    prop_num: i32,
    set_num: i32,
    sys_type: String,
    at_charge: String,
}

impl Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "@<TRIPOS>MOLECULE\n{}\n{:5}{:6}     1 0 0\n{}\n{}\n\n", 
            self.sys_name, self.at_num, self.bond_num, self.sys_type, self.at_charge)
    }
}

impl Molecule {
    pub fn new(sys_name: &str, at_num: Option<&i32>, bond_num: Option<&i32>, sub_struct_num: Option<&i32>, prop_num: Option<&i32>, 
               set_num: Option<&i32>, sys_type: Option<String>, at_charge: Option<String>) -> Molecule {
        let sys_name = sys_name.to_string();
        let at_num = *at_num.unwrap_or(&0);
        let bond_num = *bond_num.unwrap_or(&0);
        let sub_struct_num = *sub_struct_num.unwrap_or(&0);
        let prop_num = *prop_num.unwrap_or(&0);
        let set_num = *set_num.unwrap_or(&0);
        let sys_type = sys_type.unwrap_or("SMALL".to_string());
        let at_charge = at_charge.unwrap_or("USER_CHARGES".to_string());
        Molecule {
            sys_name, at_num, bond_num, sub_struct_num, prop_num, set_num, sys_type,at_charge
        }
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct Atom {
    pub atom_id: usize,
    pub atom_name: String,
    x: f64,
    y: f64,
    z: f64,
    at: String,
    sub_struct_id: i32,
    sub_struct_name: String,
    atom_charge: f64,
    pub element: String
}

impl Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:7} {:10}{:12.4}{:12.4}{:12.4} {:7}{:3} {:9}{:8.4}\n", 
            self.atom_id, self.atom_name, self.x, self.y, self.z, self.at, 
            self.sub_struct_id, self.sub_struct_name, self.atom_charge)
    }
}

impl Atom {
    fn from(line: &str) -> Atom {
        let line: Vec<&str> = line.trim().split_whitespace().collect();
        let atom_id: usize = line[0].parse().unwrap();
        let atom_name: String = line[1].to_string();
        let x: f64 = line[2].parse().unwrap();
        let y: f64 = line[3].parse().unwrap();
        let z: f64 = line[4].parse().unwrap();
        let at: String = line[5].to_string();
        // 修改残基名为MOL, 残基编号为1
        let sub_struct_id = 1;
        let sub_struct_name = "MOL".to_string();
        let atom_charge = line.get(8);
        let atom_charge = atom_charge.and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let element: Vec<&str> = at.split(".").collect();
        let element = element[0].to_string();
        Atom {
            atom_id, atom_name, x, y, z, at, sub_struct_id, sub_struct_name, atom_charge, element
        }
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct Bond {
    bond_id: usize,
    pub a1: usize,
    pub a2: usize,
    bt: String,
}

impl Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:6}{:5}{:5} {}\n", 
            self.bond_id, self.a1, self.a2, self.bt)
    }
}

impl Bond {
    fn from(line: &str) -> Bond {
        let line: Vec<&str> = line.trim().split_whitespace().collect();
        let bond_id: usize = line[0].parse().unwrap();
        let a1: usize = line[1].parse().unwrap();
        let a2: usize = line[2].parse().unwrap();
        let bt: String = line[3].to_string();
        Bond {
            bond_id, a1, a2, bt
        }
    }
}

impl Display for MOL2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut out = format!("; Created by gen-rtp (https://github.com/supernova4869/gen-rtp)\n\n");
        out.push_str(format!("{}", self.mol).as_str());
        out.push_str("@<TRIPOS>ATOM\n");
        for a in &self.atoms {
            out.push_str(format!("{}", a).as_str());
        }
        out.push_str("@<TRIPOS>BOND\n");
        for b in &self.bonds {
            out.push_str(format!("{}", b).as_str());
        }
        write!(f, "{}", out)
    }
}

impl MOL2 {
    pub fn output(&self, outfile: &str) {
        let mut file = fs::File::create(outfile).unwrap();
        file.write_fmt(format_args!("{}", self)).unwrap();
        println!("Written to {}", outfile);
    }

    pub fn to_hdb(&mut self, out: &str,
        exclude_n: &Vec<usize>, exclude_c: &Vec<usize>,
        atom_n: Option<usize>, atom_c: Option<usize>,
        n_name: &Option<String>, c_name: &Option<String>,
        atom_adjn: Option<usize>, atom_adjc: Option<usize>,
        adjn_name: &Option<String>, adjc_name: &Option<String>) {
        // 前后残基中的原子名加前缀
        for atom in &mut self.atoms {
            if let Some(atom_n) = atom_n {
                if exclude_n.contains(&atom.atom_id) {
                    if atom.atom_id != atom_n {
                        atom.atom_name = "-".to_string() + &atom.atom_name;
                    } else {
                        atom.atom_name = n_name.as_ref().unwrap().to_string();
                    }
                }
            }
            if let Some(atom_c) = atom_c {
                if exclude_c.contains(&atom.atom_id) {
                    if atom.atom_id != atom_c {
                        atom.atom_name = "+".to_string() + &atom.atom_name;
                    } else {
                        atom.atom_name = c_name.as_ref().unwrap().to_string();
                    }
                }
            }
            if let Some(atom_adjn) = atom_adjn {
                if atom.atom_id == atom_adjn {
                    atom.atom_name = adjn_name.as_ref().unwrap().to_string();
                }
            }
            if let Some(atom_adjc) = atom_adjc {
                if atom.atom_id == atom_adjc {
                    atom.atom_name = adjc_name.as_ref().unwrap().to_string();
                }
            }
        }

        // atoms layout: H--i--j--k
        let mut items: Vec<HDBItem> = vec![];
        // 1. 找到所有的非排除重原子
        let atoms_heavy: Vec<&Atom> = self.atoms.iter()
            .filter(|&h| h.element.ne("H"))
            .filter(|&h| !exclude_n.contains(&h.atom_id))
            .filter(|&h| !exclude_c.contains(&h.atom_id))
            .collect();
        for atom_i in atoms_heavy {
            // 2. 找到和i相连的重原子j
            let atom_j = get_adj_heavy_id(self, atom_i.atom_id);
            // 3. 找到和j相连的重原子k, 排除i
            let mut atom_k = get_adj_heavy_id(self, atom_j[0]);
            atom_k.retain(|&a| a != atom_i.atom_id);
            // 4. 判断重原子连接的H类型
            let htype = get_htype_from_heavy_atom(self, atom_i.atom_id);
            let hs = get_adj_h_id(self, atom_i.atom_id);
            if hs.len() == 0 {
                continue;
            } else {
                let h = &self.atoms[hs[0] as usize - 1];
                let h_basename = self.get_hbasename(h);
                let cur_h = match htype {
                    // type 1, 环H/肽H
                    1 => HDBItem::new(1, 1, h_basename, 
                        vec![
                            atom_i.atom_name.to_string(), 
                            self.atoms[atom_j[0] as usize - 1].atom_name.to_string(), 
                            self.atoms[atom_j[1] as usize - 1].atom_name.to_string()
                        ]
                    ),
                    // type 2, 羟H
                    2 => HDBItem::new(1, 2, h_basename, 
                        vec![
                            atom_i.atom_name.to_string(), 
                            self.atoms[atom_j[0] as usize - 1].atom_name.to_string(), 
                            self.atoms[atom_k[0] as usize - 1].atom_name.to_string()
                        ]
                    ),
                    // type 3, 烯H/酰胺H
                    3 => HDBItem::new(2, 3, h_basename, 
                        vec![
                            atom_i.atom_name.to_string(), 
                            self.atoms[atom_j[0] as usize - 1].atom_name.to_string(), 
                            self.atoms[atom_k[0] as usize - 1].atom_name.to_string()
                        ]
                    ),
                    // type 4, 甲基H
                    4 => HDBItem::new(3, 4, h_basename, 
                        vec![
                            atom_i.atom_name.to_string(), 
                            self.atoms[atom_j[0] as usize - 1].atom_name.to_string(), 
                            self.atoms[atom_k[0] as usize - 1].atom_name.to_string()
                        ]
                    ),
                    // type 5, 特丁基H
                    5 => HDBItem::new(1, 5, h_basename, 
                        vec![
                            atom_i.atom_name.to_string(), 
                            self.atoms[atom_j[0] as usize - 1].atom_name.to_string(), 
                            self.atoms[atom_j[1] as usize - 1].atom_name.to_string(),
                            self.atoms[atom_j[2] as usize - 1].atom_name.to_string()
                        ]
                    ),
                    // type 6, 亚甲基H
                    6 => HDBItem::new(2, 6, h_basename, 
                        vec![
                            atom_i.atom_name.to_string(), 
                            self.atoms[atom_j[0] as usize - 1].atom_name.to_string(), 
                            self.atoms[atom_j[1] as usize - 1].atom_name.to_string()
                        ]
                    ),
                    _ => HDBItem::new(1, 1, "".to_string(), vec![])
                };
                items.push(cur_h);
            }
        }
        // 5. 输出文件
        let mut outfile = fs::File::create(out).unwrap();
        outfile.write_all(format!("{:5}    {}\n", self.mol.sys_name, items.len()).as_bytes()).unwrap();
        for item in &items {
            outfile.write_all(format!("{:<7}{:<7}{:7}", item.h_num, item.h_type, item.h_atom).as_bytes()).unwrap();
            for ha in &item.heavy_atoms {
                outfile.write_all(format!("{:7}", ha).as_bytes()).unwrap();
            }
            outfile.write_all(b"\n").unwrap();
        }

        println!("Finished writing rtp file to {}", out);
    }

    pub fn get_hbasename(&self, h: &Atom) -> String {
        let heavy = &self.atoms[get_adj_heavy_id(self, h.atom_id)[0] as usize - 1];
        format!("H{}", &heavy.atom_name[heavy.element.len()..])
    }
}
