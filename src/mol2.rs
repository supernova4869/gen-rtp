use std::{fs, io::Write};
use crate::hdb::HDBItem;
use crate::hdb::{get_adj_h_id, get_adj_heavy_id, get_htype_from_heavy_atom};

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct MOL2 {
    mol: Molecule,
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
        let sys_name = Some(mol2_content[mol_ln + 1].trim().to_string());
        let num: Vec<i32> = mol2_content[mol_ln + 2].trim().split_whitespace().map(|s| s.parse().unwrap()).collect();
        let at_num = Some(num[0]);
        let bond_num = Some(num[1]);
        let sub_struct_num = Some(num[2]);
        let prop_num = Some(num[3]);
        let set_num = Some(num[4]);
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
        for &at_line in &mol2_content[atom_ln + 1 .. atom_ln + 1 + at_num.unwrap() as usize] {
            atoms.push(Atom::from(at_line));
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
        for &bond_line in &mol2_content[bond_ln + 1 .. bond_ln + 1 + bond_num.unwrap() as usize] {
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
    sys_name: String,
    at_num: i32,
    bond_num: i32,
    sub_struct_num: i32,
    prop_num: i32,
    set_num: i32,
    sys_type: String,
    at_charge: String,
}

impl Molecule {
    pub fn new(sys_name: Option<String>, at_num: Option<i32>, bond_num: Option<i32>, sub_struct_num: Option<i32>, prop_num: Option<i32>, 
               set_num: Option<i32>, sys_type: Option<String>, at_charge: Option<String>) -> Molecule {
        let sys_name = sys_name.unwrap_or("MOL".to_string());
        let at_num = at_num.unwrap_or(0);
        let bond_num = bond_num.unwrap_or(0);
        let sub_struct_num = sub_struct_num.unwrap_or(0);
        let prop_num = prop_num.unwrap_or(0);
        let set_num = set_num.unwrap_or(0);
        let sys_type = sys_type.unwrap_or("SMALL".to_string());
        let at_charge = at_charge.unwrap_or("USER_CHARGES".to_string());
        Molecule {
            sys_name, at_num, bond_num, sub_struct_num, prop_num, set_num, sys_type,at_charge
        }
    }

    fn output(&self) -> String {
        "@<TRIPOS>MOLECULE\n".to_string() + format!("{}\n{:5}{:6}     1 0 0\n{}\n{}\n\n", 
            self.sys_name, self.at_num, self.bond_num, self.sys_type, self.at_charge).as_str()
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct Atom {
    pub atom_id: i32,
    pub atom_name: String,
    x: f64,
    y: f64,
    z: f64,
    at: String,
    sub_struct_id: Option<i32>,
    sub_struct_name: Option<String>,
    atom_charge: Option<f64>,
    pub element: String
}

impl Atom {
    fn from(line: &str) -> Atom {
        let line: Vec<&str> = line.trim().split_whitespace().collect();
        let atom_id: i32 = line[0].parse().unwrap();
        let atom_name: String = line[1].to_string();
        let x: f64 = line[2].parse().unwrap();
        let y: f64 = line[3].parse().unwrap();
        let z: f64 = line[4].parse().unwrap();
        let at: String = line[5].to_string();
        let sub_struct_id: Option<i32> = Some(line[6].parse().unwrap());
        let sub_struct_name: Option<String> = Some(line[7].to_string());
        let atom_charge: Option<f64> = Some(line[8].parse().unwrap());
        let element: Vec<&str> = at.split(".").collect();
        let element = element[0].to_string();
        Atom {
            atom_id, atom_name, x, y, z, at, sub_struct_id, sub_struct_name, atom_charge, element
        }
    }

    fn output(&self) -> String {
        format!("{:7} {:10}{:12.4}{:12.4}{:12.4} {:7}{:3} {:9}{:8.4}\n", 
            self.atom_id, self.atom_name, self.x, self.y, self.z, self.at, 
            self.sub_struct_id.unwrap_or(0), self.sub_struct_name.to_owned().unwrap(), self.atom_charge.unwrap_or(0.0))
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct Bond {
    bond_id: i32,
    pub a1: i32,
    pub a2: i32,
    bt: String,
}

impl Bond {
    fn from(line: &str) -> Bond {
        let line: Vec<&str> = line.trim().split_whitespace().collect();
        let bond_id: i32 = line[0].parse().unwrap();
        let a1: i32 = line[1].parse().unwrap();
        let a2: i32 = line[2].parse().unwrap();
        let bt: String = line[3].to_string();
        Bond {
            bond_id, a1, a2, bt
        }
    }

    fn output(&self) -> String {
        format!("{:6}{:5}{:5} {}\n", 
            self.bond_id, self.a1, self.a2, self.bt)
    }
}

impl MOL2 {
    pub fn output(&self, outfile: &str) {
        // TODO: 输出文件
        let mut file = fs::File::create(outfile).unwrap();
        
        file.write_all(b"#	Created by:	gen-rtp\n#	Author:	Supernova\n\n").unwrap();
        file.write_all(self.mol.output().as_bytes()).unwrap();
        file.write_all(b"@<TRIPOS>ATOM\n").unwrap();
        for a in &self.atoms {
            file.write_fmt(format_args!("{}", a.output())).unwrap();
        }
        file.write_all(b"@<TRIPOS>BOND\n").unwrap();
        for b in &self.bonds {
            file.write_fmt(format_args!("{}", b.output())).unwrap();
        }
        
        println!("Written to {}", outfile);
    }

    pub fn top2hdb(&self, out: &str, exclude_n: Vec<i32>, exclude_c: Vec<i32>) {
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
        outfile.write_all(format!("{:5}{}\n", self.mol.sys_name, items.len()).as_bytes()).unwrap();
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
