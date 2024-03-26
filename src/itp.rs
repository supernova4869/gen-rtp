use regex::Regex;
use std::hash::{Hash, Hasher};
use std::{collections::HashSet, fs, io::Write};
use std::str::FromStr;
use std::fmt::{self, Debug, Display};

use crate::mol2::MOL2;

pub struct TopolAtomtype {
    name: String,
    mass: f64,
    charge: f64,
    ptype: String,
    sigma: f64,
    epsilon: f64,
}

impl PartialEq for TopolAtomtype {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name 
            && (self.mass - other.mass).abs() < std::f64::EPSILON
            && (self.charge - other.charge).abs() < std::f64::EPSILON
            && self.ptype == other.ptype 
            && self.sigma == other.sigma 
            && self.epsilon == other.epsilon 
    }
}

impl Eq for TopolAtomtype {}

impl Hash for TopolAtomtype {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
        self.mass.to_bits().hash(state);
        self.charge.to_bits().hash(state);
        self.ptype.hash(state);
        self.sigma.to_bits().hash(state);
        self.epsilon.to_bits().hash(state);
    }
}

#[derive(Clone, Debug)]
pub struct TopolAtom {
    pub nr: usize,
    _type: String,
    resnr: i32,
    resname: String,
    pub atom: String,
    cgnr: i32,
    charge: f64,
    mass: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct TopolBond {
    ai: TopolAtom,
    aj: TopolAtom,
    funct: i32,
    c0: Option<f64>,
    c1: Option<f64>,
}

pub struct TopolPair {
    ai: TopolAtom,
    aj: TopolAtom,
    funct: i32,
}

pub struct TopolConstraint {
    ai: TopolAtom,
    aj: TopolAtom,
    funct: i32,
    cs: Option<Vec<f64>>,
}

pub struct TopolAngle {
    ai: TopolAtom,
    aj: TopolAtom,
    ak: TopolAtom,
    funct: i32,
    c0: Option<f64>,
    c1: Option<f64>,
}

pub struct TopolDihedral {
    ai: TopolAtom,
    aj: TopolAtom,
    ak: TopolAtom,
    al: TopolAtom,
    funct: i32,
    c0: Option<f64>,
    c1: Option<f64>,
    c2: Option<f64>,
}

pub struct TopolExclusion {
    ex_atoms: Vec<TopolAtom>,
}

pub struct Topol {
    atomtypes: HashSet<TopolAtomtype>,
    moleculetype: String,
    nrexcl: i32,
    pub atoms: Vec<TopolAtom>,
    bonds: Vec<TopolBond>,
    pairs: Vec<TopolPair>,
    constraints: Vec<TopolConstraint>,
    angles: Vec<TopolAngle>,
    dihedrals: Vec<TopolDihedral>,
    exclusions: Vec<TopolExclusion>,
}

fn get_atom_from_nr(atoms: &Vec<TopolAtom>, nr: usize) -> &TopolAtom {
    atoms.iter().find(|&a| a.nr == nr).unwrap()
}


// fn get_nr_from_atom(atoms: &Vec<TopolAtom>, resnr: i32, atname: &str) -> &TopolAtom {
//     atoms.iter().find(|&a| a.resnr == resnr && a.atom == atname).unwrap()
// }

impl Topol {
    pub fn from(file: &str, mol2: &MOL2) -> Topol {
        println!("Reading topology of {}...", file);
        let lines = fs::read_to_string(file).unwrap();
        let lines: Vec<&str> = lines.split("\n").collect();
        let re = Regex::new(r"\s*;.*").unwrap();
        let lines: Vec<String> = lines.iter()
            .map(|&s| re.replace(s, "").to_string())
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();

        // topol items
        let mut attypes: HashSet<TopolAtomtype> = HashSet::new();
        let mut mol = mol2.mol.sys_name.to_string();
        let mut nrexcl = 3;
        let mut atoms: Vec<TopolAtom> = vec![];
        let mut bonds: Vec<TopolBond> = vec![];
        let mut pairs: Vec<TopolPair> = vec![];
        let mut constraints: Vec<TopolConstraint> = vec![];
        let mut angles: Vec<TopolAngle> = vec![];
        let mut dihedrals: Vec<TopolDihedral> = vec![];
        let mut exclusions: Vec<TopolExclusion> = vec![];

        let mut cur_item = "";
        let re = Regex::new(r".*\[(.*)].*").unwrap();

        for line in &lines {
            if re.is_match(line) {
                let caps = re.captures(line).unwrap();
                cur_item = caps.get(1).unwrap().as_str().trim();
            }
            else {
                if cur_item.eq("atomtypes") {
                    attypes.insert(TopolAtomtype::from(line));
                } else if cur_item == "moleculetype" {
                    let paras: Vec<&str> = line.split_whitespace().collect();
                    mol = paras[0].to_string();
                    nrexcl = paras[1].parse().unwrap();
                } else if cur_item == "atoms" {
                    // Fuck, I must change the atom name here
                    let mut a = TopolAtom::from(line);
                    a.atom = mol2.atoms[a.nr - 1].atom_name.to_string();
                    atoms.push(a);
                } else if cur_item == "bonds" {
                    bonds.push(TopolBond::from(&atoms, line));
                } else if cur_item == "pairs" {
                    pairs.push(TopolPair::from(&atoms, line));
                } else if cur_item == "constraints" {
                    constraints.push(TopolConstraint::from(&atoms, line));
                } else if cur_item == "angles" {
                    angles.push(TopolAngle::from(&atoms, line));
                } else if cur_item == "dihedrals" {
                    dihedrals.push(TopolDihedral::from(&atoms, line));
                } else if cur_item == "exclusions" {
                    exclusions.push(TopolExclusion::from(&atoms, line));
                }
            }
        }
        println!("Finished reading topology of {}\n", mol);
        Topol {
            atomtypes: attypes, 
            moleculetype: mol,
            nrexcl, atoms, bonds, pairs, constraints, angles, dihedrals, exclusions
        }
    }

    pub fn to_rtp(&mut self, outfile: &str, ff: &str, 
        exclude_n: &Vec<usize>, exclude_c: &Vec<usize>, 
        atom_n: Option<usize>, atom_c: Option<usize>,
        n_name: &Option<String>, c_name: &Option<String>) {
        // 前后残基中的原子名加前缀
        for atom in &mut self.atoms {
            if let Some(atom_n) = atom_n {
                if exclude_n.contains(&atom.nr) {
                    if atom.nr != atom_n {
                        atom.atom = "-".to_string() + &atom.atom;
                    } else {
                        atom.atom = n_name.as_ref().unwrap().to_string();
                    }
                }
            }
            if let Some(atom_c) = atom_c {
                if exclude_c.contains(&atom.nr) {
                    if atom.nr != atom_c {
                        atom.atom = "+".to_string() + &atom.atom;
                    } else {
                        atom.atom = c_name.as_ref().unwrap().to_string();
                    }
                }
            }
        }

        let mut file = fs::File::create(outfile).unwrap();
        
        file.write_all(b"; rtp created by gen-rtp\n").unwrap();
        file.write_all(format!("; converted from top of {}\n\n", self.moleculetype).as_bytes()).unwrap();
        if ff == "amber" {
            file.write_all(b"[ bondedtypes ]\n").unwrap();
            file.write_all(b"; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih\n").unwrap();
            file.write_all(b"     1       1          9          4        1         3      1     0\n\n").unwrap();
        } else if ff == "gromos" {
            file.write_all(b"[ bondedtypes ]\n; bonds  angles  dihedrals  impropers\n").unwrap();
            file.write_all(b"    2       2          1          2\n\n").unwrap();
        } else {
            println!("Error: invalid forcefield, only support amber and gromos.\n");
            return
        }
    
        // 残基名
        file.write_all(format!("[ {} ]\n", self.moleculetype).as_bytes()).unwrap();
    
        // [ atoms ]字段：记录残基中每个原子的名称、类型和电荷、电荷组
        file.write_all(b" [ atoms ]\n").unwrap();
        for atom in &self.atoms {
            if !exclude_n.is_empty() && exclude_n.contains(&atom.nr) {
                file.write_all(format!("; {}\t; previous residue\n", atom.to_rtp()).as_bytes()).unwrap();
            } else if !exclude_c.is_empty() && exclude_c.contains(&atom.nr) {
                file.write_all(format!("; {}\t; next residue\n", atom.to_rtp()).as_bytes()).unwrap();
            } else {
                file.write_all((atom.to_rtp() + "\n").as_bytes()).unwrap();
            }
        }
    
        // [ bonds ]字段：原子间的连接信息
        file.write_all(b" [ bonds ]\n").unwrap();
        self.bonds.retain(|b| !exclude_n.contains(&b.ai.nr) || !exclude_n.contains(&b.aj.nr));
        self.bonds.retain(|b| !exclude_c.contains(&b.ai.nr) || !exclude_c.contains(&b.aj.nr));
        
        match ff {      // amber力场保留-, gromos力场保留+
            "amber" => {
                self.bonds.retain(|b| !exclude_c.contains(&b.ai.nr) && !exclude_c.contains(&b.aj.nr));
                for bond in &mut self.bonds {
                    if let Some(atom_n) = atom_n {
                        if bond.ai.nr == atom_n {
                            bond.ai.atom = n_name.as_ref().unwrap().to_string();
                        } else if bond.aj.nr == atom_n {
                            bond.aj.atom = n_name.as_ref().unwrap().to_string();
                        }
                    }
                }
            },
            "gromos" => {
                self.bonds.retain(|b| !exclude_n.contains(&b.ai.nr) && !exclude_n.contains(&b.aj.nr));
                for bond in &mut self.bonds {
                    if let Some(atom_c) = atom_c {
                        if bond.ai.nr == atom_c {
                            bond.ai.atom = c_name.as_ref().unwrap().to_string();
                        } else if bond.aj.nr == atom_c {
                            bond.aj.atom = c_name.as_ref().unwrap().to_string();
                        }
                    }
                }
            },
            _ => ()
        }
        for bond in &self.bonds {
            if !exclude_n.is_empty() && exclude_n.contains(&bond.ai.nr) && exclude_n.contains(&bond.aj.nr) {
                file.write_all(format!(";- {}\n", bond.to_rtp()).as_bytes()).unwrap();
            } else if exclude_c.is_empty() && (exclude_c.contains(&bond.ai.nr) && exclude_c.contains(&bond.aj.nr)) {
                file.write_all(format!(";+ {}\n", bond.to_rtp()).as_bytes()).unwrap();
            } else {
                file.write_all((bond.to_rtp() + "\n").as_bytes()).unwrap();
            }
        }
    
        // [ angles ]字段：键角信息
        file.write_all(b" [ angles ]\n").unwrap();
        self.angles.retain(|a| !exclude_n.contains(&a.ai.nr) && !exclude_n.contains(&a.aj.nr) && !exclude_n.contains(&a.ak.nr));
        self.angles.retain(|a| !exclude_c.contains(&a.ai.nr) && !exclude_c.contains(&a.aj.nr) && !exclude_c.contains(&a.ak.nr));
        for angle in &self.angles {
            file.write_all((angle.to_rtp() + "\n").as_bytes()).unwrap();
        }
    
        // [ dihedrals ]字段：二面角信息
        file.write_all(b" [ dihedrals ]\n").unwrap();
        self.dihedrals.retain(|d| !exclude_n.contains(&d.ai.nr) && !exclude_n.contains(&d.aj.nr) && !exclude_n.contains(&d.ak.nr) && !exclude_n.contains(&d.al.nr));
        self.dihedrals.retain(|d| !exclude_c.contains(&d.ai.nr) && !exclude_c.contains(&d.aj.nr) && !exclude_c.contains(&d.ak.nr) && !exclude_c.contains(&d.al.nr));
        for dihedral in &self.dihedrals {
            // 理论上2用来描述improper, 但sobtop生成拓扑时采用2描述proper, 这里为特殊应对
            if dihedral.funct == 9 || dihedral.funct == 1 || dihedral.funct == 2 {
                file.write_all((dihedral.to_rtp(ff) + "\n").as_bytes()).unwrap();
            }
        }
    
        // [ impropers ]字段：反常二面角信息
        file.write_all(b" [ impropers ]\n").unwrap();
        for dihedral in &self.dihedrals {
            if dihedral.funct == 4 {
                file.write_all((dihedral.to_rtp(ff) + "\n").as_bytes()).unwrap();
            }
        }
    
        println!("Finished writing rtp file to {}", outfile);
    }
}

impl Display for Topol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut out = format!("; Created by gen-rtp (https://gitee.com/supernova_bingbing/SuperMDA)\n");
        // 输出原子类型
        if !self.atomtypes.is_empty() {
            out.push_str("\n[ atomtypes ]\n; name    at.num        mass       charge    ptype      sigma (nm)      epsilon (kJ/mol)\n");
            for at in &self.atomtypes {
                out.push_str(format!("{}\n", at).as_str());
            }
        }
        // 输出残基名
        out.push_str("\n[ moleculetype ]\n; name   nrexcl\n");
        out.push_str(format!("  {}    {}\n", self.moleculetype, self.nrexcl).as_str());
        // 输出原子
        out.push_str("\n[ atoms ]\n;    nr   type  resnr  resid    atom   cgnr      charge     mass\n");
        if !self.atoms.is_empty() {
            for a in &self.atoms {
                out.push_str(format!("{}\n", a).as_str());
            }
        }
        // 输出键
        out.push_str("\n[ bonds ]\n;    ai     aj    funct           c0           c1\n");
        if !self.bonds.is_empty() {
            for b in &self.bonds {
                out.push_str(format!("{}\n", b).as_str());
            }
        }
        // 输出 pair 编号
        if !self.pairs.is_empty() {
            out.push_str("\n[ pairs ]\n;    ai     aj    funct  ;  all 1-4 pairs but the ones excluded in GROMOS itp\n");
            for p in &self.pairs {
                out.push_str(format!("{}\n", p).as_str());
            }
        }
        // 输出 constraint 编号
        if !self.constraints.is_empty() {
            out.push_str("\n[ constraints ]\n");
            for c in &self.constraints {
                out.push_str(format!("{}\n", c).as_str());
            }
        }
        // 输出键角编号
        out.push_str("\n[ angles ]\n;    ai     aj     ak    funct     angle       fc\n");
        if !self.angles.is_empty() {
            for a in &self.angles {
                out.push_str(format!("{}\n", a).as_str());
            }
        }
        // 输出 proper 二面角编号
        out.push_str("\n[ dihedrals ]; proper\n;    ai     aj     ak     al    funct     phase       kd      pn   ; funct = 9 or 1\n");
        if !self.dihedrals.is_empty() {
            for d in &self.dihedrals {
                if d.funct == 1 || d.funct == 9 {
                    out.push_str(format!("{}\n", d).as_str());
                }
            }
        }
        // # 输出 improper 二面角编号
        out.push_str("\n[ dihedrals ]; improper\n;    ai     aj     ak     al    funct     phase       kd      pn   ; funct = 4\n");
        out.push_str(";    ai     aj     ak     al    funct      ksi         k           ; funct = 2\n");
        if !self.dihedrals.is_empty() {
            for d in &self.dihedrals {
                if d.funct == 4 {
                    out.push_str(format!("{}\n", d).as_str());
                }
            }
        }
        if !self.dihedrals.is_empty() {
            for d in &self.dihedrals {
                if d.funct == 2 {
                    out.push_str(format!("{}\n", d).as_str());
                }
            }
        }
        // 输出 exclusion 编号
        if !self.exclusions.is_empty() {
            out.push_str("\n[ exclusions ]\n");
            for ex in &self.exclusions {
                out.push_str(format!("{}\n", ex).as_str());
            }
        }
        write!(f, "{}", out)
    }
}

impl TopolAtomtype {
    fn from(line: &String) -> TopolAtomtype {
        let paras: Vec<&str> = line.split_whitespace().collect();
        let name = paras[0].to_string();
        // 有时候第二个是原子序号, 所以质量倒数计数
        let mass: f64 = paras[paras.len() - 5].parse().unwrap();
        let charge: f64 = paras[paras.len() - 4].parse().unwrap();
        let ptype = paras[paras.len() - 3].to_string();
        let sigma: f64 = paras[paras.len() - 2].parse().unwrap();
        let epsilon: f64 = paras.last().unwrap().parse().unwrap();
        TopolAtomtype{ name, mass, charge, ptype, sigma, epsilon }
    }
}

impl Display for TopolAtomtype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "  {:8}{:8}{:10.6}{:13.6}{:>5}{:18.6}{:16.6}", self.name, self.name, self.mass, self.charge, self.ptype, self.sigma, self.epsilon)
    }
}

impl TopolAtom {
    fn from(line: &String) -> TopolAtom {
        let paras: Vec<&str> = line.split_whitespace().collect();
        let nr: usize = paras[0].parse().unwrap();
        let _type = paras[1].to_string();
        let resnr: i32 = paras[2].parse().unwrap();
        let resname = paras[3].to_string();
        let atom = paras[4].to_string();
        let cgnr: i32 = paras[5].parse().unwrap();
        let charge: f64 = paras[6].parse().unwrap();
        let mass: Option<f64> = get_param_at(&paras, 7);
        TopolAtom{ nr, _type, resnr, resname, atom, cgnr, charge, mass }
    }

    fn to_rtp(&self) -> String {
        format!("{:>8}{:>6}{:12.6}{:5}", self.atom, self._type, self.charge, self.cgnr)
    }
}

impl Display for TopolAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let out = match self.mass {
            Some(m) => format!("{:7}{:>7}{:7}{:>7}{:>8}{:7}{:12.6}{:9.4}", 
                self.nr, self._type, self.resnr, self.resname, self.atom, self.cgnr, self.charge, m),
            None => format!("{:7}{:>7}{:7}{:>7}{:>8}{:7}{:12.6}", 
                self.nr, self._type, self.resnr, self.resname, self.atom, self.cgnr, self.charge)
        };
        write!(f, "{}", out)
    }
}

impl TopolBond {
    fn from(atoms: &Vec<TopolAtom>, line: &String) -> TopolBond {
        let paras: Vec<&str> = line.split_whitespace().collect();
        let ai: usize = paras[0].parse().unwrap();
        let aj: usize = paras[1].parse().unwrap();
        let funct: i32 = paras[2].parse().unwrap();
        let ai = get_atom_from_nr(atoms, ai).to_owned();
        let aj = get_atom_from_nr(atoms, aj).to_owned();
        let c0: Option<f64> = get_param_at(&paras, 3);
        let c1: Option<f64> = get_param_at(&paras, 4);
        TopolBond{ ai, aj, funct, c0, c1 }
    }

    fn to_rtp(&self) -> String {
        match self.c0 {
            Some(c0) => format!("{:>7}{:>7}{:13.6}{:13.6e}", self.ai.atom, self.aj.atom, c0, self.c1.unwrap()),
            None => format!("{:7}{:7}", self.ai.atom, self.aj.atom)
        }
    }
}

impl Display for TopolBond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let out = match self.c0 {
            Some(c0) => format!("{:7}{:7}{:9}{:13.6}{:13.6e}", self.ai.nr, self.aj.nr, self.funct, c0, self.c1.unwrap()),
            None => format!("{:7}{:7}{:9}", self.ai.nr, self.aj.nr, self.funct)
        };
        write!(f, "{}", out)
    }
}

impl TopolPair {
    fn from(atoms: &Vec<TopolAtom>, line: &String) -> TopolPair {
        let paras: Vec<&str> = line.split_whitespace().collect();
        let ai: usize = paras[0].parse().unwrap();
        let aj: usize = paras[1].parse().unwrap();
        let funct: i32 = paras[2].parse().unwrap();
        let ai = get_atom_from_nr(atoms, ai).to_owned();
        let aj = get_atom_from_nr(atoms, aj).to_owned();
        TopolPair{ ai, aj, funct }
    }
}

impl Display for TopolPair {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:7}{:7}{:9}", self.ai.nr, self.aj.nr, self.funct)
    }
}

impl TopolConstraint {
    fn from(atoms: &Vec<TopolAtom>, line: &String) -> TopolConstraint {
        let paras: Vec<&str> = line.split_whitespace().collect();
        let ai: usize = paras[0].parse().unwrap();
        let aj: usize = paras[1].parse().unwrap();
        let funct: i32 = paras[2].parse().unwrap();
        let ai = get_atom_from_nr(atoms, ai).to_owned();
        let aj = get_atom_from_nr(atoms, aj).to_owned();
        let cs: Option<Vec<f64>> = match paras.get(4..) {
            Some(s) => Some(s.iter().map(|&s| s.parse().unwrap()).collect()),
            None => None
        };
        TopolConstraint{ ai, aj, funct, cs }
    }
}

impl Display for TopolConstraint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut out = format!("{:7}{:7}{:9}", self.ai.nr, self.aj.nr, self.funct);
        if let Some(cs) = &self.cs {
            for c in cs {
                out.push_str(format!("{:13.4e}", c).as_str());
            }
        }
        write!(f, "{}", out)
    }
}

impl TopolAngle {
    fn from(atoms: &Vec<TopolAtom>, line: &String) -> TopolAngle {
        let paras: Vec<&str> = line.split_whitespace().collect();
        let ai: usize = paras[0].parse().unwrap();
        let aj: usize = paras[1].parse().unwrap();
        let ak: usize = paras[2].parse().unwrap();
        let funct: i32 = paras[3].parse().unwrap();
        let ai = get_atom_from_nr(atoms, ai).to_owned();
        let aj = get_atom_from_nr(atoms, aj).to_owned();
        let ak = get_atom_from_nr(atoms, ak).to_owned();
        let c0: Option<f64> = get_param_at(&paras, 4);
        let c1: Option<f64> = get_param_at(&paras, 5);
        TopolAngle{ ai, aj, ak, funct, c0, c1 }
    }

    fn to_rtp(&self) -> String {
        match self.c0 {
            Some(c0) => format!("{:>7}{:>7}{:>7}{:10.2}{:9.2}", self.ai.atom, self.aj.atom, self.ak.atom, c0, self.c1.unwrap()),
            None => format!("{:7}{:7}{:7}", self.ai.atom, self.aj.atom, self.ak.atom)
        }
    }
}

impl Display for TopolAngle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let out = match self.c0 {
            Some(c0) => format!("{:7}{:7}{:7}{:9}{:10.2}{:9.2}", self.ai.nr, self.aj.nr, self.ak.nr, self.funct, c0, self.c1.unwrap()),
            None => format!("{:7}{:7}{:7}{:9}", self.ai.nr, self.aj.nr, self.ak.nr, self.funct)
        };
        write!(f, "{}", out)
    }
}

impl TopolDihedral {
    fn from(atoms: &Vec<TopolAtom>, line: &String) -> TopolDihedral {
        let paras: Vec<&str> = line.split_whitespace().collect();
        let ai: usize = paras[0].parse().unwrap();
        let aj: usize = paras[1].parse().unwrap();
        let ak: usize = paras[2].parse().unwrap();
        let al: usize = paras[3].parse().unwrap();
        let funct: i32 = paras[4].parse().unwrap();
        let ai = get_atom_from_nr(atoms, ai).to_owned();
        let aj = get_atom_from_nr(atoms, aj).to_owned();
        let ak = get_atom_from_nr(atoms, ak).to_owned();
        let al = get_atom_from_nr(atoms, al).to_owned();
        let c0: Option<f64> = get_param_at(&paras, 5);
        let c1: Option<f64> = get_param_at(&paras, 6);
        let c2: Option<f64> = get_param_at(&paras, 7);
        TopolDihedral{ ai, aj, ak, al, funct, c0, c1, c2 }
    }
    fn to_rtp(&self, ff: &str) -> String {
        match self.c0 {
            Some(c0) => {
                match self.funct {
                    // 如果是2且是amber力场, 将funct作为第一个参数, 提示删掉
                    2 => match ff {
                        "amber" => format!("{:>7}{:>7}{:>7}{:>7}{:9}{:10.2}{:9.2}        ; Delete the default funct \"9\" before 2 after pdb2gmx!!!", 
                                    self.ai.atom, self.aj.atom, self.ak.atom, self.al.atom, self.funct, c0, self.c1.unwrap()),
                        "gromos" => format!("{:>7}{:>7}{:>7}{:>7}{:10.2}{:9.2}", 
                                    self.ai.atom, self.aj.atom, self.ak.atom, self.al.atom, c0, self.c1.unwrap()),
                        _ => "".to_string()
                    },
                    _ => format!("{:>7}{:>7}{:>7}{:>7}{:10.2}{:9.2}{:8}", 
                                self.ai.atom, self.aj.atom, self.ak.atom, self.al.atom, c0, self.c1.unwrap(), self.c2.unwrap())
                }
            },
            None => format!("{:>7}{:>7}{:>7}{:>7}", self.ai.atom, self.aj.atom, self.ak.atom, self.al.atom)
        }
    }
}

impl Display for TopolDihedral {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let out = match self.c0 {
            Some(c0) => {
                match self.funct {
                    2 => format!("{:7}{:7}{:7}{:7}{:9}{:10.2}{:9.2}", 
                                self.ai.nr, self.aj.nr, self.ak.nr, self.al.nr, self.funct, c0, self.c1.unwrap()),
                    _ => format!("{:7}{:7}{:7}{:7}{:9}{:10.2}{:9.2}{:8}", 
                                self.ai.nr, self.aj.nr, self.ak.nr, self.al.nr, self.funct, c0, self.c1.unwrap(), self.c2.unwrap())
                }
            },
            None => format!("{:7}{:7}{:7}{:7}{:9}", self.ai.nr, self.aj.nr, self.ak.nr, self.al.nr, self.funct)
        };
        write!(f, "{}", out)
    }
}

impl TopolExclusion {
    fn from(atoms: &Vec<TopolAtom>, line: &String) -> TopolExclusion {
        let paras: Vec<&str> = line.split_whitespace().collect();
        let atnums: Option<Vec<usize>> = match paras.get(..) {
            Some(s) => Some(s.iter().map(|&s| s.parse().unwrap()).collect()),
            None => None
        };
        let ex_atoms = atnums.unwrap().iter().map(|&a| get_atom_from_nr(&atoms, a).to_owned()).collect();
        TopolExclusion{ ex_atoms }
    }
}

impl Display for TopolExclusion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut out = "".to_string();
        for a in &self.ex_atoms {
            out.push_str(format!("{:8}", a.nr).as_str());
        }
        write!(f, "{}", out)
    }
}

fn get_param_at<T: FromStr>(paras: &Vec<&str>, id: usize) -> Option<T> where <T as FromStr>::Err: Debug {
    let p: Option<T> = match paras.get(id) {
        Some(&s) => Some(s.parse().unwrap()),
        None => None
    };
    return p
}

// 处理range泛型问题
// fn get_params_at<T: FromStr, U>(paras: &Vec<&str>, range: U) -> Option<Vec<T>> 
// where <T as FromStr>::Err: Debug, U: SliceIndex<Vec<&str>> + std::slice::SliceIndex<[&str]> {
//     let cs: Option<Vec<T>> = match paras.get(range) {
//         Some(s) => Some(s.iter().map(|&s| s.parse().unwrap()).collect()),
//         None => None
//     };
//     return cs
// }
