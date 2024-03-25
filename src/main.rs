mod mol2;
mod htype;
mod itp;
mod hdb;

use std::io;
use std::path::Path;
use mol2::MOL2;
use htype::get_adj_h_id;
use itp::Topol;
use regex::Regex;

fn main() {
    // 读取mol2
    println!("Input the mol2 file name:");
    let mut inp: String = String::new();
    io::stdin().read_line(&mut inp).expect("Failed to read line");
    let inp = inp.trim();
    let inp: String = match inp.starts_with("\"") {
        true => inp[1..inp.len() - 1].to_string(),
        false => inp.to_string()
    };
    let re = Regex::new(r"\\").unwrap();
    let inp = re.replace_all(&inp, "/").to_string();
    println!("Reading mol2 file: {}", inp);
    
    let mol2 = &mut MOL2::from(inp.as_str());
    // 修改H命名
    let mol2_bak = &mol2.clone();
    for a in &mol2_bak.atoms {
        if a.element.ne("H") {
            let adj_h = get_adj_h_id(mol2, a.atom_id);
            for (i, &h) in adj_h.iter().enumerate() {
                // 根据相连H数量修改H名字
                let h_basename = mol2_bak.get_hbasename(&mol2_bak.atoms[h as usize - 1]);
                mol2.atoms[h as usize - 1].atom_name = match adj_h.len() {
                    1 => h_basename,
                    _ => h_basename + (i + 1).to_string().as_str()
                };
            }
        }
    }
    
    // 输出mol2
    let parent_path = Path::new(Path::new(inp.as_str())).parent().unwrap().join("new.mol2");
    let parent_path = parent_path.to_str().unwrap();
    mol2.output(parent_path);

    // 读取itp, 选择性删除连接原子成键信息
    let mut itp = Topol::from("src/CTP.itp");
    // 输出rtp, 特殊处理2号规则
    itp.to_rtp("src/CTP.rtp", "amber", vec![], vec![59, 60, 61, 62, 63, 64, 65, 66], 0, 0);
    // 输出hdb, 根据H类型
    mol2.top2hdb("src/CTP.hdb", vec![], vec![59, 60, 61, 62, 63, 64, 65, 66]);
}
