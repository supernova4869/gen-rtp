mod mol2;
mod itp;
mod hdb;
mod utils;

use std::path::Path;
use mol2::MOL2;
use hdb::get_adj_h_id;
use itp::Topol;
use std::io;

fn main() {
    // 读取mol2
    println!("Input the mol2 file name:");
    let mol2_file = utils::read_file();
    println!("Reading mol2 file: {}", mol2_file);
    
    let mol2 = &mut MOL2::from(mol2_file.as_str());
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
    let mol2_file = Path::new(mol2_file.as_str());
    let mol2_name = mol2_file.file_name().unwrap().to_str().unwrap();
    let parent_path = mol2_file.parent().unwrap().join("new".to_string() + mol2_name);
    let parent_path = parent_path.to_str().unwrap();
    mol2.output(parent_path);

    // 读取itp, 选择性删除连接原子成键信息
    println!("Input the itp file name:");
    let itp_file = utils::read_file();
    let mut itp = Topol::from(itp_file.as_str());
    // 输出rtp, 特殊处理2号规则
    let itp_file = Path::new(Path::new(itp_file.as_str()));
    let parent_path = itp_file.parent().unwrap();
    let stem = itp_file.file_stem().unwrap();
    let rtp_name = stem.to_str().unwrap().to_string() + ".rtp";
    let new_path = parent_path.join(rtp_name);
    let new_path = new_path.to_str().unwrap();
    itp.to_rtp(new_path, "amber", vec![], vec![59, 60, 61, 62, 63, 64, 65, 66], 0, 0);
    // 输出hdb, 根据H类型
    let hdb_name = stem.to_str().unwrap().to_string() + ".hdb";
    let new_path = parent_path.join(hdb_name);
    let new_path = new_path.to_str().unwrap();
    mol2.top2hdb(new_path, vec![], vec![59, 60, 61, 62, 63, 64, 65, 66]);

    println!("Press any key to exit");
    io::stdin().read_line(&mut String::new()).expect("Failed to read line");
}
