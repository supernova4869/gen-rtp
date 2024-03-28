mod mol2;
mod itp;
mod hdb;
mod utils;

use mol2::MOL2;
use hdb::get_adj_h_id;
use itp::Topol;
use std::io;
use std::env;
use std::process::exit;

fn main() {
    // 读取mol2
    println!(" GEN-RTP: An `rtp` and `hdb` file generator to be used");
    println!(" together with the Sobtop program.");
    println!(" Website: https://github.com/supernovaZhangJiaXing/gen-rtp");
    println!(" Developed by Jiaxing Zhang, at Tianjin University");
    println!(" Contact me: zhangjiaxing7137@tju.edu.cn");
    println!(" Version 0.1.3, 2024-Mar-28");
    println!();
    println!("Input path of `mol2` file, e.g. D:\\Conan\\Haibara_Ai.mol2");
    println!("(Hint: You can directly load it by `gen-rtp Miyano_Shiho.mol2)`");
    let args: Vec<String> = env::args().collect();
    let mol2_file = match args.len() {
        1 => utils::read_file(),
        2 => args[1].to_owned(),
        _ => {
            println!("Too more arguments.");
            io::stdin().read_line(&mut String::new()).unwrap();
            exit(0)
        }
    };
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
    // 输入排除列表
    let (prev_atoms, 
        next_atoms, 
        prev_con_atom, 
        next_con_atom,
        prev_atom_name,
        next_atom_name) = utils::get_exclude_atoms();
    
    // 输出mol2
    let mol2_stem = utils::get_stemname(&mol2_file);
    let parent_path = utils::get_parent_path(&mol2_file);
    let out = parent_path.join(mol2_stem.to_string() + ".mol2");
    mol2.output(out.as_os_str().to_str().unwrap());

    // 读取itp, 更新H原子名, 选择性删除连接原子成键信息
    let itp_file = parent_path.join(mol2_stem.to_string() + ".itp");
    println!("Input the itp file name (default: {}):", itp_file.to_str().unwrap());
    let inp = utils::get_input(itp_file.to_str().unwrap().to_string());
    let itp_file = match inp.is_empty() {
        true => utils::read_file(),
        false => inp
    };
    let mut itp = Topol::from(itp_file.as_str(), mol2);
    // 输出rtp, 特殊处理2号规则
    let itp_stem = utils::get_stemname(&itp_file);
    let parent_path = utils::get_parent_path(&itp_file);
    let rtp_name = itp_stem.to_string() + ".rtp";
    let out = &parent_path.join(rtp_name);
    let out = out.as_os_str().to_str().unwrap();
    itp.to_rtp(out, 
        "amber", 
        &prev_atoms, 
        &next_atoms, 
        prev_con_atom, 
        next_con_atom,
        &prev_atom_name, 
        &next_atom_name);
    // 输出hdb, 根据H类型
    let hdb_name = itp_stem + ".hdb";
    let out = parent_path.join(hdb_name);
    let out = out.as_os_str().to_str().unwrap();
    mol2.to_hdb(out,
        &prev_atoms, 
        &next_atoms, 
        prev_con_atom, 
        next_con_atom,
        &prev_atom_name, 
        &next_atom_name);

    println!("Press any key to exit");
    io::stdin().read_line(&mut String::new()).expect("Failed to read line");
}
