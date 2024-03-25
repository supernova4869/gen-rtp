mod mol2;
mod htype;
mod itp;
mod hdb;

use mol2::MOL2;
use htype::get_adj_h_id;
use itp::Topol;

fn main() {
    // 读取mol2
    let mol2 = &mut MOL2::from("src/CTP.mol2");
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
                    _ => format!("{}{}", h_basename, i + 1)
                };
            }
        }
    }
    
    // 输出mol2
    mol2.output("new.mol2");

    // 读取itp, 选择性删除连接原子成键信息
    let mut itp = Topol::from("src/CTP.itp");
    // 输出rtp, 特殊处理2号规则
    itp.to_rtp("src/CTP.rtp", "amber", vec![], vec![59, 60, 61, 62, 63, 64, 65, 66], 0, 0);
    // 输出hdb, 根据H类型
    mol2.top2hdb("src/CTP.hdb", vec![], vec![59, 60, 61, 62, 63, 64, 65, 66]);
}
