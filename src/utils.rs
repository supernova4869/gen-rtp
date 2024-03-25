use std::io;
use regex::Regex;

pub fn read_file() -> String {
    let mut inp: String = String::new();
    io::stdin().read_line(&mut inp).expect("Failed to read line");
    let inp = inp.trim();
    let inp: String = match inp.starts_with("\"") {
        true => inp[1..inp.len() - 1].to_string(),
        false => inp.to_string()
    };
    let re = Regex::new(r"\\").unwrap();
    re.replace_all(&inp, "/").to_string()
}
