use std::{fs::File, io::Write};


#[no_mangle]
pub extern "C" fn hello_from_rust() {
    println!("Hello from rust");
    let mut path = std::env::current_dir().unwrap();
    path.push("test");
    path.set_extension("dat");

    let buf = format!("LULE");

    let mut file = File::create(&path).expect("LULE1");
    file.write_all(buf.as_bytes()).expect("LULE2");
}