use std::env;
use std::path::{Path, PathBuf};

fn main() {
    // Target directory for downloaded files
    let mut out_dir = PathBuf::from("target")
        .canonicalize()
        .expect("Cannot canonicalize path");
    out_dir = out_dir.join("aacgm_v2");
    if !Path::new(&out_dir).is_dir() {
        std::fs::create_dir_all(out_dir.clone()).expect("Cannot create target directory");
    }

    // Get the aacgm_v2 C code files
    let code_dir = "c_aacgm_v2.7";
    let code_path = out_dir.join(code_dir);
    let tar_path = out_dir.join(format!("{code_dir}.tar"));
    let tar_file = match std::fs::File::open(&tar_path) {
        Ok(f) => f,
        Err(_) => {
            // Download the aacgm_v2 code
            let code_url = format!("https://superdarn.thayer.dartmouth.edu/aacgm/{code_dir}.tar");
            let mut response =
                reqwest::blocking::get(code_url).expect("Unable to get library tarball");
            // Put the response (the file itself) into archive_file
            let mut archive_file =
                std::fs::File::create(&tar_path).expect("Unable to open library tarball");
            response
                .copy_to(&mut archive_file)
                .expect("Could not copy code tarball");
            std::fs::File::open(&tar_path).expect("Unable to open newly-written code tarball")
        }
    };
    // Unpack the tarball
    tar::Archive::new(tar_file)
        .unpack(&code_path)
        .expect("Unable to unpack library tarball");

    // Get the coefficient files
    let coeffs_dir = "aacgm_coeffs-14";
    let coeffs_path = out_dir.join(coeffs_dir);
    let tar_path = out_dir.join(format!("{coeffs_dir}.tar"));
    let tar_file = match std::fs::File::open(&tar_path) {
        Ok(f) => f,
        Err(_) => {
            // Download the aacgm_v2 coefficients
            let coeffs_url =
                format!("https://superdarn.thayer.dartmouth.edu/aacgm/{coeffs_dir}.tar");
            let mut response =
                reqwest::blocking::get(coeffs_url).expect("Unable to get coefficients tarball");
            // Put the response (the file itself) into archive_file
            let mut archive_file =
                std::fs::File::create(&tar_path).expect("Unable to open coefficients tarball");
            response
                .copy_to(&mut archive_file)
                .expect("Could not copy coefficients tarball");
            std::fs::File::open(&tar_path)
                .expect("Unable to open newly-written coefficients tarball")
        }
    };
    // Unpack the tarball
    tar::Archive::new(tar_file)
        .unpack(&coeffs_path)
        .expect("Unable to unpack coefficients tarball");

    // export necessary environment variables to any crates which depend on this one
    println!("cargo:coeffs_dat_prefix={}", coeffs_path.join("aacgm_coeffs-14-").display());
    println!("cargo:igrf_coeffs={}", code_path.join("magmodel_1590-2025.txt").display());

    // This is the path to the C header file
    let header_path = code_path.join("aacgmlib_v2.h");
    let header_path_str = header_path.to_str().expect("Path is not a valid string");

    // Tell cargo to look for shared libraries in the specified directory
    println!("cargo:rustc-link-search={}", code_path.display());

    // Tell cargo to tell rustc to link the aacgm shared library
    println!("cargo:rustc-link-lib=aacgm_v2");

    // Tell cargo to invalidate the built crate whenever the header file changes
    println!("cargo:rerun-if-changed={}", header_path_str);

    // Run `clang` to compile the `aacgmlib_v2.c` file into a `aacgmlib_v2.o` object file.
    // Unwrap if not possible to spawn the process.
    let pretty_code_path = code_path.clone().into_os_string().into_string().expect("cannot convert code_path to String");
    println!("clang -c -o {pretty_code_path}/aacgmlib_v2.o {pretty_code_path}/aacgmlib_v2.c");
    println!("clang -c -o {pretty_code_path}/igrflib.o {pretty_code_path}/igrflib.c");
    println!("clang -c -o {pretty_code_path}/astalglib.o {pretty_code_path}/astalglib.c");
    println!("clang -c -o {pretty_code_path}/mlt_v2.o {pretty_code_path}/mlt_v2.c");
    println!("clang -c -o {pretty_code_path}/rtime.o {pretty_code_path}/rtime.c");

    if std::fs::File::open(code_path.join("aacgmlib_v2.c")).is_err() {
        panic!("C code file missing!")
    }

    // Compile the aacgmlib_v2 library to an archive file `libaacgmlib_v2.a`
    cc::Build::new()
        .file(code_path.join("aacgmlib_v2.c"))
        .file(code_path.join("astalglib.c"))
        .file(code_path.join("igrflib.c"))
        .file(code_path.join("mlt_v2.c"))
        .file(code_path.join("rtime.c"))
        .warnings(false)
        .compile("aacgm_v2");

    let bindings = bindgen::Builder::default()
        // Throw in stdio.h first so all types are available when generating aacgmlib_v2.h header
        .header_contents("_stdio.h", "#include <stdio.h>\n")
        // The input header to generate bindings for
        .header(format!("{pretty_code_path}/aacgmlib_v2.h"))
        // Tell cargo to invalidate the built crate whenever any of the included header files change
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings
        .generate()
        // Unwrap the result and panic on failure
        .expect("Unable to generate bindings");

    // Write the bindings to $OUT_DIR/bindings.rs
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
