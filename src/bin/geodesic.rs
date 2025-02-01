use clap::Parser;
use geodesics::{Args, geodesics};

pub fn main() {
    let args = Args::parse();
    let scene = pars3d::load(&args.input).expect("Failed to load mesh");

    let mesh = scene.into_flattened_mesh();
    assert!(args.src_idx < mesh.v.len(), "Source vertex out of bounds");
    geodesics(args.src_idx, &mesh.v, &mesh.f);
}
