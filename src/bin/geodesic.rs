use clap::Parser;
use geodesics::{Args, geodesics};

use pars3d::coloring::magma;
use pars3d::visualization::vertex_scalar_coloring;

pub fn main() {
    let args = Args::parse();
    let scene = pars3d::load(&args.input).expect("Failed to load mesh");

    let mut mesh = scene.into_flattened_mesh();
    mesh.triangulate();
    assert!(args.src_idx < mesh.v.len(), "Source vertex out of bounds");

    let out_dists = geodesics(args.src_idx, &mesh.v, &mesh.f);
    assert_eq!(out_dists.len(), mesh.v.len());

    let out_colors = vertex_scalar_coloring(&out_dists, magma, 0.1, [0., 0., 1.]);
    assert_eq!(out_colors.len(), mesh.v.len());
    mesh.vert_colors = out_colors;
    let out_scene = mesh.into_scene();

    pars3d::save(&args.output, &out_scene).expect("Failed to save output mesh");
}
