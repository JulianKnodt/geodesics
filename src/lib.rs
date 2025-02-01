#![feature(cmp_minmax)]

use std::cmp::minmax;
use std::collections::HashMap;
use std::collections::hash_map::Entry;

use clap::Parser;
use pars3d::FaceKind;

pub type F = f32;

#[derive(Debug, Parser)]
pub struct Args {
    /// Path to input mesh
    #[arg(short, long, required = true)]
    pub input: String,

    /// Index of source vertex.
    #[arg(short, long, required = true)]
    pub src_idx: usize,
}

pub fn geodesics(src_idx: usize, vs: &[[F; 3]], fs: &[FaceKind]) {
    let src = vs[src_idx];

    let mut edge_len_sum = 0.;
    let mut num_edges = 0;

    let mut edge_face_adj: HashMap<[usize; 2], Vec<usize>> = HashMap::new();
    let mut vert_adjs: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut edge_lens: HashMap<[usize; 2], F> = HashMap::new();
    let mut vert_face_adj: HashMap<usize, Vec<usize>> = HashMap::new();

    for (fi, f) in fs.iter().enumerate() {
        let f = f.as_slice();
        for &vi in f {
            vert_face_adj.entry(vi).or_default().push(fi);
        }
        for i in 0..f.len() {
            let n = (i + 1) % f.len();
            let edge_len = dist(vs[f[i]], vs[f[n]]);

            let e = minmax(i, n);
            match edge_face_adj.entry(e) {
                Entry::Occupied(mut o) => {
                    o.get_mut().push(fi);
                }
                Entry::Vacant(v) => {
                    edge_len_sum += edge_len;
                    num_edges += 1;
                    edge_lens.insert(e, edge_len);
                    v.insert(vec![fi]);

                    vert_adjs.entry(i).or_default().push(n);
                    vert_adjs.entry(n).or_default().push(i);
                }
            }
        }
    }

    let avg_edge_len = edge_len_sum / num_edges as F;
    let inv_avg_edge_len = avg_edge_len.recip();

    for el in edge_lens.values_mut() {
        *el *= inv_avg_edge_len;
    }

    let mut iter = 0;

    let mut face_to_extra_dist = vec![0.; fs.len()];
    let mut face_to_center_dist = vec![F::INFINITY; fs.len()];
    let mut edge_to_vert_dist = HashMap::new();

    use std::collections::VecDeque;
    let mut queue_from = VecDeque::new();
    let mut queue_to = VecDeque::new();

    for &adj_face in &vert_face_adj[&src_idx] {
        let centroid = fs[adj_face]
            .as_slice()
            .into_iter()
            .fold([0.; 3], |acc, &vi| add(acc, vs[vi]));

        face_to_extra_dist[adj_face] = 0.;
        face_to_center_dist[adj_face] = dist(centroid, src);

        for [e0, e1] in fs[adj_face].edges() {
          edge_to_vert_dist.insert((minmax(e0, e1), adj_face), dist(src, vs[e1]));
        }

        queue_from.extend(
            fs[adj_face]
                .edges()
                .filter(|e| !e.contains(&src_idx))
                .map(|e| (e, adj_face, false)),
        );
    }

    while !queue_from.is_empty() {
        iter += 1;

        queue_to.clear();

        while let Some(([e0, e1] @ og_ord_edge, src_f, is_behind)) = queue_from.pop_front() {
            assert!(face_to_center_dist[src_f].is_finite());
            assert!(face_to_extra_dist[src_f] >= 0.);

            let edge = minmax(e0, e1);
            let e_len = edge_lens[&edge];

            // for each adjacent face to this edge propagate distances
            for &tgt_f in &edge_face_adj[&edge] {
                if tgt_f == src_f {
                    continue;
                }

                let ref_d_t = face_to_center_dist[tgt_f];
                let prev_sigma_t = face_to_extra_dist[src_f];

                let mut els = fs[tgt_f]
                    .edges()
                    .map(|[a, b]| minmax(a, b))
                    .filter(|&e| e != edge)
                    .map(|e| edge_lens[&e])
                    .map(sqr);
                let e2 = els.next().unwrap();

                let d1_sqr = edge_to_vertex_dist[&(edge, src_f)];
                let d2_sqr =

                let px = (sqr(e_len) + (e2 - els.sum::<F>())) / (e_len + e_len);
                let py = (e2 - sqr(px)).max(0.).sqrt();

                let cx = (px + e_len) / 3.;
                let cy = py / 3.;

                let sigma_t = prev_sigma_t;
            }
        }

        std::mem::swap(&mut queue_from, &mut queue_to);
    }
}

pub fn sqr(x: F) -> F {
    x * x
}

pub fn add<const N: usize>(a: [F; N], b: [F; N]) -> [F; N] {
    std::array::from_fn(|i| a[i] + b[i])
}

pub fn sub<const N: usize>(a: [F; N], b: [F; N]) -> [F; N] {
    std::array::from_fn(|i| a[i] - b[i])
}

pub fn length_sq<const N: usize>(v: [F; N]) -> F {
    v.into_iter().map(|a| a * a).sum::<F>()
}

pub fn length<const N: usize>(v: [F; N]) -> F {
    length_sq(v).sqrt()
}
pub fn dist<const N: usize>(a: [F; N], b: [F; N]) -> F {
    length(sub(a, b))
}
