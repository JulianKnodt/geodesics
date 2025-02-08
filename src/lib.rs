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

    /// Path to colored PLY output.
    #[arg(short, long, default_value = "out.ply")]
    pub output: String,
}

/*
enum EdgeFaceAdj {
  Boundary,
  Manifold(usize),
  NonManifold(Vec<usize>)
}

impl EdgeFaceAdj {
  pub fn as_slice(&self) -> &[usize] {
    match self {
      Self::Boundary => &[],
      Self::Manifold(v) => std::slice::from_ref(v),
      Self::NonManifold(a) => a.as_slice(),
    }
  }
}
*/

pub fn geodesics(src_idx: usize, vs: &[[F; 3]], fs: &[FaceKind]) -> Vec<F> {
    let src = vs[src_idx];

    let start = std::time::Instant::now();

    let mut edge_len_sum = 0.;
    let mut num_edges = 0;

    let mut edge_face_adj: HashMap<[usize; 2], Vec<usize>> = HashMap::new();
    let mut vert_adjs: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut face_edge_lens: Vec<[F; 3]> = vec![[0.; 3]; fs.len()];

    let mut vert_face_adj: Vec<Vec<usize>> = vec![vec![]; vs.len()];

    for (fi, f) in fs.iter().enumerate() {
        let f = f.as_slice();
        for &vi in f {
            vert_face_adj[vi].push(fi);
        }
        for i in 0..f.len() {
            let n = (i + 1) % f.len();
            let edge_len = dist(vs[f[i]], vs[f[n]]);
            face_edge_lens[fi][i] = edge_len;

            let e = minmax(f[i], f[n]);
            match edge_face_adj.entry(e) {
                Entry::Occupied(mut o) => {
                    o.get_mut().push(fi);
                }
                Entry::Vacant(v) => {
                    edge_len_sum += edge_len;
                    num_edges += 1;
                    v.insert(vec![fi]);

                    vert_adjs.entry(i).or_default().push(n);
                    vert_adjs.entry(n).or_default().push(i);
                }
            }
        }
    }

    let avg_edge_len = edge_len_sum / num_edges as F;
    let inv_avg_edge_len = avg_edge_len.recip();

    for els in face_edge_lens.iter_mut() {
        for el in els.iter_mut() {
            *el *= inv_avg_edge_len;
        }
    }

    let mut iters = 0;
    let mut updates = 0;
    let mut expansions = 0;

    let mut face_to_extra_dist = vec![0.; fs.len()];
    let mut face_to_center_dist = vec![F::INFINITY; fs.len()];
    let mut vert_to_src_dist: Vec<[F; 3]> = vec![[F::INFINITY; 3]; fs.len()];

    use std::collections::VecDeque;
    let mut queue_from = VecDeque::new();
    let mut queue_to = VecDeque::new();

    for &adj_face in &vert_face_adj[src_idx] {
        let fk = &fs[adj_face];
        let f = fk.as_slice();
        assert!(f.contains(&src_idx));
        let centroid = f.iter().fold([0.; 3], |acc, &vi| add(acc, vs[vi]));

        face_to_center_dist[adj_face] = dist(centroid, src);
        let si = f.iter().position(|&i| i == src_idx).unwrap();
        vert_to_src_dist[adj_face][si] = 0.;

        for (vi, &v) in f.iter().enumerate() {
            if v == src_idx {
                continue;
            }
            let dist = face_edge_lens[adj_face][vi];
            vert_to_src_dist[adj_face][vi] = sqr(dist);
        }

        queue_from.extend(
            fk.edge_idxs()
                .filter(|&[e0, e1]| f[e0] != src_idx && f[e1] != src_idx)
                .map(|e| (e, adj_face, false)),
        );
    }

    while !queue_from.is_empty() {
        iters += 1;

        queue_to.clear();

        while let Some(([eii0, eii1], src_f, is_behind)) = queue_from.pop_front() {
            debug_assert!(face_to_center_dist[src_f].is_finite());
            debug_assert!(face_to_extra_dist[src_f] >= 0.);
            let ei0 = fs[src_f].as_slice()[eii0];
            let ei1 = fs[src_f].as_slice()[eii1];
            debug_assert_eq!(fs[src_f].next(ei0), ei1);

            expansions += 1;

            let e1 = face_edge_lens[src_f][eii0];

            let d1_sqr = vert_to_src_dist[src_f][eii1];
            debug_assert!(d1_sqr.is_finite());
            let d2_sqr = vert_to_src_dist[src_f][eii0];
            debug_assert!(d2_sqr.is_finite());

            let sx = (sqr(e1) + (d1_sqr - d2_sqr)) / (e1 + e1);
            let sy_neg = (d1_sqr - sqr(sx)).max(0.).sqrt();

            let prev_sigma_t = face_to_extra_dist[src_f];

            // for each adjacent face to this edge propagate distances
            for &tgt_f in &edge_face_adj[&minmax(ei0, ei1)] {
                if tgt_f == src_f {
                    continue;
                }
                assert!(tgt_f < face_edge_lens.len());

                let Some([v0, v1, v2]) = fs[tgt_f].as_tri() else {
                    continue;
                };
                let opp_i = other_i([v0, v1, v2], ei0, ei1);

                let e2 = face_edge_lens[tgt_f][opp_i];
                let prev_i = opp_i.checked_sub(1).unwrap_or(2);
                let e3 = face_edge_lens[tgt_f][prev_i];
                let next_i = (opp_i + 1) % 3;

                let px = (sqr(e1) + sqr(e2) - sqr(e3)) / (e1 + e1);
                let py = (sqr(e2) - sqr(px)).max(0.).sqrt();

                let cx = (px + e1) / 3.;
                let cy = py / 3.;

                let d_s_1 = (sqr(sx) + sqr(sy_neg)).sqrt();
                let d_s_2 = (sqr(sx - e1) + sqr(sy_neg)).sqrt();

                let d_c_1 = (sqr(cx) + sqr(cy)).sqrt();
                let d_c_2 = (sqr(cx - e1) + sqr(cy)).sqrt();

                let mut h2_behind = false;
                let mut h3_behind = false;
                let (d_a, d_b, d_c, sigma_t, d_t) = if is_behind {
                    let dis1 = d_c_1 + d_s_1;
                    let dis2 = d_c_2 + d_s_2;

                    if dis1 < dis2 {
                        let sigma_t = prev_sigma_t + d_s_1;
                        (0., sqr(e1), sqr(e2), sigma_t, sigma_t + d_c_1)
                    } else {
                        let sigma_t = prev_sigma_t + d_s_2;
                        (sqr(e1), 0., sqr(e3), sigma_t, sigma_t + d_c_2)
                    }
                } else {
                    let [bend_left, bend_right] = data_driven_bending_heuristic(
                        e1,
                        &[e1, e2, e3],
                        px,
                        py,
                        cx,
                        cy,
                        sx,
                        sy_neg,
                    );

                    if bend_left {
                        let sigma_t = prev_sigma_t + d_s_1;
                        (0., sqr(e1), sqr(e2), sigma_t, sigma_t + d_c_1)
                    } else if bend_right {
                        let sigma_t = prev_sigma_t + d_s_2;
                        (sqr(e1), 0., sqr(e3), sigma_t, sigma_t + d_c_2)
                    } else {
                        h2_behind = py * sx + px * sy_neg < 0.;
                        h3_behind = py * (e1 - sx) - (px - e1) * sy_neg < 0.;
                        let d_t = prev_sigma_t + (sqr(cx - sx) + sqr(cy + sy_neg)).sqrt();
                        (
                            d1_sqr,
                            d2_sqr,
                            sqr(px - sx) + sqr(py + sy_neg),
                            prev_sigma_t,
                            d_t,
                        )
                    }
                };

                if d_t < face_to_center_dist[tgt_f] {
                    updates += 1;

                    face_to_center_dist[tgt_f] = d_t;
                    face_to_extra_dist[tgt_f] = sigma_t;

                    vert_to_src_dist[tgt_f][opp_i] = d_c;
                    vert_to_src_dist[tgt_f][next_i] = d_a;
                    vert_to_src_dist[tgt_f][prev_i] = d_b;

                    let next_queue = if d_t <= iters as F {
                        &mut queue_from
                    } else {
                        &mut queue_to
                    };
                    next_queue.push_back(([opp_i, next_i], tgt_f, h3_behind));
                    next_queue.push_back(([prev_i, opp_i], tgt_f, h2_behind));
                }
            }
        }

        std::mem::swap(&mut queue_from, &mut queue_to);
    }
    println!("{:?}", start.elapsed());

    println!("Iters = {iters}");
    println!("Expansions = {expansions}");
    println!("Updates = {updates}");

    let mut out_dists = vec![F::INFINITY; vs.len()];
    // for each edge, check face_to_extra_dist + vert_to_src_dist
    // also need to multiply by avg edge length
    for (fi, dists) in vert_to_src_dist.iter().enumerate() {
        for (vi, &dist) in dists.iter().enumerate() {
          if !dist.is_finite() {
              continue;
          }
          let sigma = face_to_extra_dist[fi];
          let vi = fs[fi].as_slice()[vi];
          out_dists[vi] = out_dists[vi].min(sigma + dist.sqrt());
        }
    }

    for d in out_dists.iter_mut() {
        *d *= avg_edge_len;
    }

    out_dists
}

pub fn data_driven_bending_heuristic(
    e_len: F,
    es: &[F],

    px: F,
    py: F,

    cx: F,
    cy: F,

    sx: F,
    sy_neg: F,
) -> [bool; 2] {
    const THRESH_C: F = 5.1424;
    const THRESH_G: F = 4.20638;
    const THRESH_H: F = 0.504201;
    const THRESH_HG: F = 2.84918;
    const LAMBDA: [F; 16] = [
        0.320991, 0.446887, 0.595879, 0.270094, 0.236679, 0.159685, 0.0872932, 0.434132, 1.0,
        0.726262, 0.0635997, 0.0515979, 0.56903, 0.0447586, 0.0612103, 0.718198,
    ];

    let [max_e, min_e] = es
        .iter()
        .fold([F::INFINITY, F::NEG_INFINITY], |[l, h], &n| {
            [l.min(n), h.max(n)]
        });

    let b0 = max_e > THRESH_C * e_len;
    let b1 = max_e > THRESH_G * min_e;

    let b2 = py < THRESH_H * e_len;
    let b3 = py < THRESH_HG * max_e;

    let idx = b0 as usize + (b1 as usize * 2) + (b2 as usize * 4) + (b3 as usize * 8);
    let l = unsafe { LAMBDA.get_unchecked(idx) };

    let qx = px * (1. - l) + cx * l;
    let qy = py * (1. - l) + cy * l;

    let ttx = qx * sy_neg + sx * qy;
    [ttx < 0., ttx > e_len * (qy + sy_neg)]
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
    dist_sq(a, b).sqrt()
}

pub fn dist_sq<const N: usize>(a: [F; N], b: [F; N]) -> F {
    length_sq(sub(a, b))
}

pub fn other(v: [usize; 3], e0: usize, e1: usize) -> (usize, [usize; 2]) {
    match v {
        [a, b, c] | [c, a, b] | [b, c, a] if b == e0 && c == e1 => (a, [b, c]),
        [a, c, b] | [b, a, c] | [c, b, a] if b == e0 && c == e1 => (a, [c, b]),
        _ => unreachable!(),
    }
}

pub fn other_i(v: [usize; 3], e0: usize, e1: usize) -> usize {
    match v {
        [_, b, c] | [_, c, b] if b == e0 && c == e1 => 0,
        [b, _, c] | [c, _, b] if b == e0 && c == e1 => 1,
        [b, c, _] | [c, b, _] if b == e0 && c == e1 => 2,
        _ => unreachable!(),
    }
}
