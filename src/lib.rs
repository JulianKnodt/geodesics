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

pub fn geodesics(src_idx: usize, vs: &[[F; 3]], fs: &[FaceKind]) -> Vec<F> {
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

            let e = minmax(f[i], f[n]);
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
    let mut updates = 0;

    let mut face_to_extra_dist = vec![0.; fs.len()];
    let mut face_to_center_dist = vec![F::INFINITY; fs.len()];
    let mut vert_to_src_dist = vec![F::NEG_INFINITY; vs.len()];
    vert_to_src_dist[src_idx] = 0.;

    use std::collections::VecDeque;
    let mut queue_from = VecDeque::new();
    let mut queue_to = VecDeque::new();

    for &adj_face in &vert_face_adj[&src_idx] {
        let f = fs[adj_face].as_slice();
        let centroid = f.into_iter().fold([0.; 3], |acc, &vi| add(acc, vs[vi]));

        face_to_extra_dist[adj_face] = 0.;
        face_to_center_dist[adj_face] = dist(centroid, src);

        for &v in f {
            if v == src_idx || vert_to_src_dist[v].is_finite() {
                continue;
            }
            vert_to_src_dist[v] = dist_sq(src, vs[v]);
        }

        queue_from.extend(
            fs[adj_face]
                .edges()
                .filter(|e| !e.contains(&src_idx))
                .map(|e| (minmax(e[0], e[1]), adj_face, false)),
        );
    }

    while !queue_from.is_empty() {
        iter += 1;

        queue_to.clear();

        while let Some((edge @ [ei0, ei1], src_f, is_behind)) = queue_from.pop_front() {
            assert!(face_to_center_dist[src_f].is_finite());
            assert!(face_to_extra_dist[src_f] >= 0.);
            assert!(ei0 < ei1);

            let e1 = edge_lens[&edge];

            let d1_sqr = vert_to_src_dist[ei0];
            let d2_sqr = vert_to_src_dist[ei1];

            let sx = (sqr(e1) + (d1_sqr - d2_sqr)) / (e1 + e1);
            let sy_neg = (d1_sqr - sqr(sx)).max(0.).sqrt();

            let prev_sigma_t = face_to_extra_dist[src_f];

            // for each adjacent face to this edge propagate distances
            for &tgt_f in &edge_face_adj[&edge] {
                if tgt_f == src_f {
                    continue;
                }

                let [v0, v1, v2] = fs[tgt_f].as_tri().expect("Only tris supported currently");
                let opp = other([v0, v1, v2], ei0, ei1);

                let eo_1 = minmax(ei0, opp);
                let eo_2 = minmax(ei1, opp);

                let e2 = edge_lens[&eo_1];
                let e3 = edge_lens[&eo_2];

                let px = (sqr(e1) + (sqr(e2) - sqr(e3))) / (e1 + e1);
                let py = (e2 - sqr(px)).max(0.).sqrt();

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
                        // TODO this is e3
                        (sqr(e1), 0., sqr(e3), sigma_t, sigma_t + d_c_2)
                    }
                } else {
                    let [bend_left, bend_right] = data_driven_bending_heuristic(
                        e1,
                        fs[tgt_f]
                            .edges()
                            .map(|[a, b]| minmax(a, b))
                            .map(|e| edge_lens[&e]),
                        px,
                        py,
                        cx,
                        cy,
                        sx,
                        sy_neg,
                    );

                    if bend_left {
                        let sigma_t = prev_sigma_t + d_s_1;
                        (sqr(e1), 0., sqr(e3), sigma_t, sigma_t + d_c_1)
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

                    vert_to_src_dist[v0] = d_b;
                    vert_to_src_dist[v1] = d_a;
                    vert_to_src_dist[v2] = d_c;

                    let next_queue = if d_t < iter as F {
                        &mut queue_from
                    } else {
                        &mut queue_to
                    };
                    next_queue.push_back((eo_1, tgt_f, h2_behind));
                    next_queue.push_back((eo_2, tgt_f, h3_behind));
                }
            }
        }

        std::mem::swap(&mut queue_from, &mut queue_to);
    }

    println!("{updates}");

    let mut out_dists = vec![F::INFINITY; vs.len()];
    // for each edge, check face_to_extra_dist + vert_to_src_dist
    // also need to multiply by avg edge length
    for (e, adj_fs) in edge_face_adj.iter() {
      for &adj_f in adj_fs {
        let sigma = face_to_extra_dist[adj_f];
        for &v in e {
          let dist = vert_to_src_dist[v].sqrt();
          if !dist.is_finite() {
            continue;
          }
          out_dists[v] = out_dists[v].min(sigma + dist);
        }
      }
    }

    for d in out_dists.iter_mut() {
      *d *= avg_edge_len;
    }

    out_dists
}

pub fn data_driven_bending_heuristic(
    e_len: F,
    es: impl Iterator<Item = F>,
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

    let [max_e, min_e] = es.fold([F::INFINITY, F::NEG_INFINITY], |[l, h], n| {
        [l.min(n), h.max(n)]
    });

    let b0 = max_e > THRESH_C * e_len;
    let b1 = max_e > THRESH_G * min_e;

    let b2 = py < THRESH_H * e_len;
    let b3 = py < THRESH_HG * max_e;

    let idx = b0 as usize + (b1 as usize * 2) + (b2 as usize * 4) + (b3 as usize * 8);
    let l = LAMBDA[idx];

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

pub fn other(v: [usize; 3], e0: usize, e1: usize) -> usize {
    match v {
        [a, b, c] | [c, a, b] | [b, c, a] | [a, c, b] | [b, a, c] | [c, b, a]
            if b == e0 && c == e1 =>
        {
            a
        }
        _ => unreachable!(),
    }
}
