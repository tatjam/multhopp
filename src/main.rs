extern crate approx as approx;
extern crate faer as fa;
use core::{f64, num};

use approx::*;
use fa::{prelude::SpSolver, FaerMat};
use plotters::prelude::*;

mod plot;

fn cosine_samples(num_points: usize, b: f64) -> (Vec<f64>, Vec<f64>) {
    // We use the mapping y = b/2 * cos(theta)
    // Uniformly sample in cosine space (equivalent to cosine sampling in y space)
    // We don't include the extremes, as they are useless
    let cos_thetas: Vec<f64> = (0..num_points)
        .into_iter()
        .map(|i| ((i as f64 + 1.0) / (num_points as f64 + 1.0) - 0.5) * 2.0)
        .collect();
    let thetas: Vec<f64> = cos_thetas.iter().map(|c| c.acos()).collect();
    let ys: Vec<f64> = cos_thetas.iter().map(|c| b / 2.0 * c).collect();

    return (thetas, ys);
}

// Returns A coefficients
fn multhopp(
    cuerda: fn(f64) -> f64,
    alpha: fn(f64) -> f64,
    clalpha: fn(f64) -> f64,
    num_points: usize,
    b: f64,
) -> fa::Mat<f64> {
    let (thetas, ys) = cosine_samples(num_points, b);

    let c_samples: Vec<f64> = ys.iter().map(|y| cuerda(*y)).collect();
    let alpha_samples: Vec<f64> = ys.iter().map(|y| alpha(*y)).collect();
    let clalpha_samples: Vec<f64> = ys.iter().map(|y| clalpha(*y)).collect();

    let mut mat = fa::Mat::<f64>::zeros(num_points, num_points);
    let mut rhs = fa::Mat::<f64>::zeros(num_points, 1);

    for j in 0..num_points {
        // theta ranges from nearly 0 to nearly pi
        let theta = thetas[j];
        let c = c_samples[j];
        let clalpha = clalpha_samples[j];
        let alpha = alpha_samples[j];
        let k = c / b;

        for n in 1..(num_points + 1) {
            let t1 = 2.0 / (k * clalpha);
            let t2 = (n as f64) / (2.0 * theta.sin());
            let t = t1 + t2;
            mat[(j, n - 1)] = t * (n as f64 * theta).sin();
        }

        rhs[(j, 0)] = alpha;
    }

    println!("{:?}", mat);
    println!("{:?}", rhs);

    let sln = mat.full_piv_lu().solve(rhs);

    return sln;
}

const B: f64 = 12.0;

fn c_rectangular(y: f64) -> f64 {
    10.0
}

fn main() {
    let root_area = BitMapBackend::new("out/test.png", (1024, 768)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let sln = multhopp(c_rectangular, |_| 0.15, |_| 5.1, 15, B);
    let b_range = (-B / 2.0)..(B / 2.0);

    println!("{:?}", sln);

    let mut cc = ChartBuilder::on(&root_area)
        .margin(5)
        .set_all_label_area_size(50)
        .build_cartesian_2d(b_range.clone(), -0.2f32..0.2)
        .unwrap();

    cc.configure_mesh()
        .x_labels(20)
        .y_labels(15)
        .draw()
        .unwrap();

    let b_range_step = b_range.step(0.001);

    cc.draw_series(LineSeries::new(
        b_range_step.values().map(|y| {
            let theta = (2.0 * y / B).acos();
            let mut sum = 0.0;
            for i in 0..sln.nrows() {
                let n = i + 1;
                sum += sln[(i, 0)] * ((n as f64) * theta).sin();
            }
            (y, sum as f32)
        }),
        &RED,
    ))
    .unwrap()
    .label("cl");

    cc.draw_series(LineSeries::new(
        b_range_step
            .values()
            .map(|y| (y, c_rectangular(y) as f32 * 0.001)),
        &BLUE,
    ))
    .unwrap()
    .label("cl");

    root_area.present().unwrap();
}
