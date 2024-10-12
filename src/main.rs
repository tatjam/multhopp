#[macro_use]
extern crate approx as approx;

extern crate faer as fa;
use core::{f64, num};

use fa::{prelude::SpSolver, FaerMat};
use plotters::prelude::{DiscreteRanged, IntoLinspace};
use std::{fs::File, io::Write};

mod plot;

fn cosine_samples(num_points: usize, b: f64) -> (Vec<f64>, Vec<f64>) {
    // We use the mapping y = b/2 * cos(theta)
    // Uniformly sample in cosine space (equivalent to cosine sampling in y space)
    // We don't include the extremes, as they are useless
    /*let cos_thetas: Vec<f64> = (0..num_points)
        .into_iter()
        .map(|i| ((i as f64 + 1.0) / (num_points as f64 + 1.0) - 0.5) * 2.0)
        .collect();
    let thetas: Vec<f64> = cos_thetas.iter().map(|c| c.acos()).collect();
    let ys: Vec<f64> = cos_thetas.iter().map(|c| b / 2.0 * c).collect();*/

    // Linear sampling in theta space, from nearly 0 to nearly PI
    let thetas: Vec<f64> = (0..num_points)
        .into_iter()
        .map(|i| (i as f64 + 1.0) / (num_points as f64 + 1.0) * std::f64::consts::PI)
        .collect();
    // Ys thus go from nearly b / 2.0 * cos(0) to nearly b / 2.0 * cos(PI)
    let ys: Vec<f64> = thetas.iter().map(|c| b / 2.0 * c.cos()).collect();

    return (thetas, ys);
}

// Returns A coefficients, given functions take non dimensional
// spanwise position that ranges from -1 to 1
fn multhopp(
    cuerda: fn(f64) -> f64,
    alpha: fn(f64) -> f64,
    clalpha: fn(f64) -> f64,
    num_points: usize,
    b: f64,
) -> fa::Mat<f64> {
    let (thetas, ys) = cosine_samples(num_points, b);
    let ys_adim: Vec<f64> = ys.iter().map(|y| 2.0 * y / b).collect();

    let c_samples: Vec<f64> = ys_adim.iter().map(|y| cuerda(*y)).collect();
    let alpha_samples: Vec<f64> = ys_adim.iter().map(|y| alpha(*y)).collect();
    let clalpha_samples: Vec<f64> = ys_adim.iter().map(|y| clalpha(*y)).collect();

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
            /*let t1 = 2.0 / (k * clalpha);
            let t2 = (n as f64) / (2.0 * theta.sin());
            let t = t1 + t2;
            mat[(j, n - 1)] = t * (n as f64 * theta).sin();*/
            mat[(j, n - 1)] = (4.0 * b / c) * ((n as f64) * theta).sin()
                + clalpha * (n as f64) * ((n as f64) * theta).sin() / theta.sin();
        }

        rhs[(j, 0)] = clalpha * alpha;
    }

    println!("{:?}{:?}", mat, rhs);

    let sln = mat.full_piv_lu().solve(rhs);

    return sln;
}

// We exploit the fact that AR = b^2 / S
// S is int_(-b/2)^(b/2) c(y) dy
// by changing variable to y' = 2 y / b we find that
// b^2 / S = b^2 / b int_(-1)^(1) c(y') dy' = 2 b / int_(-0.5)^(0.5) c(y') dy'
// We use rectangle integration with heaps of points to get accurate
fn estimate_AR(cuerda: fn(f64) -> f64, b: f64) -> f64 {
    let dyp = 0.001;
    let mut int = 0.0;
    let mut yp = -1.0;
    while yp <= 1.0 {
        int += cuerda(yp) * dyp;
        yp += dyp;
    }
    return 2.0 * b / int;
}

fn run_case(
    name: &str,
    cuerda: fn(f64) -> f64,
    alpha: fn(f64) -> f64,
    clalpha: fn(f64) -> f64,
    num_points: usize,
    b: f64,
) {
    let fname = [name, ".txt"].concat();
    let mut ofile = File::options()
        .write(true)
        .append(false)
        .truncate(true)
        .create(true)
        .open(fname)
        .unwrap();

    let sln = multhopp(cuerda, alpha, clalpha, num_points, b);

    writeln!(&mut ofile, "{:?}", sln).unwrap();
    let k = {
        let mut sum = 0.0;
        let base = sln[(0, 0)] * sln[(0, 0)];
        // i = 1 means n = 2!
        for i in 1..sln.nrows() {
            let n = (i + 1) as f64;
            sum += n * sln[(i, 0)] * sln[(i, 0)] / base;
        }
        1.0 + sum
    };
    writeln!(&mut ofile, "ostwald k = {}", k).unwrap();

    let plotname = [name, ".png"].concat();
    let mut ctx = plot::make_default_plot(&plotname, b);
    plot::add_cl(&mut ctx, &sln);
    plot::add_geom("cuerda", &mut ctx, cuerda, 0.05);
    plot::finish(ctx);
}

fn main() {
    run_case(
        "out/rectangular",
        |_| 1.0,
        |_| 0.1,
        |_| 2.0 * std::f64::consts::PI,
        45,
        1.0,
    );
    run_case(
        "out/elliptical",
        |y| {
            // y goes from -1 to 1
            (1.0 - y * y).sqrt()
        },
        |_| 0.1,
        |_| 2.0 * std::f64::consts::PI,
        45,
        1.0,
    );

    run_case(
        "out/enunciado",
        |y| 2.5 - y.abs() * (2.5 - 1.2),
        |y| 0.5,
        |y| 5.1,
        40,
        12.0,
    );
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn ar_rectangles() {
        // In rectangles, AR is simply b / c
        assert_relative_eq!(estimate_AR(|_| 1.0, 1.0), 1.0, max_relative = 0.001);
        assert_relative_eq!(estimate_AR(|_| 2.0, 1.0), 0.5, max_relative = 0.001);
        assert_relative_eq!(estimate_AR(|_| 1.0, 2.0), 2.0, max_relative = 0.001);
    }

    #[test]
    fn ar_trapezoid() {
        // Mean chord is that such that b * SMC = S
        // Thus in a trapezoid, with S = b * c_end + 0.5 * b * (c_start - c_end)
        // we have that SMC = c_end + 0.5 * (c_start - c_end)
        // (AR = b / SMC)
        assert_relative_eq!(
            estimate_AR(|yp| 1.0 - yp.abs() * 0.5, 1.0),
            1.0 / (0.5 + 0.5 * (1.0 - 0.5)),
            max_relative = 0.001
        );
        assert_relative_eq!(
            estimate_AR(|yp| 2.0 - yp.abs(), 1.0),
            1.0 / (1.0 + 0.5 * (2.0 - 1.0)),
            max_relative = 0.001
        );
        assert_relative_eq!(
            estimate_AR(|yp| 3.0 - yp.abs() * 0.5, 2.0),
            2.0 / (2.5 + 0.5 * (3.0 - 2.5)),
            max_relative = 0.001
        );
    }

    #[test]
    fn ar_ellipses() {
        // The surface area of a semi-ellipse of semiaxes A, B is 0.5 * pi * A * B, thus
        // SMC = 0.5 * pi * A * B / b. Now, B = b / 2, thus
        // SMC = 0.25 * pi * A
        // Hence for an ellipse AR = b / (0.25 * pi * A)
        assert_relative_eq!(
            estimate_AR(|yp| (1.0 - yp * yp).sqrt(), 1.0),
            1.0 / (0.25 * std::f64::consts::PI * 1.0),
            max_relative = 0.001
        );
        assert_relative_eq!(
            estimate_AR(|yp| (1.0 - yp * yp).sqrt() * 3.0, 1.0),
            1.0 / (0.25 * std::f64::consts::PI * 3.0),
            max_relative = 0.001
        );
        assert_relative_eq!(
            estimate_AR(|yp| (1.0 - yp * yp).sqrt(), 2.0),
            2.0 / (0.25 * std::f64::consts::PI * 1.0),
            max_relative = 0.001
        );
    }
}
