#[allow(unused_imports)]
#[macro_use]
extern crate approx as approx;

extern crate faer as fa;
use core::f64;

use fa::{prelude::SpSolver, FaerMat};
use std::{fs::File, io::Write};

mod plot;

const DEG_TO_RAD: f64 = 0.0174532925199;

// Returns a delta alpha of attack to be added
// yp ranges from -1 to 1
fn ailerons(yp: f64, extend: f64, position: f64) -> f64 {
    if yp > 1.0 - extend {
        position
    } else if yp < -1.0 + extend {
        -position
    } else {
        0.0
    }
}

fn linear_law(yp: f64, root: f64, tip: f64) -> f64 {
    let tip_factor = yp.abs();
    return root * (1.0 - tip_factor) + tip * tip_factor;
}

fn cosine_samples(num_points: usize, b: f64) -> (Vec<f64>, Vec<f64>) {
    // We use the mapping y = b/2 * cos(theta)
    // Uniformly sample in cosine space (equivalent to cosine sampling in y space)
    // We don't include the extremes, as they are useless

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
            let t1 = 2.0 / (k * clalpha);
            let t2 = (n as f64) / (2.0 * theta.sin());
            let t = t1 + t2;
            mat[(j, n - 1)] = t * (n as f64 * theta).sin();
            // This is equivalent, as done in my TFG:
            /*mat[(j, n - 1)] = (4.0 * b / c) * ((n as f64) * theta).sin()
            + clalpha * (n as f64) * ((n as f64) * theta).sin() / theta.sin();*/
        }

        rhs[(j, 0)] = alpha;
        //rhs[(j, 0)] = clalpha * alpha;
    }

    let sln = mat.full_piv_lu().solve(rhs);

    return sln;
}

// S is int_(-b/2)^(b/2) c(y) dy
// by changing variable to y' = 2 y / b we find that
// S = b / 2 int_(-1)^(1) c(y') dy'
fn wing_area(cuerda: fn(f64) -> f64, b: f64) -> f64 {
    let dyp = 0.001;
    let mut int = 0.0;
    let mut yp = -1.0;
    while yp <= 1.0 {
        int += cuerda(yp) * dyp;
        yp += dyp;
    }
    0.5 * b * int
}

// Mean chord is that such that b * mean_chord = S
fn mean_chord(cuerda: fn(f64) -> f64, b: f64) -> f64 {
    wing_area(cuerda, b) / b
}

// We exploit the fact that AR = b^2 / S
// We use rectangle integration with heaps of points to get accurate
fn estimate_ar(cuerda: fn(f64) -> f64, b: f64) -> f64 {
    return b * b / wing_area(cuerda, b);
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

    // We assume maximum lift happens midpoint, we exceed by 20% just in case
    let mag = {
        let mut sum = 0.0;
        for i in 0..sln.nrows() {
            let n = (i + 1) as f64;
            sum += sln[(i, 0)] * (n * std::f64::consts::PI * 0.5).sin();
        }
        2.0 * b * sum / mean_chord(cuerda, b) * 1.2
    };

    {
        let plotname = [name, "-cl.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b, -mag, mag);
        plot::add_cl(&mut ctx, &sln, mean_chord(cuerda, b));
        plot::finish(ctx);
    }
    /*{
        let plotname = [name, "-geom.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b);
        plot::add_fn("cuerda", &mut ctx, cuerda, 0.05, true);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-alpha.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b);
        plot::add_fn("alpha", &mut ctx, alpha, 0.05, false);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-cd.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b);
        plot::add_cd(&mut ctx, &sln);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-cmx.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b);
        plot::add_cmx(&mut ctx, &sln);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-cmy.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b);
        plot::add_cmy(&mut ctx, &sln);
        plot::finish(ctx);
    }*/
}

fn main() {
    const AOA: f64 = 0.0 * DEG_TO_RAD;
    run_case(
        "out/rectangular",
        |_| 1.0,
        |_| AOA,
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
        |_| AOA,
        |_| 2.0 * std::f64::consts::PI,
        45,
        1.0,
    );

    run_case(
        "out/enunciado",
        |y| linear_law(y, 2.5, 1.2),
        |y| {
            AOA + linear_law(y, 2.0 * DEG_TO_RAD, -2.0 * DEG_TO_RAD)
                + ailerons(y, 0.1, 4.0 * DEG_TO_RAD)
        },
        |_| 5.1,
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
        assert_relative_eq!(estimate_ar(|_| 1.0, 1.0), 1.0, max_relative = 0.001);
        assert_relative_eq!(estimate_ar(|_| 2.0, 1.0), 0.5, max_relative = 0.001);
        assert_relative_eq!(estimate_ar(|_| 1.0, 2.0), 2.0, max_relative = 0.001);
    }

    #[test]
    fn ar_trapezoid() {
        // Mean chord is that such that b * SMC = S
        // Thus in a trapezoid, with S = b * c_end + 0.5 * b * (c_start - c_end)
        // we have that SMC = c_end + 0.5 * (c_start - c_end)
        // (AR = b / SMC)
        assert_relative_eq!(
            estimate_ar(|yp| 1.0 - yp.abs() * 0.5, 1.0),
            1.0 / (0.5 + 0.5 * (1.0 - 0.5)),
            max_relative = 0.001
        );
        assert_relative_eq!(
            estimate_ar(|yp| 2.0 - yp.abs(), 1.0),
            1.0 / (1.0 + 0.5 * (2.0 - 1.0)),
            max_relative = 0.001
        );
        assert_relative_eq!(
            estimate_ar(|yp| 3.0 - yp.abs() * 0.5, 2.0),
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
            estimate_ar(|yp| (1.0 - yp * yp).sqrt(), 1.0),
            1.0 / (0.25 * std::f64::consts::PI * 1.0),
            max_relative = 0.001
        );
        assert_relative_eq!(
            estimate_ar(|yp| (1.0 - yp * yp).sqrt() * 3.0, 1.0),
            1.0 / (0.25 * std::f64::consts::PI * 3.0),
            max_relative = 0.001
        );
        assert_relative_eq!(
            estimate_ar(|yp| (1.0 - yp * yp).sqrt(), 2.0),
            2.0 / (0.25 * std::f64::consts::PI * 1.0),
            max_relative = 0.001
        );
    }
}
