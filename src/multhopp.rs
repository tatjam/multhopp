use fa::{prelude::SpSolver, FaerMat};

// Returns A coefficients, given functions take non dimensional
// spanwise position that ranges from -1 to 1
pub fn multhopp(
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
        // If TFG method is used, write this instead:
        //rhs[(j, 0)] = clalpha * alpha;
    }

    let sln = mat.full_piv_lu().solve(rhs);

    return sln;
}

// S is int_(-b/2)^(b/2) c(y) dy
// by changing variable to y' = 2 y / b we find that
// S = b / 2 int_(-1)^(1) c(y') dy'
pub fn wing_area(cuerda: fn(f64) -> f64, b: f64) -> f64 {
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
pub fn mean_chord(cuerda: fn(f64) -> f64, b: f64) -> f64 {
    wing_area(cuerda, b) / b
}

// We exploit the fact that AR = b^2 / S
// We use rectangle integration with heaps of points to get accurate
pub fn estimate_ar(cuerda: fn(f64) -> f64, b: f64) -> f64 {
    return b * b / wing_area(cuerda, b);
}

pub fn cosine_samples(num_points: usize, b: f64) -> (Vec<f64>, Vec<f64>) {
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

pub struct AeroCoefficients {
    cl: f64,
    cd: f64,
    cmx: f64,
    cmz: f64,
}

pub fn calculate_coefficients(sln: &fa::Mat<f64>, ar: f64) -> AeroCoefficients {
    use std::f64::consts::PI;
    let cl = PI * ar / 2.0 * sln[(0, 0)];
    let cd = cl * cl / (PI * ar)
        * (1.0 + {
            let mut sum = 0.0;
            for i in 1..sln.nrows() {
                let n = (i + 1) as f64;
                sum += n * sln[(i, 0)] * sln[(i, 0)];
            }
            sum / (sln[(0, 0)] * sln[(0, 0)])
        });

    let cmx = PI * ar / 8.0 * sln[(1, 0)];
    let cmz = -PI * ar / 16.0 * {
        let mut sum = 0.0;
        for i in 0..(sln.nrows() - 1) {
            let n = (i + 1) as f64;
            sum += (2.0 * n + 1.0) * sln[(i, 0)] * sln[(i + 1, 0)];
        }
        sum
    };

    return AeroCoefficients {
        cl: cl,
        cd: cd,
        cmx: cmx,
        cmz: cmz,
    };
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
