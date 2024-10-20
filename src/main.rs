#[allow(unused_imports)]
#[macro_use]
extern crate approx as approx;

extern crate faer as fa;
use core::f64;

mod elliptical_sweep;
mod multhopp;
mod plot;
mod rectangular_prandtl;
mod run_case;
mod wing;

const DEG_TO_RAD: f64 = 0.0174532925199;

fn main() {
    const AOA: f64 = 1.0 * DEG_TO_RAD;
    elliptical_sweep::elliptical_sweep();
    rectangular_prandtl::rectangular_prandtl();
    /*
    run_case::run_case(
        "out/enunciado",
        |y| wing::linear_law(y, 2.5, 1.2),
        |y| {
            AOA + wing::linear_law(y, 2.0 * DEG_TO_RAD, -2.0 * DEG_TO_RAD)
                + wing::ailerons(y, 0.1, 0.0 * DEG_TO_RAD)
        },
        |_| 5.1,
        40,
        12.0,
    );
    */
}
