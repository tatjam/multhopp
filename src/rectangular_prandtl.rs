use core::fmt;

use full_palette::{BROWN, GREEN_400, LIGHTBLUE, LIGHTGREEN};
use plotters::prelude::*;

use super::*;

pub fn rectangular_prandtl() {
    const AOA: f64 = 5.0 * DEG_TO_RAD;
    const CL_ALPHA: f64 = 2.0 * std::f64::consts::PI;

    const NPOINTS: [usize; 5] = [1, 3, 5, 7, 15];
    const AR_OVER_2PI: [f64; 7] = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75];

    // num_points, ar, cl, cdi
    let mut ar_pairs: Vec<(usize, Vec<(f64, f64, f64)>)> = Vec::new();

    // CL and CDi sweep, by AR and by number of points
    for num_points in NPOINTS {
        let mut vecn: Vec<(f64, f64, f64)> = Vec::new();

        for ar_over_2pi in AR_OVER_2PI {
            let ar = ar_over_2pi * (2.0 * std::f64::consts::PI);
            // AR of rectangular wing with chord = 1 is simply
            // b = ar
            let b = ar;
            let sln = multhopp::multhopp(|y| 1.0, |_| AOA, |_| CL_ALPHA, num_points, b);
            let ar_estim = multhopp::estimate_ar(|y| 1.0, b);
            assert_relative_eq!(ar_estim, ar, max_relative = 0.01);

            let cl = run_case::estimate_cl(&sln, ar_estim);
            let cdi = run_case::estimate_cdi(&sln, ar_estim);
            vecn.push((ar, cl, cdi));
        }
        ar_pairs.push((num_points, vecn));
    }

    // and nother of cl and cdi over ar, with different colors for each num_points
    {
        /*let root_area =
        BitMapBackend::new("out/rectangular-sweep.png", (1024, 768)).into_drawing_area();*/
        let root_area =
            SVGBackend::new("out/rectangular-sweep.svg", (1024, 768)).into_drawing_area();
        root_area.fill(&WHITE).unwrap();

        let mut cc = ChartBuilder::on(&root_area)
            .margin(5)
            .set_all_label_area_size(50)
            .build_cartesian_2d(AR_OVER_2PI[0]..AR_OVER_2PI[AR_OVER_2PI.len() - 1], 0.0..0.6)
            .unwrap()
            .set_secondary_coord(
                AR_OVER_2PI[0]..AR_OVER_2PI[AR_OVER_2PI.len() - 1],
                0.004..0.014,
            );

        cc.configure_mesh()
            .x_labels(5)
            .y_labels(20)
            .x_desc("AR / (2 pi)")
            .y_desc("CL")
            .draw()
            .unwrap();

        cc.configure_secondary_axes().y_desc("CDi").draw().unwrap();

        for (j, v) in ar_pairs.iter().enumerate() {
            let i = j + 1;
            cc.draw_series(v.1.iter().map(|(ar, cl, cdi)| {
                Circle::new((*ar / (2.0 * std::f64::consts::PI), *cl), i as f64, &BLUE)
            }))
            .unwrap()
            .label(format!("cl num_points = {}", v.0))
            .legend(move |(x, y)| Circle::new((x + 10, y), i as f64, &BLUE));
            cc.draw_secondary_series(v.1.iter().map(|(ar, cl, cdi)| {
                Circle::new((*ar / (2.0 * std::f64::consts::PI), *cdi), i as f64, &GREEN)
            }))
            .unwrap()
            .label(format!("cdi num_points = {}", v.0))
            .legend(move |(x, y)| Circle::new((x + 10, y), i as f64, &GREEN));
        }

        // From "The Elements of Aerofoil and Airscrew Theory", H. Glauert
        let glauert = [
            (0.25, 0.426, 0.007),
            (0.5, 0.587, 0.019),
            (0.75, 0.675, 0.034),
            (1.0, 0.729, 0.049),
            (1.25, 0.767, 0.063),
            (1.5, 0.794, 0.076),
            (1.75, 0.815, 0.088),
        ];

        cc.draw_series(glauert.iter().map(|(ar_over_2pi, a_over_a0, delta)| {
            let ar = 2.0 * std::f64::consts::PI * ar_over_2pi;
            let cl = 2.0 * std::f64::consts::PI * a_over_a0 * AOA;
            Cross::new((*ar_over_2pi, cl), 8.0, &BLACK)
        }))
        .unwrap()
        .label(format!("cl glauert"))
        .legend(move |(x, y)| Cross::new((x + 10, y), 5.0, &BLACK));

        cc.draw_secondary_series(glauert.iter().map(|(ar_over_2pi, a_over_a0, delta)| {
            let ar = 2.0 * std::f64::consts::PI * ar_over_2pi;
            let cl = 2.0 * std::f64::consts::PI * a_over_a0 * AOA;
            let cd = cl * cl / (std::f64::consts::PI * ar) * (1.0 + delta);
            Cross::new((*ar_over_2pi, cd), 8.0, &BLACK)
        }))
        .unwrap()
        .label(format!("cdi glauert"))
        .legend(move |(x, y)| Cross::new((x + 10, y), 5.0, &BLACK));

        cc.configure_series_labels()
            .border_style(&BLACK)
            .background_style(&WHITE.mix(0.8))
            .position(SeriesLabelPosition::LowerLeft)
            .draw()
            .unwrap();
    }

    {}
}
