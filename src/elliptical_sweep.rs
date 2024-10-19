use core::fmt;

use full_palette::{BROWN, LIGHTBLUE, LIGHTGREEN};
use plotters::prelude::*;

use super::*;

pub fn elliptical_sweep() {
    const AOA: f64 = 5.0 * DEG_TO_RAD;
    const CL_ALPHA: f64 = 2.0 * std::f64::consts::PI;

    const NPOINTS: [usize; 5] = [3, 5, 15, 25, 45];
    const AR: [f64; 8] = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 12.0];

    // num_points, k
    let mut k_pairs: Vec<(usize, f64)> = Vec::new();
    // num_points, ar, cl, cdi
    let mut ar_pairs: Vec<(usize, Vec<(f64, f64, f64)>)> = Vec::new();

    // k sweep
    for num_points in NPOINTS {
        let sln = multhopp::multhopp(
            |y| (1.0 - y * y).sqrt(),
            |_| AOA,
            |_| CL_ALPHA,
            num_points,
            1.0,
        );

        let k = run_case::estimate_ostwald_factor(&sln);
        k_pairs.push((num_points, k));
    }
    println!("k_pairs: {:#?}", k_pairs);

    // CL and CDi sweep, by AR and by number of points
    for num_points in NPOINTS {
        let mut vecn: Vec<(f64, f64, f64)> = Vec::new();

        for ar in AR {
            // AR of an elliptical wing is b^2 / S
            // but S = 0.5 * pi * (b / 2) * 1, thus
            // (We fix chord to 1)
            // AR = b / (0.25 * pi * 1), solving for b
            let b = ar * 0.25 * std::f64::consts::PI;
            let sln = multhopp::multhopp(
                |y| (1.0 - y * y).sqrt(),
                |_| AOA,
                |_| CL_ALPHA,
                num_points,
                b,
            );
            let ar_estim = multhopp::estimate_ar(|y| (1.0 - y * y).sqrt(), b);
            assert_relative_eq!(ar_estim, ar, max_relative = 0.01);

            let cl = run_case::estimate_cl(&sln, ar_estim);
            let cdi = run_case::estimate_cdi(&sln, ar_estim);
            vecn.push((ar, cl, cdi));
        }
        ar_pairs.push((num_points, vecn));
    }

    println!("{:#?}", ar_pairs);

    // and nother of cl and cdi over ar, with different colors for each num_points
    {
        let root_area =
            BitMapBackend::new("out/elliptical-sweep.png", (1024, 768)).into_drawing_area();
        root_area.fill(&WHITE).unwrap();

        let mut cc = ChartBuilder::on(&root_area)
            .margin(5)
            .set_all_label_area_size(50)
            .build_cartesian_2d(AR[0]..AR[AR.len() - 1], 0.0..1.5)
            .unwrap();

        cc.configure_mesh()
            .x_labels(20)
            .y_labels(15)
            .draw()
            .unwrap();

        for (j, v) in ar_pairs.iter().enumerate() {
            let i = j + 1;
            cc.draw_series(
                v.1.iter()
                    .map(|(ar, cl, cdi)| Circle::new((*ar, *cl), i as f64, &BLUE)),
            )
            .unwrap()
            .label(format!("cl num_points = {}", v.0))
            .legend(move |(x, y)| Circle::new((x + 10, y), i as f64, &BLUE));
            cc.draw_series(
                v.1.iter()
                    .map(|(ar, cl, cdi)| Circle::new((*ar, *cdi * 100.0), i as f64, &GREEN)),
            )
            .unwrap()
            .label(format!("cdi num_points = {}", v.0))
            .legend(move |(x, y)| Circle::new((x + 10, y), i as f64, &GREEN));
        }

        // Analytical series
        cc.draw_series(LineSeries::new(
            (AR[0]..AR[AR.len() - 1]).step(0.1).values().map(|ar| {
                (
                    ar,
                    AOA * CL_ALPHA / (1.0 + CL_ALPHA / (std::f64::consts::PI * ar)),
                )
            }),
            &LIGHTBLUE,
        ))
        .unwrap()
        .label("CL (ideal)")
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &LIGHTBLUE));

        cc.draw_series(LineSeries::new(
            (AR[0]..AR[AR.len() - 1]).step(0.1).values().map(|ar| {
                (
                    ar,
                    100.0
                        * (AOA * CL_ALPHA / (1.0 + CL_ALPHA / (std::f64::consts::PI * ar))).powi(2)
                        / (std::f64::consts::PI * ar),
                )
            }),
            &LIGHTGREEN,
        ))
        .unwrap()
        .label("CDi (ideal)")
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &LIGHTGREEN));

        cc.configure_series_labels()
            .border_style(&BLACK)
            .background_style(&WHITE.mix(0.8))
            .position(SeriesLabelPosition::UpperRight)
            .draw()
            .unwrap();
    }

    {}
}
