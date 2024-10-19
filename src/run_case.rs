use super::*;
use std::{fs::File, io::Write};

pub fn run_case(
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

    let sln = multhopp::multhopp(cuerda, alpha, clalpha, num_points, b);

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

    let mean_c = multhopp::mean_chord(cuerda, b);

    // We assume maximum lift happens midpoint, we exceed by 20% just in case
    let mag = {
        let mut sum = 0.0;
        for i in 0..sln.nrows() {
            let n = (i + 1) as f64;
            sum += sln[(i, 0)] * (n * std::f64::consts::PI * 0.5).sin();
        }
        2.0 * b * sum / mean_c * 1.2
    };

    {
        let plotname = [name, "-cl.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b, -mag, mag);
        plot::add_cl(&mut ctx, &sln, mean_c);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-geom.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b, -mean_c, mean_c);
        plot::add_fn("cuerda", &mut ctx, cuerda, 0.05, true);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-alpha.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b, -0.2, 0.2);
        plot::add_fn("alpha", &mut ctx, alpha, 0.05, false);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-cd.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b, -mag * 0.03, 0.01);
        plot::add_cd(&mut ctx, &sln, mean_c);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-cmx.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b, -mag, mag);
        plot::add_cmx(&mut ctx, &sln, mean_c);
        plot::finish(ctx);
    }
    {
        let plotname = [name, "-cmy.png"].concat();
        let mut ctx = plot::make_default_plot(&plotname, b, -mag * 0.03, mag * 0.03);
        plot::add_cmy(&mut ctx, &sln, mean_c);
        plot::finish(ctx);
    }
}
