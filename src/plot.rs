use plotters::coord::*;
use plotters::prelude::*;
use types::RangedCoordf32;
use types::RangedCoordf64;

pub struct PlotContext<'a> {
    cc: ChartContext<'a, BitMapBackend<'a>, Cartesian2d<RangedCoordf64, RangedCoordf32>>,
    root_area: DrawingArea<BitMapBackend<'a>, Shift>,
    b_range: std::ops::Range<f64>,
}

pub fn make_default_plot<'a>(name: &'a str, b: f64, y0: f64, y1: f64) -> PlotContext<'a> {
    let b_range = (-b / 2.0)..(b / 2.0);

    let root_area = BitMapBackend::new(name, (1024, 768)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut cc = ChartBuilder::on(&root_area)
        .margin(5)
        .set_all_label_area_size(50)
        .build_cartesian_2d(b_range.clone(), (y0 as f32)..(y1 as f32))
        .unwrap();

    cc.configure_mesh()
        .x_labels(20)
        .y_labels(15)
        .draw()
        .unwrap();

    return PlotContext {
        cc: cc,
        root_area: root_area,
        b_range: b_range,
    };
}

fn cl(sln: &fa::Mat<f64>, theta: f64, b: f64, mean_c: f64) -> f64 {
    let mut sum = 0.0;
    // l = rho * V_infty * gamma,
    // gamma = b V_infty sum(A_n * sin(n theta))
    // substituting
    // l = rho * b * V_infty^2 * sum(A_n * sin(n theta))
    // but c_l = l / (0.5 * rho * V_infty^2 * c)
    // c_l = 2.0 * b / c * sum(A_n * sin(n theta))
    // Note that we use the mean chord, and not each sections' chord
    for i in 0..sln.nrows() {
        let n = i + 1;
        sum += sln[(i, 0)] * ((n as f64) * theta).sin();
    }
    // We now multiply by the other terms
    sum *= 2.0 * b / mean_c;
    return sum;
}

fn cm(sln: &fa::Mat<f64>, theta: f64, b: f64, mean_c: f64) -> f64 {
    let mut sum = 0.0;
    // cl * induced_alpha =
    // cL * - 1 / (2 * sin(theta)) * sum(n An sin(n theta))
    // =  -2.0 * b / c / (2 * sin(theta)) * sum(A_n * sin(n theta)) * sum(A_n * sin(n theta))
    for i in 0..sln.nrows() {
        let n = i + 1;
        sum += sln[(i, 0)] * ((n as f64) * theta).sin();
    }
    // square it
    sum *= sum;
    // We now multiply by the other terms
    sum *= -2.0 * b / mean_c / (2.0 * theta.sin());
    return sum;
}

pub fn add_cl(ctx: &mut PlotContext, sln: &fa::Mat<f64>, mean_c: f64) {
    let b_range_step = ctx.b_range.clone().step(0.001);
    let b = ctx.b_range.end - ctx.b_range.start;
    ctx.cc
        .draw_series(LineSeries::new(
            b_range_step.values().map(|y| {
                let theta = (2.0 * y / b).acos();
                (y, cl(sln, theta, b, mean_c) as f32)
            }),
            &RED,
        ))
        .unwrap()
        .label("cl");
}

pub fn add_cd(ctx: &mut PlotContext, sln: &fa::Mat<f64>, mean_c: f64) {
    let b_range_step = ctx.b_range.clone().step(0.001);
    let b = ctx.b_range.end - ctx.b_range.start;
    ctx.cc
        .draw_series(LineSeries::new(
            b_range_step.values().map(|y| {
                let theta = (2.0 * y / b).acos();
                (y, cm(sln, theta, b, mean_c) as f32)
            }),
            &RED,
        ))
        .unwrap()
        .label("dl");
}

// Simply the cl times arm, further by span (TODO: or mean chord) to make non-dimensional
pub fn add_cmx(ctx: &mut PlotContext, sln: &fa::Mat<f64>, mean_c: f64) {
    let b_range_step = ctx.b_range.clone().step(0.001);
    let b = ctx.b_range.end - ctx.b_range.start;
    ctx.cc
        .draw_series(LineSeries::new(
            b_range_step.values().map(|y| {
                let theta = (2.0 * y / b).acos();
                (y, (cl(sln, theta, b, mean_c) * y / b) as f32)
            }),
            &RED,
        ))
        .unwrap()
        .label("cmx");
}

pub fn add_cmy(ctx: &mut PlotContext, sln: &fa::Mat<f64>, mean_c: f64) {
    let b_range_step = ctx.b_range.clone().step(0.001);
    let b = ctx.b_range.end - ctx.b_range.start;
    ctx.cc
        .draw_series(LineSeries::new(
            b_range_step.values().map(|y| {
                let theta = (2.0 * y / b).acos();
                (y, (cm(sln, theta, b, mean_c) * y / b) as f32)
            }),
            &RED,
        ))
        .unwrap()
        .label("cmy");
}

pub fn add_fn(
    name: &str,
    ctx: &mut PlotContext,
    cuerda: fn(f64) -> f64,
    scale: f64,
    quarter: bool,
) {
    let b_range_step = ctx.b_range.clone().step(0.001);
    let b = ctx.b_range.end - ctx.b_range.start;

    let (fseg, sseg) = if quarter { (0.25, 0.75) } else { (1.0, 0.0) };

    // We place the top line such that it's 1/4 of the chord
    ctx.cc
        .draw_series(LineSeries::new(
            b_range_step.values().map(|y| {
                let adimy = 2.0 * y / b;
                (y, (cuerda(adimy) * scale * fseg * b) as f32)
            }),
            &BLUE,
        ))
        .unwrap()
        .label(name);

    if quarter {
        // And the bottom line such that it's the remainder
        ctx.cc
            .draw_series(LineSeries::new(
                b_range_step.values().map(|y| {
                    let adimy = 2.0 * y / b;
                    (y, (-cuerda(adimy) * scale * sseg * b) as f32)
                }),
                &BLUE,
            ))
            .unwrap()
            .label(name);
    }
}

// This "consumes" the plot context
pub fn finish(ctx: PlotContext) {
    ctx.root_area.present().unwrap();
}
