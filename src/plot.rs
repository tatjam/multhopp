use plotters::coord::*;
use plotters::prelude::*;
use types::RangedCoordf32;
use types::RangedCoordf64;

pub struct PlotContext<'a> {
    cc: ChartContext<'a, BitMapBackend<'a>, Cartesian2d<RangedCoordf64, RangedCoordf32>>,
    root_area: DrawingArea<BitMapBackend<'a>, Shift>,
    b_range: std::ops::Range<f64>,
}

pub fn make_default_plot<'a>(name: &'a str, b: f64) -> PlotContext<'a> {
    let b_range = (-b / 2.0)..(b / 2.0);

    let root_area = BitMapBackend::new(name, (1024, 768)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut cc = ChartBuilder::on(&root_area)
        .margin(5)
        .set_all_label_area_size(50)
        .build_cartesian_2d(b_range.clone(), -(0.2 * b) as f32..(0.2 * b) as f32)
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

pub fn add_cl(ctx: &mut PlotContext, sln: &fa::Mat<f64>) {
    let b_range_step = ctx.b_range.clone().step(0.001);
    let b = ctx.b_range.end - ctx.b_range.start;
    ctx.cc
        .draw_series(LineSeries::new(
            b_range_step.values().map(|y| {
                let theta = (2.0 * y / b).acos();
                let mut sum = 0.0;
                for i in 0..sln.nrows() {
                    let n = i + 1;
                    sum += sln[(i, 0)] * ((n as f64) * theta).sin();
                }
                // We assume incoming air speed of 1m/s
                sum *= 1.225 * b;
                (y, sum as f32)
            }),
            &RED,
        ))
        .unwrap()
        .label("cl");
}

pub fn add_geom(name: &str, ctx: &mut PlotContext, cuerda: fn(f64) -> f64, scale: f64) {
    let b_range_step = ctx.b_range.clone().step(0.001);
    let b = ctx.b_range.end - ctx.b_range.start;
    // We place the top line such that it's 1/4 of the chord
    ctx.cc
        .draw_series(LineSeries::new(
            b_range_step.values().map(|y| {
                let adimy = 2.0 * y / b;
                (y, (cuerda(adimy) * scale * 0.25 * b) as f32)
            }),
            &BLUE,
        ))
        .unwrap()
        .label(name);

    // And the bottom line such that it's the remainder
    ctx.cc
        .draw_series(LineSeries::new(
            b_range_step.values().map(|y| {
                let adimy = 2.0 * y / b;
                (y, (-cuerda(adimy) * scale * 0.75 * b) as f32)
            }),
            &BLUE,
        ))
        .unwrap()
        .label(name);
}

// This "consumes" the plot context
pub fn finish(ctx: PlotContext) {
    ctx.root_area.present().unwrap();
}
