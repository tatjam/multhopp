fn elliptical_sweep() {
    run_case(
        "out/elliptical",
        |y| {
            // y goes from -1 to 1
            (1.0 - y * y).sqrt()
        },
        |_| AOA + 5.0 * DEG_TO_RAD,
        |_| 2.0 * std::f64::consts::PI,
        45,
        1.0,
    );
}
