// Returns a delta alpha of attack to be added
// yp ranges from -1 to 1
pub fn ailerons(yp: f64, extend: f64, position: f64) -> f64 {
    if yp > 1.0 - extend {
        position
    } else if yp < -1.0 + extend {
        -position
    } else {
        0.0
    }
}

pub fn linear_law(yp: f64, root: f64, tip: f64) -> f64 {
    let tip_factor = yp.abs();
    return root * (1.0 - tip_factor) + tip * tip_factor;
}
