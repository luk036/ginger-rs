use std::sync::OnceLock;

/// Number of precomputed Van der Corput sequence values
pub const VDC_TABLE_SIZE: usize = 1000;

struct LdsTables {
    vdc2: Vec<f64>,
    circle_x: Vec<f64>,
    circle_y: Vec<f64>,
    cos_pi_vdc2: Vec<f64>,
}

static TABLES: OnceLock<LdsTables> = OnceLock::new();

fn get_tables() -> &'static LdsTables {
    TABLES.get_or_init(|| {
        // Generate VdCorput base-2 sequence
        let mut vgen = lds_rs::VdCorput::new(2);
        let vdc2: Vec<f64> = (0..VDC_TABLE_SIZE).map(|_| vgen.pop()).collect();

        let two_pi = std::f64::consts::TAU;
        let mut circle_x = Vec::with_capacity(VDC_TABLE_SIZE);
        let mut circle_y = Vec::with_capacity(VDC_TABLE_SIZE);
        let mut cos_pi_vdc2 = Vec::with_capacity(VDC_TABLE_SIZE);
        for &v in &vdc2 {
            let theta = v * two_pi;
            circle_x.push(theta.cos());
            circle_y.push(theta.sin());
            cos_pi_vdc2.push((std::f64::consts::PI * v).cos());
        }

        LdsTables {
            vdc2,
            circle_x,
            circle_y,
            cos_pi_vdc2,
        }
    })
}

/// Access the precomputed Van der Corput base-2 value at the given index
#[inline]
pub fn vdc2_table(index: usize) -> f64 {
    get_tables().vdc2[index]
}

/// Access the precomputed Circle base-2 x-coordinate at the given index
#[inline]
pub fn circle2_table_x(index: usize) -> f64 {
    get_tables().circle_x[index]
}

/// Access the precomputed Circle base-2 y-coordinate at the given index
#[inline]
pub fn circle2_table_y(index: usize) -> f64 {
    get_tables().circle_y[index]
}

/// Access the precomputed cos(π * vdc2_table\[i\]) value at the given index
#[inline]
pub fn cos_pi_vdc2(index: usize) -> f64 {
    get_tables().cos_pi_vdc2[index]
}
