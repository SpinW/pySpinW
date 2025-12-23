/// Mathematical constants used in the routines.
use faer::Scale;

use crate::C64;

// for convenience
// `Scale` is the faer scalar type used for matrix-scalar multiplication
pub static J: C64 = C64::new(0., 1.);
pub static SCALAR_J: Scale<C64> = Scale(J);

// Bohr magneton in units meV/T
pub static MU_B: f64 = 0.05788382;
