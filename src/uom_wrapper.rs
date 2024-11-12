use std::marker::PhantomData;

extern crate uom;
pub use uom::*;
use si::{angle::radian, f64::Angle, ISQ, Quantity, SI};
use typenum::{N1, N2, P1, P3, Z0};

// TODO: Switch to uom's LuminousFlux once https://github.com/iliekturtles/uom/pull/313 has been
// merged and released.
pub type LuminousFlux = Quantity<ISQ<Z0, Z0, Z0, Z0, Z0, Z0, P1>, SI<f64>, f64>;

// TODO: Create astronomical_unit_per_day velocity unit to uom library
pub const MPS_TO_AUPD: f64 = 86400. / 1.495_979_E11;

pub const G: Quantity<ISQ<P3, N1, N2, Z0, Z0, Z0, Z0>, SI<f64>, f64> = Quantity {
    dimension: PhantomData,
    units: PhantomData,
    value: 6.67430e-11, // m³/kg/s²
};

// TODO: Once https://github.com/iliekturtles/uom/pull/494 has been merged and released, switch to
// using Angle::rem_euclid()
pub fn rem_euclid(angle: Angle, modulus: Angle) -> Angle {
    Angle::new::<radian>(angle.value.rem_euclid(modulus.value))
}
