use nalgebra::{ArrayStorage, Const, Vector};

pub const TOL: f64 = 1e-8;

pub trait Sizeable {
    fn abs(&self) -> f64;
}

impl Sizeable for f64 {
    fn abs(&self) -> f64 {
        f64::abs(*self)
    }
}

impl<const R: usize> Sizeable for Vector<f64, Const<R>, ArrayStorage<f64, R, 1>> {
    fn abs(&self) -> f64 {
        self.magnitude()
    }
}

#[macro_export]
macro_rules! assert_abs_eq {
    ($actual:expr, $expected:expr $(,)?) => {
        match (&$actual, &$expected) {
            (act_val, exp_val) => {
                if $crate::test::Sizeable::abs(&(act_val - exp_val)) > $crate::test::TOL {
                    panic!(
                        "{:?} != {:?} within absolute tolerance {:?}",
                        act_val, exp_val, $crate::test::TOL
                    )
                }
            }
        }
    };
}
pub use crate::assert_abs_eq;

#[macro_export]
macro_rules! assert_rel_eq {
    ($actual:expr, $expected:expr $(,)?) => {
        match (&$actual, &$expected) {
            (act_val, exp_val) => {
                let scaled_tol = f64::max(
                    f64::MIN_POSITIVE, $crate::test::TOL * $crate::test::Sizeable::abs(exp_val)
                );

                if $crate::test::Sizeable::abs(&(act_val - exp_val)) > scaled_tol {
                    panic!(
                        "{:?} != {:?} within relative tolerance {:?}",
                        act_val, exp_val, $crate::test::TOL
                    )
                }
            }
        }
    };
}
pub use crate::assert_rel_eq;
