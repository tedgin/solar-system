// TODO (deferred) replace with an existing library from crates.io

use std::f64::consts;
use std::ops::{Add, Div, Mul, Neg, Sub, SubAssign};

#[derive(Clone, Copy, PartialEq, PartialOrd)]
pub struct Angle {
    rad: f64,
}

impl Angle {
    const DEG_PER_ROT: f64 = 360.;
    const RAD_PER_DEG: f64 = consts::TAU / Self::DEG_PER_ROT;

    pub const NAN: Self = Self { rad: f64::NAN };

    pub fn from_deg(deg: f64) -> Self {
        Self { rad: Self::RAD_PER_DEG * deg }
    }

    pub fn from_rad(rad: f64) -> Self {
        Self { rad }
    }

    pub fn from_rot(rot: f64) -> Self {
        Self { rad: consts::TAU * rot }
    }

    pub fn to_rad(self) -> f64 {
        self.rad
    }

    pub fn reduce(self) -> Self {
        let reduced = self.rad.rem_euclid(consts::TAU);
        if reduced > consts::PI {
            Self { rad: reduced - consts::TAU }
        } else {
            Self { rad: reduced }
        }
    }

    pub fn abs(&self) -> Self {
        Self { rad: self.rad.abs() }
    }

    pub fn cos(self) -> f64 {
        self.rad.cos()
    }

    pub fn sin(self) -> f64 {
        self.rad.sin()
    }
}

impl Neg for Angle {
    type Output = Self;

    fn neg(self) -> Self {
        if self.rad < consts::PI {
            Self { rad: -self.rad }
        } else {
            self
        }
    }
}

impl Add for Angle {
    type Output = Self;

    fn add(self, addend: Self) -> Self {
        Self { rad: self.rad + addend.rad }
    }
}

impl Sub for Angle {
    type Output = Self;

    fn sub(self, subtrahend: Self) -> Self {
        Self { rad: self.rad - subtrahend.rad }
    }
}

impl SubAssign for Angle {
    fn sub_assign(&mut self, subtrahend: Self) {
        self.rad -= subtrahend.rad;
    }
}

impl Mul<Angle> for f64 {
    type Output = Angle;

    fn mul(self, multiplicand: Angle) -> Angle {
        Angle { rad: self * multiplicand.rad }
    }
}

#[derive(Clone, Copy)]
pub struct Displacement {
    m: f64,
}

impl Displacement {
    const M_PER_KM: f64 = 1e3;
    const M_PER_GM: f64 = 1e9;

    pub fn from_m(m: f64) -> Self {
        Self { m }
    }

    pub fn from_km(km: f64) -> Self {
        Self { m: Self::M_PER_KM * km }
    }

    pub fn from_gm(gm: f64) -> Self {
        Self { m: Self::M_PER_GM * gm }
    }

    pub fn to_m(self) -> f64 {
        self.m
    }

    pub fn max(self, other: Self) -> Self {
        Self { m: self.m.max(other.m) }
    }
}

impl Add for Displacement {
    type Output = Self;

    fn add(self, addend: Self) -> Self {
        Self { m: self.m + addend.m }
    }
}

impl Mul<f64> for Displacement {
    type Output = Self;

    fn mul(self, multiplicand: f64) -> Self {
        Self { m: self.m * multiplicand }
    }
}

impl Div<f64> for Displacement {
    type Output = Self;

    fn div(self, divisor: f64) -> Self {
        Self { m : self.m / divisor }
    }
}

#[derive(Clone, Copy)]
pub struct Mass {
    kg: f64,
}

impl Mass {
    pub fn from_kg(kg: f64) -> Self {
        Self { kg }
    }

    pub fn to_kg(self) -> f64 {
        self.kg
    }
}

impl Add for Mass {
    type Output = Self;

    fn add(self, addend: Self) -> Self {
        Self { kg: self.kg + addend.kg }
    }
}

#[derive(Clone, Copy)]
pub struct Time {
    s: f64,
}

impl Time {
    const HR_PER_DAY: f64 = 24.;
    const MIN_PER_HR: f64 = 60.;
    const S_PER_MIN: f64 = 60.;
    const S_PER_DAY: f64 = Self::S_PER_MIN * Self::MIN_PER_HR * Self::HR_PER_DAY;

    pub const NAN: Self = Self { s: f64::NAN };

    pub fn from_day(day: f64) -> Self {
        Self { s: Self::S_PER_DAY * day }
    }

    pub fn from_s(s: f64) -> Self {
        Self { s }
    }
}

impl Sub for Time {
    type Output = Self;

    fn sub(self, subtrahend: Self) -> Self {
        Self {
            s: self.s - subtrahend.s,
        }
    }
}

impl Div for Time {
    type Output = f64;

    fn div(self, divisor: Self) -> f64 {
        self.s / divisor.s
    }
}

pub struct Velocity {
    mps: f64,
}

impl Velocity {
    pub fn from_m_per_s(mps: f64) -> Velocity {
        Self { mps }
    }

    pub fn to_m_per_s(self) -> f64 {
        self.mps
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts;

    use super::*;

    #[test]
    fn test_angle_from_deg() {
        assert_eq!(Angle::from_deg(180.).rad, consts::PI);
    }

    #[test]
    fn test_angle_from_rot() {
        assert_eq!(Angle::from_rot(1.).rad, consts::TAU);
    }

    #[test]
    fn test_angle_reduce() {
        assert_eq!(Angle::from_deg(-180.).reduce().rad, consts::PI);
        assert_eq!(Angle::from_deg(360.).reduce().rad, 0.);
        assert!(Angle::from_deg(180. * (1. + f64::EPSILON)).reduce().rad < 0.);
        assert!(Angle::from_deg(-180. * (1. - f64::EPSILON)).reduce().rad < 0.);
    }

    #[test]
    fn test_angle_abs() {
        assert_eq!(Angle::from_rad(2.).abs().rad, 2.);
        assert_eq!(Angle::from_rad(-1.).abs().rad, 1.);
    }

    #[test]
    fn test_angle_cos() {
        assert_eq!(Angle::from_rad(0.).cos(), 1.)
    }

    #[test]
    fn test_angle_sin() {
        assert_eq!(Angle::from_rad(consts::PI).sin(), consts::PI.sin())
    }

    #[test]
    fn test_angle_neg() {
        assert_eq!((-Angle::from_deg(90.)).rad, -consts::FRAC_PI_2);
        assert_eq!((-Angle::from_deg(180.)).rad, consts::PI);
    }

    #[test]
    fn test_angle_sub() {
        assert_eq!((Angle::from_rad(1.) - Angle::from_rad(-3.)).rad, 4.)
    }

    #[test]
    fn test_angle_sub_assign() {
        let mut a = Angle::from_rad(1.);
        a -= Angle::from_rad(-3.);
        assert_eq!(a.rad, 4.)
    }

    #[test]
    fn test_displacement_max() {
        assert_eq!(Displacement::from_m(-2.).max(Displacement::from_m(1.)).to_m(), 1.)
    }

    #[test]
    fn test_displacement_add() {
        assert_eq!((Displacement::from_m(1.) + Displacement::from_m(2.)).m, 3.)
    }

    #[test]
    fn test_displacement_mul() {
        assert_eq!((Displacement::from_m(2.) * 3.).m, 6.)
    }

    #[test]
    fn test_displacement_div() {
        assert_eq!((Displacement::from_m(4.) / 2.).to_m(), 2.)
    }

    #[test]
    fn test_mass_add() {
        assert_eq!((Mass::from_kg(1.) + Mass::from_kg(2.)).to_kg(), 3.)
    }

    #[test]
    fn test_time_from_day() {
        assert_eq!(Time::from_day(0.).s, 0.);
        assert_eq!(Time::from_day(2.).s, 172800.);
    }

    #[test]
    fn test_time_div() {
        assert_eq!(Time::from_s(3.) / Time::from_s(2.), 1.5);
    }

    #[test]
    fn test_time_sub() {
        assert_eq!((Time::from_s(4.) - Time::from_s(3.)).s, 1.);
    }
}
