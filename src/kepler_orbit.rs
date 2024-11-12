use std::f64::consts;

use nalgebra::Vector2;
use crate::uom_wrapper::{
    G,
    rem_euclid,
    si::{
        angle::{radian, revolution},
        f64::{Angle, Length, Mass, Time, Velocity},
        length::meter,
        velocity::meter_per_second,
    },
    typenum::P3,
};

const TOL: f64 = 1e-8;

pub fn apsis(eccentricity: f64, semimajor_axis: Length) -> Length {
    (1. + eccentricity) * semimajor_axis
}

pub fn period(primary_mass: Mass, satellite_mass: Mass, semimajor_axis: Length) -> Time {
    let tot_mass = primary_mass + satellite_mass;
    consts::TAU * (semimajor_axis.powi(P3::new()) / (G * tot_mass)).sqrt()
}

pub fn mean_anomaly(period: Time, periapsis_time: Time, current_time: Time) -> Angle {
    rem_euclid(
        Angle::new::<revolution>(((current_time - periapsis_time) / period).into()),
        Angle::FULL_TURN,
    )
}

// See http://www.stargazing.net/kepler/mean.html
pub fn eccentric_anomaly(eccentricity: f64, mean_anomaly: Angle) -> Angle {
    let tol = Angle::new::<radian>(TOL);
    let e = eccentricity;
    let ma = mean_anomaly;
    let mut ea = ma;
    let mut d = Angle::new::<radian>(f64::INFINITY);
    while d.abs() >= tol {
        d = ea - Angle::new::<radian>((e * ea.sin()).into()) - ma;
        ea -= d / (1. - f64::from(e * ea.cos()));
    }
    rem_euclid(ea, Angle::FULL_TURN)
}

// See https://en.wikipedia.org/wiki/True_anomaly
pub fn true_anomaly(eccentricity: f64, eccentric_anomaly: Angle) -> Angle {
    let e = eccentricity;
    let ea = eccentric_anomaly;
    let b = e / (1. + (1. - e.powi(2)).sqrt());
    rem_euclid(
        ea + 2. * (b * ea.sin() / (1. - f64::from(b * ea.cos()))).atan(),
        Angle::FULL_TURN,
    )
}

pub fn radial_distance(
    semimajor_axis: Length,
    eccentricity: f64,
    eccentric_anomaly: Angle,
) -> Length {
    semimajor_axis * (1. - f64::from(eccentricity * eccentric_anomaly.cos()))
}

pub fn speed(
    primary_mass: Mass,
    satellite_mass: Mass,
    semimajor_axis: Length,
    radial_distance: Length,
) -> Velocity {
    let m = primary_mass + satellite_mass;
    (G * m * (2. / radial_distance - 1. / semimajor_axis)).sqrt()
}

// TODO: Replace Vector2<Length> once uom PRs accepted.
pub fn position_m(radial_distance: Length, true_anomaly: Angle) -> Vector2<f64> {
    let r = radial_distance;
    let nu = true_anomaly;
    Vector2::new((r * nu.cos()).get::<meter>(), (r * nu.sin()).get::<meter>())
}

// TODO: Replace Vector2<Velocity> once uom PRs accepted.
pub fn velocity_mps(eccentricity: f64, speed: Velocity, position: &Vector2<f64>) -> Vector2<f64> {
    let v_dir =
        Vector2::new(-position[1], (1. - eccentricity.powi(2)) * position[0]).normalize();
    v_dir.scale(speed.get::<meter_per_second>())
}

#[cfg(test)]
mod tests {
    use std::f64::consts;

    use crate::uom_wrapper::{
        G, si::{angle::degree, f64::{Angle, Mass, Time}, mass::kilogram, time::second}
    };

    use crate::test::{assert_abs_eq, assert_rel_eq};

    use super::*;

    #[test]
    fn test_period() {
        let act_per = period(
                Mass::new::<kilogram>(7. / G.value),
                Mass::new::<kilogram>(1. / G.value),
                Length::new::<meter>(2.),
        );
        assert_eq!(f64::from(act_per / Time::new::<second>(1.)), consts::TAU)
    }

    #[test]
    fn test_mean_anomaly() {
        let act_ma = mean_anomaly(
                Time::new::<second>(2.),
                Time::new::<second>(1.),
                Time::new::<second>(2.),
        );
        assert_abs_eq!(act_ma.get::<radian>(), consts::PI)
    }

    #[test]
    fn test_eccentric_anomaly_circle() {
        assert_eq!(
                eccentric_anomaly(0., Angle::new::<radian>(consts::FRAC_PI_2)).get::<radian>(),
                consts::FRAC_PI_2
        )
    }

    #[test]
    fn test_eccentric_anomaly_mean_relation() {
        let e = 0.5;
        let ma = Angle::new::<revolution>(0.25);
        let ea = eccentric_anomaly(e, ma);
        assert_abs_eq!(ma.get::<radian>(), ea.get::<radian>() - f64::from(e * ea.sin()))
    }

    #[test]
    fn test_true_anomaly_circle() {
        assert_eq!(
                true_anomaly(0., Angle::new::<radian>(consts::FRAC_PI_2)).get::<radian>(),
                consts::FRAC_PI_2
        )
    }

    #[test]
    fn test_true_anomaly_eccentric_relation() {
        let e = 0.5;
        let ea = Angle::new::<revolution>(0.25);
        let ta = true_anomaly(e, ea);
        assert_abs_eq!(
                f64::from(ta.cos()),
                (f64::from(ea.cos()) - e) / (1. - f64::from(e * ea.cos()))
        )
    }

    #[test]
    fn test_radial_distance_extremes() {
        let d = radial_distance(Length::new::<meter>(2.), 0.5, Angle::new::<revolution>(0.));
        assert_rel_eq!(d.get::<meter>(), 1.);
        let d = radial_distance(Length::new::<meter>(2.), 0.5, Angle::new::<revolution>(0.5));
        assert_rel_eq!(d.get::<meter>(), 3.);
    }

    #[test]
    fn test_radial_distance_true_anomaly_relation() {
        let a = Length::new::<meter>(10.);
        let e = 0.5f64;
        let ta = Angle::new::<revolution>(0.25);
        let exp_r = a * (1. - e.powi(2)) / (1. + f64::from(e * ta.cos()));

        let ea = Angle::new::<radian>(f64::atan2(
                (f64::sqrt(1. - e.powi(2)) * ta.sin()).into(),
                e + f64::from(ta.cos()),
        ));
        let act_r = radial_distance(a, e, ea);
        assert_rel_eq!(act_r.get::<meter>(), exp_r.get::<meter>())
    }

    #[test]
    fn test_speed() {
        let act_speed = speed(
                Mass::new::<kilogram>(5. / G.value),
                Mass::new::<kilogram>(1. / G.value),
                Length::new::<meter>(3.),
                Length::new::<meter>(2.),
        );
        assert_rel_eq!(act_speed.get::<meter_per_second>(), 2.)
    }

    #[test]
    fn test_position() {
        let exp_pos = Vector2::new(f64::sqrt(3.), 1.);
        let act_pos = position_m(Length::new::<meter>(2.), Angle::new::<degree>(30.));
        assert_rel_eq!(act_pos, exp_pos)
    }
}
