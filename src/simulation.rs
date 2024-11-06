use core::f64;
use std::{
    cmp::Eq,
    collections::HashMap,
    fmt::{Debug, Formatter, Result},
    marker::PhantomData,
};

use nalgebra::{Rotation3, Vector3};

use uom::{
    si::{
        angle::{degree, radian},
        f64::{Angle, Length, LuminousIntensity, Mass, SolidAngle, Time},
        length::{gigameter, kilometer, meter},
        luminous_intensity::candela,
        mass::kilogram,
        solid_angle::steradian,
        time::{day, second},
        Quantity, ISQ, SI,
    },
    typenum::{N1, N2, P1, P3, Z0},
};

// TODO: Once https://github.com/iliekturtles/uom/pull/494 has been merged and released, switch to
// using Angle::rem_euclid()
fn rem_euclid(angle: Angle, modulus: Angle) -> Angle {
    Angle::new::<radian>(angle.value.rem_euclid(modulus.value))
}

// TODO: propose this as a new derived quantity to uom
pub type LuminousFlux = Quantity<ISQ<Z0, Z0, Z0, Z0, Z0, Z0, P1>, SI<f64>, f64>;

const G: Quantity<ISQ<P3, N1, N2, Z0, Z0, Z0, Z0>, SI<f64>, f64> = Quantity {
    dimension: PhantomData,
    units: PhantomData,
    value: 6.67430e-11, // m³/kg/s²
};

const TOL: f64 = 1e-8;

#[cfg(test)]
mod asserts {
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
                    use crate::simulation::TOL;

                    if Sizeable::abs(&(act_val - exp_val)) > TOL {
                        panic!(
                            "{:?} != {:?} within absolute tolerance {:?}",
                            act_val, exp_val, TOL
                        )
                    }
                }
            }
        };
    }
    use nalgebra::{ArrayStorage, Const, Vector};

    pub use crate::assert_abs_eq;

    #[macro_export]
    macro_rules! assert_rel_eq {
        ($actual:expr, $expected:expr $(,)?) => {
            match (&$actual, &$expected) {
                (act_val, exp_val) => {
                    use crate::simulation::TOL;

                    let scaled_tol = f64::max(f64::MIN_POSITIVE, TOL * Sizeable::abs(exp_val));

                    if Sizeable::abs(&(act_val - exp_val)) > scaled_tol {
                        panic!(
                            "{:?} != {:?} within relative tolerance {:?}",
                            act_val, exp_val, TOL
                        )
                    }
                }
            }
        };
    }
    pub use crate::assert_rel_eq;
}


mod kepler_orbit {
    use core::f64;
    use std::f64::consts;

    use nalgebra::Vector2;
    use uom::{
        si::{
            angle::{radian, revolution},
            f64::{Angle, Length, Mass, Time, Velocity},
            length::meter,
            ratio::ratio,
            velocity::meter_per_second,
        },
        typenum::P3,
    };

    pub fn apsis(eccentricity: f64, semimajor_axis: Length) -> Length {
        (1. + eccentricity) * semimajor_axis
    }

    pub fn period(primary_mass: Mass, satellite_mass: Mass, semimajor_axis: Length) -> Time {
        let tot_mass = primary_mass + satellite_mass;
        consts::TAU * (semimajor_axis.powi(P3::new()) / (super::G * tot_mass)).sqrt()
    }

    pub fn mean_anomaly(period: Time, periapsis_time: Time, current_time: Time) -> Angle {
        super::rem_euclid(
            Angle::new::<revolution>(((current_time - periapsis_time) / period).get::<ratio>()),
            Angle::FULL_TURN,
        )
    }

    // See http://www.stargazing.net/kepler/mean.html
    pub fn eccentric_anomaly(eccentricity: f64, mean_anomaly: Angle) -> Angle {
        let tol = Angle::new::<radian>(super::TOL);
        let e = eccentricity;
        let ma = mean_anomaly;
        let mut ea = ma;
        let mut d = Angle::new::<radian>(f64::INFINITY);
        while d.abs() >= tol {
            d = ea - Angle::new::<radian>(e * ea.sin().get::<ratio>()) - ma;
            ea -= d / (1. - e * ea.cos().get::<ratio>());
        }
        super::rem_euclid(ea, Angle::FULL_TURN)
    }

    // See https://en.wikipedia.org/wiki/True_anomaly
    pub fn true_anomaly(eccentricity: f64, eccentric_anomaly: Angle) -> Angle {
        let e = eccentricity;
        let ea = eccentric_anomaly;
        let b = e / (1. + (1. - e.powi(2)).sqrt());
        super::rem_euclid(
            ea + 2. * (b * ea.sin() / (1. - b * ea.cos().get::<ratio>())).atan(),
            Angle::FULL_TURN,
        )
    }

    pub fn radial_distance(
        semimajor_axis: Length,
        eccentricity: f64,
        eccentric_anomaly: Angle,
    ) -> Length {
        semimajor_axis * (1. - eccentricity * eccentric_anomaly.cos().get::<ratio>())
    }

    pub fn speed(
        primary_mass: Mass,
        satellite_mass: Mass,
        semimajor_axis: Length,
        radial_distance: Length,
    ) -> Velocity {
        let m = primary_mass + satellite_mass;
        (super::G * m * (2. / radial_distance - 1. / semimajor_axis)).sqrt()
    }

    pub fn position(radial_distance: Length, true_anomaly: Angle) -> Vector2<f64> {
        let r = radial_distance.get::<meter>();
        let nu = true_anomaly;
        Vector2::new(r * nu.cos().get::<ratio>(), r * nu.sin().get::<ratio>())
    }

    pub fn velocity(eccentricity: f64, speed: Velocity, position: &Vector2<f64>) -> Vector2<f64> {
        let v_dir =
            Vector2::new(-position[1], (1. - eccentricity.powi(2)) * position[0]).normalize();
        v_dir.scale(speed.get::<meter_per_second>())
    }

    #[cfg(test)]
    mod tests {
        use std::f64::consts;

        use uom::si::angle::degree;
        use uom::si::f64::{Angle, Mass, Time};
        use uom::si::time::second;

        use super::*;
        use crate::simulation;
        use crate::simulation::asserts::*;
        use crate::simulation::*;

        #[test]
        fn test_period() {
            let act_per = period(
                Mass::new::<kilogram>(7. / simulation::G.value),
                Mass::new::<kilogram>(1. / simulation::G.value),
                Length::new::<meter>(2.),
            );
            assert_eq!(
                (act_per / Time::new::<second>(1.)).get::<ratio>(),
                consts::TAU
            )
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
            assert_abs_eq!(
                ma.get::<radian>(),
                ea.get::<radian>() - e * ea.sin().get::<ratio>()
            )
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
                ta.cos().get::<ratio>(),
                (ea.cos().get::<ratio>() - e) / (1. - e * ea.cos().get::<ratio>())
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
            let exp_r = a * (1. - e.powi(2));

            let ea = Angle::new::<revolution>(f64::atan2(
                f64::sqrt(1. - e.powi(2)) * ta.sin().get::<ratio>(),
                e + ta.cos().get::<ratio>(),
            ));
            let act_r = radial_distance(a, e, ea);
            assert_rel_eq!(act_r.get::<meter>(), exp_r.get::<meter>())
        }

        #[test]
        fn test_speed() {
            let act_speed = speed(
                Mass::new::<kilogram>(5. / simulation::G.value),
                Mass::new::<kilogram>(1. / simulation::G.value),
                Length::new::<meter>(3.),
                Length::new::<meter>(2.),
            );
            assert_rel_eq!(act_speed.get::<meter_per_second>(), 2.)
        }

        #[test]
        fn test_position() {
            let exp_pos = Vector2::new(f64::sqrt(3.), 1.);
            let act_pos = position(Length::new::<meter>(2.), Angle::new::<degree>(30.));
            assert_rel_eq!(act_pos, exp_pos)
        }
    }
}

// See https://en.wikipedia.org/wiki/Orbital_elements
fn orbit_to_ecliptic(
    inclination: Angle,
    ascending_node: Angle,
    periapsis_argument: Angle,
    vector: &Vector3<f64>,
) -> Vector3<f64> {
    let lon_rot = Rotation3::from_axis_angle(&Vector3::z_axis(), -ascending_node.get::<radian>());
    let inc_rot = Rotation3::from_axis_angle(&Vector3::x_axis(), -inclination.get::<radian>());
    let orb_rot =
        Rotation3::from_axis_angle(&Vector3::z_axis(), -periapsis_argument.get::<radian>());
    (lon_rot * inc_rot * orb_rot).transform_vector(vector)
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Body {
    Sun,
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Moon,
}

pub struct BodyProperties {
    luminosity: LuminousFlux,
    mass: Mass,
    radius: Length,
    primary: Option<Box<Self>>,
    eccentricity: f64,
    semimajor_axis: Length,
    inclination: Angle,
    ascending_node: Angle,
    periapsis_argument: Angle,
    periapsis_time: Time, // JD
    epoch: Time,
}

impl Default for BodyProperties {
    fn default() -> Self {
        Self {
            luminosity: LuminousFlux::default(),
            mass: Mass::default(),
            radius: Length::default(),
            primary: Option::default(),
            eccentricity: f64::NAN,
            semimajor_axis: Length::default(),
            inclination: Angle::new::<radian>(f64::NAN),
            ascending_node: Angle::new::<radian>(f64::NAN),
            periapsis_argument: Angle::new::<radian>(f64::NAN),
            periapsis_time: Time::new::<second>(f64::NAN),
            epoch: Time::default(),
        }
    }
}

// All property values are correct as of 2023/01/01.
// J2000 reference frame for Sun and planets and ECI for Moon both using the
// ecliptic plane.
impl BodyProperties {
    fn sun(epoch: Time) -> Self {
        Self {
            luminosity: LuminousIntensity::new::<candela>(3.6E28)
                * SolidAngle::new::<steradian>(1.),
            mass: Mass::new::<kilogram>(1.988_5e30),
            radius: Length::new::<kilometer>(695_700.),
            epoch,
            ..Default::default()
        }
    }

    fn earth(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(5.972_17e24),
            radius: Length::new::<kilometer>(6_371.0),
            primary: Some(Box::new(BodyProperties::sun(epoch))),
            eccentricity: 0.016_708_6,
            semimajor_axis: Length::new::<kilometer>(149_598_023.),
            inclination: Angle::new::<degree>(0.000_05),
            ascending_node: rem_euclid(Angle::new::<degree>(-11.260_64), Angle::FULL_TURN),
            periapsis_argument: Angle::new::<degree>(114.207_83),
            periapsis_time: Time::new::<day>(2_459_947.368_234_879_337),
            epoch,
            ..Default::default()
        }
    }

    fn moon(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(7.342e22),
            radius: Length::new::<kilometer>(1_737.4),
            primary: Some(Box::new(BodyProperties::earth(epoch))),
            eccentricity: 0.054_9,
            semimajor_axis: Length::new::<kilometer>(384_399.),
            inclination: Angle::new::<degree>(5.145),
            ascending_node: Angle::new::<degree>(101.502_922_218_058_2),
            periapsis_argument: Angle::new::<degree>(323.885_283_505_282_2),
            periapsis_time: Time::new::<day>(2_459_912.416_812_194_511),
            epoch,
            ..Default::default()
        }
    }

    fn jupiter(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(1.898_2e27),
            radius: Length::new::<kilometer>(69_911.),
            primary: Some(Box::new(BodyProperties::sun(epoch))),
            eccentricity: 0.048_9,
            semimajor_axis: Length::new::<gigameter>(778.479),
            inclination: Angle::new::<degree>(1.303),
            ascending_node: Angle::new::<degree>(100.464),
            periapsis_argument: Angle::new::<degree>(273.867),
            periapsis_time: Time::new::<day>(2_459_751.897_397_325_840),
            epoch,
            ..Default::default()
        }
    }

    fn mars(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(6.417_1e23),
            radius: Length::new::<kilometer>(3_389.5),
            primary: Some(Box::new(BodyProperties::sun(epoch))),
            eccentricity: 0.093_4,
            semimajor_axis: Length::new::<kilometer>(227_939_366.),
            inclination: Angle::new::<degree>(1.850),
            ascending_node: Angle::new::<degree>(49.578_54),
            periapsis_argument: Angle::new::<degree>(286.5),
            periapsis_time: Time::new::<day>(2_459_751.897_397_325_840),
            epoch,
            ..Default::default()
        }
    }

    fn mercury(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(3.301_1e23),
            radius: Length::new::<kilometer>(2_439.7),
            primary: Some(Box::new(BodyProperties::sun(epoch))),
            eccentricity: 0.205_630,
            semimajor_axis: Length::new::<kilometer>(57_909_050.),
            inclination: Angle::new::<degree>(7.005),
            ascending_node: Angle::new::<degree>(48.331),
            periapsis_argument: Angle::new::<degree>(29.124),
            periapsis_time: Time::new::<day>(2_459_947.345_508_896_280),
            epoch,
            ..Default::default()
        }
    }

    fn neptune(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(1.024_13e26),
            radius: Length::new::<kilometer>(24_622.),
            primary: Some(Box::new(BodyProperties::sun(epoch))),
            eccentricity: 0.008_678,
            semimajor_axis: Length::new::<kilometer>(4.50e9),
            inclination: Angle::new::<degree>(1.770),
            ascending_node: Angle::new::<degree>(131.783),
            periapsis_argument: Angle::new::<degree>(273.187),
            periapsis_time: Time::new::<day>(2_464_955.570_929_014_124),
            epoch,
            ..Default::default()
        }
    }

    fn saturn(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(5.683_4e26),
            radius: Length::new::<kilometer>(58_232.),
            primary: Some(Box::new(BodyProperties::sun(epoch))),
            eccentricity: 0.056_5,
            semimajor_axis: Length::new::<kilometer>(1_433.53e6),
            inclination: Angle::new::<degree>(2.485),
            ascending_node: Angle::new::<degree>(113.665),
            periapsis_argument: Angle::new::<degree>(339.392),
            periapsis_time: Time::new::<day>(2_459_751.897_397_325_840),
            epoch,
            ..Default::default()
        }
    }

    fn uranus(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(8.681_0e25),
            radius: Length::new::<kilometer>(25_362.),
            primary: Some(Box::new(BodyProperties::sun(epoch))),
            eccentricity: 0.047_17,
            semimajor_axis: Length::new::<gigameter>(2_870.972),
            inclination: Angle::new::<degree>(0.773),
            ascending_node: Angle::new::<degree>(74.006),
            periapsis_argument: Angle::new::<degree>(96.998_857),
            periapsis_time: Time::new::<day>(2_469_819.223_219_580_948),
            epoch,
            ..Default::default()
        }
    }

    fn venus(epoch: Time) -> Self {
        Self {
            mass: Mass::new::<kilogram>(4.867_5e24),
            radius: Length::new::<kilometer>(6_051.8),
            primary: Some(Box::new(BodyProperties::sun(epoch))),
            eccentricity: 0.006_772,
            semimajor_axis: Length::new::<kilometer>(108_208_000.),
            inclination: Angle::new::<degree>(3.394_58),
            ascending_node: Angle::new::<degree>(76.680),
            periapsis_argument: Angle::new::<degree>(54.884),
            periapsis_time: Time::new::<day>(2_460_051.982_227_367_815),
            epoch,
            ..Default::default()
        }
    }

    fn solar_system_bodies(epoch: Time) -> HashMap<Body, Self> {
        let mut bodies = HashMap::new();
        bodies.insert(Body::Sun, BodyProperties::sun(epoch));
        bodies.insert(Body::Mercury, BodyProperties::mercury(epoch));
        bodies.insert(Body::Venus, BodyProperties::venus(epoch));
        bodies.insert(Body::Earth, BodyProperties::earth(epoch));
        bodies.insert(Body::Mars, BodyProperties::mars(epoch));
        bodies.insert(Body::Jupiter, BodyProperties::jupiter(epoch));
        bodies.insert(Body::Saturn, BodyProperties::saturn(epoch));
        bodies.insert(Body::Uranus, BodyProperties::uranus(epoch));
        bodies.insert(Body::Neptune, BodyProperties::neptune(epoch));
        bodies.insert(Body::Moon, BodyProperties::moon(epoch));
        bodies
    }

    pub fn luminosity(&self) -> LuminousFlux {
        self.luminosity
    }

    pub fn radius(&self) -> Length {
        self.radius
    }

    pub fn apsis(&self) -> Length {
        match &self.primary {
            None => Length::new::<meter>(0.),
            Some(_) => kepler_orbit::apsis(self.eccentricity, self.semimajor_axis),
        }
    }

    fn eccentric_anomaly(&self) -> Angle {
        match &self.primary {
            None => Angle::new::<radian>(f64::NAN),
            Some(primary) => {
                let t = kepler_orbit::period(primary.mass, self.mass, self.semimajor_axis);
                kepler_orbit::eccentric_anomaly(
                    self.eccentricity,
                    kepler_orbit::mean_anomaly(t, self.periapsis_time, self.epoch),
                )
            }
        }
    }

    fn orbit_to_ecliptic(&self, vector: &Vector3<f64>) -> Vector3<f64> {
        match &self.primary {
            None => *vector,
            Some(_) => orbit_to_ecliptic(
                self.inclination,
                self.ascending_node,
                self.periapsis_argument,
                vector,
            ),
        }
    }

    fn true_anomaly(&self) -> Angle {
        match &self.primary {
            None => Angle::new::<radian>(f64::NAN),
            Some(_) => kepler_orbit::true_anomaly(self.eccentricity, self.eccentric_anomaly()),
        }
    }

    fn orbital_distance(&self) -> Length {
        match &self.primary {
            None => Length::new::<meter>(0.),
            Some(_) => kepler_orbit::radial_distance(
                self.semimajor_axis,
                self.eccentricity,
                self.eccentric_anomaly(),
            ),
        }
    }

    fn orbital_position(&self) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(_) => {
                let pos_2 = kepler_orbit::position(self.orbital_distance(), self.true_anomaly());
                Vector3::new(pos_2[0], pos_2[1], 0.)
            }
        }
    }

    fn primary_ecliptic_position(&self) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(_) => self.orbit_to_ecliptic(&self.orbital_position()),
        }
    }

    fn sun_ecliptic_position(&self, time: Time) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(primary) => primary.sun_ecliptic_position(time) + self.primary_ecliptic_position(),
        }
    }

    fn orbital_velocity(&self) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(primary) => {
                let orbital_distance = self.orbital_distance();
                let speed = kepler_orbit::speed(
                    primary.mass,
                    self.mass,
                    self.semimajor_axis,
                    orbital_distance,
                );
                let pos_2 = kepler_orbit::position(orbital_distance, self.true_anomaly());
                let vel_2 = kepler_orbit::velocity(self.eccentricity, speed, &pos_2);
                Vector3::new(vel_2[0], vel_2[1], 0.)
            }
        }
    }

    fn primary_ecliptic_velocity(&self) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(_) => self.orbit_to_ecliptic(&self.orbital_velocity()),
        }
    }

    fn sun_ecliptic_velocity(&self, time: Time) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(primary) => primary.sun_ecliptic_velocity(time) + self.primary_ecliptic_velocity(),
        }
    }
}

pub struct OrbitalState {
    body: Body,
    mass: Mass,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
}

impl OrbitalState {
    fn new(body: &Body, mass: Mass, position: &Vector3<f64>, velocity: &Vector3<f64>) -> Self {
        Self {
            body: *body,
            mass,
            position: *position,
            velocity: *velocity,
        }
    }

    pub fn position(&self) -> &Vector3<f64> {
        &self.position
    }

    pub fn velocity(&self) -> &Vector3<f64> {
        &self.velocity
    }

    fn apply_force(&mut self, force: &Vector3<f64>, dt: f64) {
        self.velocity += force * dt / self.mass.get::<kilogram>();
        self.position += self.velocity * dt;
    }
}

impl Debug for OrbitalState {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:?}\t{:?}", self.position, self.velocity)
    }
}

pub struct SolarSystem {
    bodies: Vec<Body>,
    body_properties: HashMap<Body, BodyProperties>,
    body_states: Vec<OrbitalState>,
    time: f64,
}

impl SolarSystem {
    pub fn init(start_time: Time) -> Self {
        let bodies = [
            Body::Earth,
            Body::Moon,
            Body::Sun,
            Body::Mercury,
            Body::Venus,
            Body::Mars,
            Body::Jupiter,
            Body::Saturn,
            Body::Uranus,
            Body::Neptune,
        ];
        let body_properties = BodyProperties::solar_system_bodies(start_time);

        let mut states: Vec<OrbitalState> = Vec::new();
        for body in bodies.iter() {
            let props = body_properties.get(body).unwrap();
            states.push(OrbitalState::new(
                body,
                props.mass,
                &props.sun_ecliptic_position(start_time),
                &props.sun_ecliptic_velocity(start_time),
            ));
        }

        Self {
            bodies: bodies.to_vec(),
            body_properties,
            body_states: states,
            time: 0.,
        }
    }

    pub fn advance_time(&mut self, dt: f64) {
        self.time += dt;

        let mut net_forces = vec![Vector3::zeros(); self.body_states.len()];

        for i in 0..self.body_states.len() {
            for j in (i + 1)..self.body_states.len() {
                let gmm = G * self.body_states[i].mass * self.body_states[j].mass;
                let r = self.body_states[i].position - self.body_states[j].position;
                let force = -gmm.value / f64::powi(r.magnitude(), 3) * r;

                net_forces[i] += force;
                net_forces[j] -= force;
            }
        }

        for i in 0..self.body_states.len() {
            self.body_states[i].apply_force(&net_forces[i], dt);
        }
    }

    pub fn position_of(&self, body: Body) -> &Vector3<f64> {
        for state in self.body_states.iter() {
            if state.body == body {
                return state.position();
            }
        }
        panic!("unknown celestial body");
    }

    pub fn velocity_of(&self, body: Body) -> &Vector3<f64> {
        for state in self.body_states.iter() {
            if state.body == body {
                return state.velocity();
            }
        }
        panic!("unknown celestial body");
    }

    // Return the properties for a requested body
    pub fn properties_of(&self, body: Body) -> &BodyProperties {
        self.body_properties.get(&body).unwrap()
    }
}

impl Debug for SolarSystem {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}:\n", self.time)?;
        for i in 0..self.bodies.len() {
            write!(f, "\t{:?}\t{:?}\n", self.bodies[i], self.body_states[i])?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use uom::si::angle::revolution;

    use super::asserts::*;
    use super::*;

    fn epoch() -> Time {
        Time::new::<day>(2_459_945.5)
    }

    #[test]
    fn test_orbit_to_ecliptic_periapsis_argument() {
        let orbit = Vector3::new(1f64, 0., 0.);
        let act = orbit_to_ecliptic(
            Angle::new::<revolution>(0.),
            Angle::new::<revolution>(0.),
            Angle::new::<revolution>(0.25),
            &orbit,
        );
        assert_rel_eq!(act, Vector3::new(0f64, -1., 0.))
    }

    #[test]
    fn test_orbit_to_ecliptic_inclination() {
        let orbit = Vector3::new(0f64, -1., 0.);
        let act = orbit_to_ecliptic(
            Angle::new::<revolution>(0.5),
            Angle::new::<revolution>(0.),
            Angle::new::<revolution>(0.),
            &orbit,
        );
        assert_rel_eq!(act, Vector3::new(0f64, 1., 0.))
    }

    #[test]
    fn test_orbit_to_ecliptic_ascending_node() {
        let orbit = Vector3::new(0f64, 1., 0.);
        let act = orbit_to_ecliptic(
            Angle::new::<revolution>(0.),
            Angle::new::<revolution>(0.25),
            Angle::new::<revolution>(0.),
            &orbit,
        );
        assert_rel_eq!(act, Vector3::new(1f64, 0., 0.))
    }

    #[test]
    fn test_orbit_to_ecliptic_combined() {
        let orbit = Vector3::new(1f64, 0., 0.);
        let act = orbit_to_ecliptic(
            Angle::new::<revolution>(0.5),
            Angle::new::<revolution>(0.25),
            Angle::new::<revolution>(0.25),
            &orbit,
        );
        assert_rel_eq!(act, orbit)
    }

    #[test]
    fn test_body_properties_eccentric_anomaly_sun() {
        let sun = BodyProperties::sun(epoch());
        assert!(sun.eccentric_anomaly().is_nan())
    }

    #[test]
    fn test_body_properties_eccentric_anomaly_not_sun() {
        let ss = BodyProperties::solar_system_bodies(epoch());
        let earth = ss.get(&Body::Earth).unwrap();
        assert!(!earth.eccentric_anomaly().is_nan())
    }

    #[test]
    fn test_body_properties_orbital_distance_sun() {
        let sun = BodyProperties::sun(epoch());
        assert_eq!(sun.orbital_distance().get::<meter>(), 0.)
    }

    #[test]
    fn test_body_properties_orbital_distance_not_sun() {
        let ss = BodyProperties::solar_system_bodies(epoch());
        let earth = ss.get(&Body::Earth).unwrap();
        assert_ne!(earth.orbital_distance().get::<meter>(), 0.)
    }

    #[test]
    fn test_body_properties_orbit_to_ecliptic_sun() {
        let sun = BodyProperties::sun(epoch());
        let vector = Vector3::new(1., 2., 3.);
        assert_eq!(sun.orbit_to_ecliptic(&vector), vector)
    }

    #[test]
    fn test_body_properties_orbit_to_ecliptic_not_sun() {
        let ss = BodyProperties::solar_system_bodies(epoch());
        let earth = ss.get(&Body::Earth).unwrap();
        let vector = Vector3::new(1., 2., 3.);
        assert_ne!(earth.orbit_to_ecliptic(&vector), vector)
    }

    #[test]
    fn test_body_properties_orbital_position_sun() {
        let sun = BodyProperties::sun(epoch());
        assert_eq!(sun.orbital_position(), Vector3::zeros());
    }

    #[test]
    fn test_body_properties_orbital_position_not_sun() {
        let ss = BodyProperties::solar_system_bodies(epoch());
        let earth = ss.get(&Body::Earth).unwrap();
        assert_ne!(earth.orbital_position(), Vector3::zeros())
    }

    #[test]
    fn test_body_properties_primary_ecliptic_position_sun() {
        let sun = BodyProperties::sun(epoch());
        assert_eq!(sun.primary_ecliptic_position(), Vector3::zeros())
    }

    #[test]
    fn test_body_properties_primary_ecliptic_position_not_sun() {
        let ss = BodyProperties::solar_system_bodies(epoch());
        let earth = ss.get(&Body::Earth).unwrap();
        assert_ne!(earth.primary_ecliptic_position(), Vector3::zeros())
    }
}
