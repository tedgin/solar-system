use crate::measures::*;
use nalgebra::{Rotation3, Vector3};
use std::cmp::Eq;
use std::collections::HashMap;
use std::fmt::{Debug, Formatter, Result};
use std::rc::{Rc, Weak};

const G: f64 = 6.67430e-11; // m³/kg/s²
const TOL: f64 = 1e-8;

#[cfg(test)]
mod asserts {
    use nalgebra::{ArrayStorage, Const, Vector};

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
    use nalgebra::Vector2;
    use std::f64::consts;

    use crate::measures::{Angle, Displacement, Mass, Time, Velocity};

    pub fn apsis(eccentricity: f64, semimajor_axis: Displacement) -> Displacement {
        (1. + eccentricity) * semimajor_axis
    }

    pub fn periapsis(eccentricity: f64, semimajor_axis: Displacement) -> Displacement {
        (1. - eccentricity) * semimajor_axis
    }

    pub fn period(primary_mass: Mass, satellite_mass: Mass, semimajor_axis: Displacement) -> Time {
        let m = (primary_mass + satellite_mass).to_kg();
        let a = semimajor_axis.to_m();
        Time::from_s(consts::TAU * f64::sqrt(a.powi(3) / (super::G * m)))
    }

    pub fn mean_anomaly(period: Time, periapsis_time: Time, current_time: Time) -> Angle {
        Angle::from_rot((current_time - periapsis_time) / period).reduce()
    }

    // See http://www.stargazing.net/kepler/mean.html
    pub fn eccentric_anomaly(eccentricity: f64, mean_anomaly: Angle) -> Angle {
        let e = eccentricity;
        let ma = mean_anomaly.to_rad();

        let mut ea = ma;
        let mut d = f64::INFINITY;
        while d.abs() >= super::TOL {
            d = ea - e * ea.sin() - ma;
            ea -= d / (1. - e * ea.cos());
        }

        Angle::from_rad(ea).reduce()
    }

    // See https://en.wikipedia.org/wiki/True_anomaly
    pub fn true_anomaly(eccentricity: f64, eccentric_anomaly: Angle) -> Angle {
        let e = eccentricity;
        let ea = eccentric_anomaly.to_rad();
        let b = e / (1. + f64::sqrt(1. - e.powi(2)));
        Angle::from_rad(ea + 2. * f64::atan(b * ea.sin() / (1. - b * ea.cos()))).reduce()
    }

    pub fn radial_distance(
        semimajor_axis: Displacement,
        eccentricity: f64,
        eccentric_anomaly: Angle,
    ) -> Displacement {
        semimajor_axis * (1. - eccentricity * eccentric_anomaly.cos())
    }

    pub fn speed(
        primary_mass: Mass,
        satellite_mass: Mass,
        semimajor_axis: Displacement,
        radial_distance: Displacement,
    ) -> Velocity {
        let m = (primary_mass + satellite_mass).to_kg();
        let r = radial_distance.to_m();
        let a = semimajor_axis.to_m();
        Velocity::from_m_per_s(f64::sqrt(super::G * m * (2. / r - 1. / a)))
    }

    pub fn position(radial_distance: Displacement, true_anomaly: Angle) -> Vector2<f64> {
        let r = radial_distance.to_m();
        let nu = true_anomaly;
        Vector2::new(r * nu.cos(), r * nu.sin())
    }

    pub fn velocity(eccentricity: f64, speed: Velocity, position: &Vector2<f64>) -> Vector2<f64> {
        let v_dir =
            Vector2::new(-position[1], (1. - eccentricity.powi(2)) * position[0]).normalize();
        v_dir.scale(speed.to_m_per_s())
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::simulation::asserts::*;
        use crate::simulation::*;

        use crate::measures::{Displacement, Mass};
        use std::f64::consts;

        #[test]
        fn test_period() {
            let act_per = period(
                Mass::from_kg(7. / G),
                Mass::from_kg(1. / G),
                Displacement::from_m(2.),
            );
            assert_eq!(act_per / Time::from_s(1.), consts::TAU)
        }

        #[test]
        fn test_mean_anomaly() {
            let act_ma = mean_anomaly(Time::from_s(2.), Time::from_s(1.), Time::from_s(2.));
            assert_abs_eq!(act_ma.to_rad(), consts::PI)
        }

        #[test]
        fn test_eccentric_anomaly_circle() {
            assert_eq!(
                eccentric_anomaly(0., Angle::from_rad(consts::FRAC_PI_2)).to_rad(),
                consts::FRAC_PI_2
            )
        }

        #[test]
        fn test_eccentric_anomaly_mean_relation() {
            let e = 0.5;
            let ma = Angle::from_rot(0.25);
            let ea = eccentric_anomaly(e, ma);
            assert_abs_eq!(ma.to_rad(), ea.to_rad() - e * ea.to_rad().sin())
        }

        #[test]
        fn test_true_anomaly_circle() {
            assert_eq!(
                true_anomaly(0., Angle::from_rad(consts::FRAC_PI_2)).to_rad(),
                consts::FRAC_PI_2
            )
        }

        #[test]
        fn test_true_anomaly_eccentric_relation() {
            let e = 0.5;
            let ea = Angle::from_rot(0.25);
            let ta = true_anomaly(e, ea);
            assert_abs_eq!(ta.cos(), (ea.cos() - e) / (1. - e * ea.cos()))
        }

        #[test]
        fn test_radial_distance_extremes() {
            assert_rel_eq!(
                radial_distance(Displacement::from_m(2.), 0.5, Angle::from_rot(0.)).to_m(),
                1.
            );
            assert_rel_eq!(
                radial_distance(Displacement::from_m(2.), 0.5, Angle::from_rot(0.5)).to_m(),
                3.
            );
        }

        #[test]
        fn test_radial_distance_true_anomaly_relation() {
            let a = Displacement::from_m(10.);
            let e = 0.5f64;
            let ta = Angle::from_rot(0.25);
            let exp_r = a * (1. - e.powi(2));

            let ea = Angle::from_rad(f64::atan2(
                f64::sqrt(1. - e.powi(2)) * ta.sin(),
                e + ta.cos(),
            ));
            let act_r = radial_distance(a, e, ea);
            assert_rel_eq!(act_r.to_m(), exp_r.to_m())
        }

        #[test]
        fn test_speed() {
            let act_speed = speed(
                Mass::from_kg(5. / G),
                Mass::from_kg(1. / G),
                Displacement::from_m(3.),
                Displacement::from_m(2.),
            );
            assert_rel_eq!(act_speed.to_m_per_s(), 2.)
        }

        #[test]
        fn test_position() {
            let exp_pos = Vector2::new(f64::sqrt(3.), 1.);
            let act_pos = position(Displacement::from_m(2.), Angle::from_deg(30.));
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
    let lon_rot = Rotation3::from_axis_angle(&Vector3::z_axis(), -ascending_node.to_rad());
    let inc_rot = Rotation3::from_axis_angle(&Vector3::x_axis(), -inclination.to_rad());
    let orb_rot = Rotation3::from_axis_angle(&Vector3::z_axis(), -periapsis_argument.to_rad());
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
    mass: Mass,
    radius: Displacement,
    primary: Weak<Self>,
    eccentricity: f64,
    semimajor_axis: Displacement,
    inclination: Angle,
    ascending_node: Angle,
    periapsis_argument: Angle,
    periapsis_time: Time, // JD
    epoch: Time,
}

// All property values are correct as of 2023/01/01.
// J2000 reference frame for Sun and planets and ECI for Moon both using the
// ecliptic plane.
impl BodyProperties {
    fn sun(epoch: Time) -> Self {
        Self {
            mass: Mass::from_kg(1.988_5e30),
            radius: Displacement::from_km(695_700.),
            primary: Weak::new(),
            eccentricity: f64::NAN,
            semimajor_axis: Displacement::from_m(0.),
            inclination: Angle::NAN,
            ascending_node: Angle::NAN,
            periapsis_argument: Angle::NAN,
            periapsis_time: Time::NAN,
            epoch,
        }
    }

    fn earth(epoch: Time, sun: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(5.972_17e24),
            radius: Displacement::from_km(6_371.0),
            primary: sun,
            eccentricity: 0.016_708_6,
            semimajor_axis: Displacement::from_km(149_598_023.),
            inclination: Angle::from_deg(0.000_05),
            ascending_node: Angle::from_deg(-11.260_64),
            periapsis_argument: Angle::from_deg(114.207_83),
            periapsis_time: Time::from_day(2_459_947.368_234_879_337),
            epoch,
        }
    }

    fn moon(epoch: Time, earth: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(7.342e22),
            radius: Displacement::from_km(1_737.4),
            primary: earth,
            eccentricity: 0.054_9,
            semimajor_axis: Displacement::from_km(384_399.),
            inclination: Angle::from_deg(5.145),
            ascending_node: Angle::from_deg(101.502_922_218_058_2),
            periapsis_argument: Angle::from_deg(323.885_283_505_282_2).reduce(),
            periapsis_time: Time::from_day(2_459_912.416_812_194_511),
            epoch,
        }
    }

    fn jupiter(epoch: Time, sun: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(1.898_2e27),
            radius: Displacement::from_km(69_911.),
            primary: sun,
            eccentricity: 0.048_9,
            semimajor_axis: Displacement::from_gm(778.479),
            inclination: Angle::from_deg(1.303),
            ascending_node: Angle::from_deg(100.464),
            periapsis_argument: Angle::from_deg(273.867).reduce(),
            periapsis_time: Time::from_day(2_459_751.897_397_325_840),
            epoch,
        }
    }

    fn mars(epoch: Time, sun: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(6.417_1e23),
            radius: Displacement::from_km(3_389.5),
            primary: sun,
            eccentricity: 0.093_4,
            semimajor_axis: Displacement::from_km(227_939_366.),
            inclination: Angle::from_deg(1.850),
            ascending_node: Angle::from_deg(49.578_54),
            periapsis_argument: Angle::from_deg(286.5).reduce(),
            periapsis_time: Time::from_day(2_459_751.897_397_325_840),
            epoch,
        }
    }

    fn mercury(epoch: Time, sun: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(3.301_1e23),
            radius: Displacement::from_km(2_439.7),
            primary: sun,
            eccentricity: 0.205_630,
            semimajor_axis: Displacement::from_km(57_909_050.),
            inclination: Angle::from_deg(7.005),
            ascending_node: Angle::from_deg(48.331),
            periapsis_argument: Angle::from_deg(29.124),
            periapsis_time: Time::from_day(2_459_947.345_508_896_280),
            epoch,
        }
    }

    fn neptune(epoch: Time, sun: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(1.024_13e26),
            radius: Displacement::from_km(24_622.),
            primary: sun,
            eccentricity: 0.008_678,
            semimajor_axis: Displacement::from_km(4.50e9),
            inclination: Angle::from_deg(1.770),
            ascending_node: Angle::from_deg(131.783),
            periapsis_argument: Angle::from_deg(273.187).reduce(),
            periapsis_time: Time::from_day(2_464_955.570_929_014_124),
            epoch,
        }
    }

    fn saturn(epoch: Time, sun: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(5.683_4e26),
            radius: Displacement::from_km(58_232.),
            primary: sun,
            eccentricity: 0.056_5,
            semimajor_axis: Displacement::from_km(1_433.53e6),
            inclination: Angle::from_deg(2.485),
            ascending_node: Angle::from_deg(113.665),
            periapsis_argument: Angle::from_deg(339.392).reduce(),
            periapsis_time: Time::from_day(2_459_751.897_397_325_840),
            epoch,
        }
    }

    fn uranus(epoch: Time, sun: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(8.681_0e25),
            radius: Displacement::from_km(25_362.),
            primary: sun,
            eccentricity: 0.047_17,
            semimajor_axis: Displacement::from_gm(2_870.972),
            inclination: Angle::from_deg(0.773),
            ascending_node: Angle::from_deg(74.006),
            periapsis_argument: Angle::from_deg(96.998_857),
            periapsis_time: Time::from_day(2_469_819.223_219_580_948),
            epoch,
        }
    }

    fn venus(epoch: Time, sun: Weak<Self>) -> Self {
        Self {
            mass: Mass::from_kg(4.867_5e24),
            radius: Displacement::from_km(6_051.8),
            primary: sun,
            eccentricity: 0.006_772,
            semimajor_axis: Displacement::from_km(108_208_000.),
            inclination: Angle::from_deg(3.394_58),
            ascending_node: Angle::from_deg(76.680),
            periapsis_argument: Angle::from_deg(54.884),
            periapsis_time: Time::from_day(2_460_051.982_227_367_815),
            epoch,
        }
    }

    fn solar_system_bodies(epoch: Time) -> HashMap<Body, Rc<Self>> {
        let mut bodies = HashMap::new();

        bodies.insert(Body::Sun, Rc::new(BodyProperties::sun(epoch)));

        let sun = Rc::downgrade(bodies.get(&Body::Sun).unwrap());
        bodies.insert(
            Body::Mercury,
            Rc::new(BodyProperties::mercury(epoch, sun.clone())),
        );
        bodies.insert(
            Body::Venus,
            Rc::new(BodyProperties::venus(epoch, sun.clone())),
        );
        bodies.insert(
            Body::Earth,
            Rc::new(BodyProperties::earth(epoch, sun.clone())),
        );
        bodies.insert(
            Body::Mars,
            Rc::new(BodyProperties::mars(epoch, sun.clone())),
        );
        bodies.insert(
            Body::Jupiter,
            Rc::new(BodyProperties::jupiter(epoch, sun.clone())),
        );
        bodies.insert(
            Body::Saturn,
            Rc::new(BodyProperties::saturn(epoch, sun.clone())),
        );
        bodies.insert(
            Body::Uranus,
            Rc::new(BodyProperties::uranus(epoch, sun.clone())),
        );
        bodies.insert(
            Body::Neptune,
            Rc::new(BodyProperties::neptune(epoch, sun.clone())),
        );

        let earth = Rc::downgrade(bodies.get(&Body::Earth).unwrap());
        bodies.insert(
            Body::Moon,
            Rc::new(BodyProperties::moon(epoch, earth.clone())),
        );

        bodies
    }

    pub fn radius(&self) -> Displacement {
        self.radius
    }

    pub fn apsis(&self) -> Displacement {
        match self.primary.upgrade() {
            None => Displacement::from_m(0.),
            Some(_) => kepler_orbit::apsis(self.eccentricity, self.semimajor_axis),
        }
    }

    pub fn periapsis(&self) -> Displacement {
        match self.primary.upgrade() {
            None => Displacement::from_m(0.),
            Some(_) => kepler_orbit::periapsis(self.eccentricity, self.semimajor_axis),
        }
    }

    fn eccentric_anomaly(&self) -> Angle {
        match self.primary.upgrade() {
            None => Angle::NAN,
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
        match self.primary.upgrade() {
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
        match self.primary.upgrade() {
            None => Angle::NAN,
            Some(_) => kepler_orbit::true_anomaly(self.eccentricity, self.eccentric_anomaly()),
        }
    }

    fn orbital_distance(&self) -> Displacement {
        match self.primary.upgrade() {
            None => Displacement::from_m(0.),
            Some(_) => kepler_orbit::radial_distance(
                self.semimajor_axis,
                self.eccentricity,
                self.eccentric_anomaly(),
            ),
        }
    }

    fn orbital_position(&self) -> Vector3<f64> {
        match self.primary.upgrade() {
            None => Vector3::zeros(),
            Some(_) => {
                let pos_2 = kepler_orbit::position(self.orbital_distance(), self.true_anomaly());
                Vector3::new(pos_2[0], pos_2[1], 0.)
            }
        }
    }

    fn primary_ecliptic_position(&self) -> Vector3<f64> {
        match self.primary.upgrade() {
            None => Vector3::zeros(),
            Some(_) => self.orbit_to_ecliptic(&self.orbital_position()),
        }
    }

    fn sun_ecliptic_position(&self, time: Time) -> Vector3<f64> {
        match self.primary.upgrade() {
            None => Vector3::zeros(),
            Some(primary) => primary.sun_ecliptic_position(time) + self.primary_ecliptic_position(),
        }
    }

    fn orbital_velocity(&self) -> Vector3<f64> {
        match self.primary.upgrade() {
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
        match self.primary.upgrade() {
            None => Vector3::zeros(),
            Some(_) => self.orbit_to_ecliptic(&self.orbital_velocity()),
        }
    }

    fn sun_ecliptic_velocity(&self, time: Time) -> Vector3<f64> {
        match self.primary.upgrade() {
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

    fn apply_force(&mut self, force: &Vector3<f64>, dt: f64) {
        self.velocity += force * dt / self.mass.to_kg();
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
    body_properties: HashMap<Body, Rc<BodyProperties>>,
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
                let gmm = G * self.body_states[i].mass.to_kg() * self.body_states[j].mass.to_kg();
                let r = self.body_states[i].position - self.body_states[j].position;
                let force = -gmm / f64::powi(r.magnitude(), 3) * r;

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

    // Return the properties for a requested body
    pub fn properties_of(&self, body: Body) -> &BodyProperties {
        self.body_properties.get(&body).unwrap().as_ref()
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
    use super::asserts::*;
    use super::*;

    use crate::measures::Time;

    fn epoch() -> Time {
        Time::from_day(2_459_945.5)
    }

    #[test]
    fn test_orbit_to_ecliptic_periapsis_argument() {
        let orbit = Vector3::new(1f64, 0., 0.);
        let act = orbit_to_ecliptic(
            Angle::from_rot(0.),
            Angle::from_rot(0.),
            Angle::from_rot(0.25),
            &orbit,
        );
        assert_rel_eq!(act, Vector3::new(0f64, -1., 0.))
    }

    #[test]
    fn test_orbit_to_ecliptic_inclination() {
        let orbit = Vector3::new(0f64, -1., 0.);
        let act = orbit_to_ecliptic(
            Angle::from_rot(0.5),
            Angle::from_rot(0.),
            Angle::from_rot(0.),
            &orbit,
        );
        assert_rel_eq!(act, Vector3::new(0f64, 1., 0.))
    }

    #[test]
    fn test_orbit_to_ecliptic_ascending_node() {
        let orbit = Vector3::new(0f64, 1., 0.);
        let act = orbit_to_ecliptic(
            Angle::from_rot(0.),
            Angle::from_rot(0.25),
            Angle::from_rot(0.),
            &orbit,
        );
        assert_rel_eq!(act, Vector3::new(1f64, 0., 0.))
    }

    #[test]
    fn test_orbit_to_ecliptic_combined() {
        let orbit = Vector3::new(1f64, 0., 0.);
        let act = orbit_to_ecliptic(
            Angle::from_rot(0.5),
            Angle::from_rot(0.25),
            Angle::from_rot(0.25),
            &orbit,
        );
        assert_rel_eq!(act, orbit)
    }

    #[test]
    fn test_body_properties_earth() {
        let sun = Rc::new(BodyProperties::sun(epoch()));
        let earth = BodyProperties::earth(epoch(), Rc::downgrade(&sun));
        assert!(earth.primary.ptr_eq(&Rc::downgrade(&sun)))
    }

    #[test]
    fn test_body_properties_jupiter() {
        let sun = Rc::new(BodyProperties::sun(epoch()));
        let jupiter = BodyProperties::jupiter(epoch(), Rc::downgrade(&sun));
        assert!(jupiter.primary.ptr_eq(&Rc::downgrade(&sun)))
    }

    #[test]
    fn test_body_properties_mars() {
        let sun = Rc::new(BodyProperties::sun(epoch()));
        let mars = BodyProperties::mars(epoch(), Rc::downgrade(&sun));
        assert!(mars.primary.ptr_eq(&Rc::downgrade(&sun)))
    }

    #[test]
    fn test_body_properties_mercury() {
        let sun = Rc::new(BodyProperties::sun(epoch()));
        let mercury = BodyProperties::mercury(epoch(), Rc::downgrade(&sun));
        assert!(mercury.primary.ptr_eq(&Rc::downgrade(&sun)))
    }

    #[test]
    fn test_body_properties_moon() {
        let primary = Rc::new(BodyProperties::sun(epoch()));
        let moon = BodyProperties::moon(epoch(), Rc::downgrade(&primary));
        assert!(moon.primary.ptr_eq(&Rc::downgrade(&primary)))
    }

    #[test]
    fn test_body_properties_neptune() {
        let sun = Rc::new(BodyProperties::sun(epoch()));
        let neptune = BodyProperties::neptune(epoch(), Rc::downgrade(&sun));
        assert!(neptune.primary.ptr_eq(&Rc::downgrade(&sun)))
    }

    #[test]
    fn test_body_properties_saturn() {
        let sun = Rc::new(BodyProperties::sun(epoch()));
        let saturn = BodyProperties::saturn(epoch(), Rc::downgrade(&sun));
        assert!(saturn.primary.ptr_eq(&Rc::downgrade(&sun)))
    }

    #[test]
    fn test_body_properties_sun() {
        assert!(BodyProperties::sun(epoch()).primary.upgrade().is_none())
    }

    #[test]
    fn test_body_properties_uranus() {
        let sun = Rc::new(BodyProperties::sun(epoch()));
        let uranus = BodyProperties::uranus(epoch(), Rc::downgrade(&sun));
        assert!(uranus.primary.ptr_eq(&Rc::downgrade(&sun)))
    }

    #[test]
    fn test_body_properties_venus() {
        let sun = Rc::new(BodyProperties::sun(epoch()));
        let venus = BodyProperties::venus(epoch(), Rc::downgrade(&sun));
        assert!(venus.primary.ptr_eq(&Rc::downgrade(&sun)))
    }

    #[test]
    fn test_body_properties_solar_system_entities() {
        let ss = BodyProperties::solar_system_bodies(epoch());
        let sun = ss.get(&Body::Sun).unwrap();

        let planets = [
            Body::Mercury,
            Body::Venus,
            Body::Earth,
            Body::Mars,
            Body::Jupiter,
            Body::Saturn,
            Body::Uranus,
            Body::Neptune,
        ];
        for planet in planets.iter() {
            assert!(
                ss.get(planet).unwrap().primary.ptr_eq(&Rc::downgrade(sun)),
                "{:?} should have Sun as its primary",
                planet
            );
        }

        let moon = ss.get(&Body::Moon).unwrap();
        let earth = ss.get(&Body::Earth).unwrap();
        assert!(
            moon.primary.ptr_eq(&Rc::downgrade(earth)),
            "Moon should have Earth as its primary"
        );
    }

    #[test]
    fn test_body_properties_eccentric_anomaly_sun() {
        let sun = BodyProperties::sun(epoch());
        assert!(sun.eccentric_anomaly().to_rad().is_nan())
    }

    #[test]
    fn test_body_properties_eccentric_anomaly_not_sun() {
        let ss = BodyProperties::solar_system_bodies(epoch());
        let earth = ss.get(&Body::Earth).unwrap();
        assert!(!earth.eccentric_anomaly().to_rad().is_nan())
    }

    #[test]
    fn test_body_properties_orbital_distance_sun() {
        let sun = BodyProperties::sun(epoch());
        assert_eq!(sun.orbital_distance().to_m(), 0.)
    }

    #[test]
    fn test_body_properties_orbital_distance_not_sun() {
        let ss = BodyProperties::solar_system_bodies(epoch());
        let earth = ss.get(&Body::Earth).unwrap();
        assert_ne!(earth.orbital_distance().to_m(), 0.)
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
