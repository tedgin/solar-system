use core::f64;
use std::{
    cmp::Eq,
    collections::{HashMap, HashSet},
    fmt::{Debug, Formatter, Result},
};

use nalgebra::{Rotation3, Vector3};

use strum::VariantArray;

use crate::uom_wrapper::{
    fmt::DisplayStyle::Description,
    G,
    LuminousFlux,
    rem_euclid,
    si::{
        angle::{degree, radian},
        f64::{Angle, Length, LuminousIntensity, Mass, SolidAngle, Time},
        length::{gigameter, kilometer, meter},
        luminous_intensity::candela,
        mass::kilogram,
        solid_angle::steradian,
        time::{day, second},
    },
};

use crate::kepler_orbit as kepler;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, VariantArray)]
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
}

// All property values are correct as of 2023/01/01.
// J2000 reference frame for Sun and planets and ECI for Moon both using the
// ecliptic plane.
impl BodyProperties {
    fn sun() -> Self {
        Self {
            luminosity: LuminousIntensity::new::<candela>(3.6E28)
                * SolidAngle::new::<steradian>(1.),
            mass: Mass::new::<kilogram>(1.988_5e30),
            radius: Length::new::<kilometer>(695_700.),
            ..Default::default()
        }
    }

    fn earth() -> Self {
        Self {
            mass: Mass::new::<kilogram>(5.972_17e24),
            radius: Length::new::<kilometer>(6_371.0),
            primary: Some(Box::new(Self::sun())),
            eccentricity: 0.016_708_6,
            semimajor_axis: Length::new::<kilometer>(149_598_023.),
            inclination: Angle::new::<degree>(0.000_05),
            ascending_node: rem_euclid(Angle::new::<degree>(-11.260_64), Angle::FULL_TURN),
            periapsis_argument: Angle::new::<degree>(114.207_83),
            periapsis_time: Time::new::<day>(2_459_947.368_234_879_337),
            ..Default::default()
        }
    }

    fn moon() -> Self {
        Self {
            mass: Mass::new::<kilogram>(7.342e22),
            radius: Length::new::<kilometer>(1_737.4),
            primary: Some(Box::new(Self::earth())),
            eccentricity: 0.054_9,
            semimajor_axis: Length::new::<kilometer>(384_399.),
            inclination: Angle::new::<degree>(5.145),
            ascending_node: Angle::new::<degree>(101.502_922_218_058_2),
            periapsis_argument: Angle::new::<degree>(323.885_283_505_282_2),
            periapsis_time: Time::new::<day>(2_459_912.416_812_194_511),
            ..Default::default()
        }
    }

    fn jupiter() -> Self {
        Self {
            mass: Mass::new::<kilogram>(1.898_2e27),
            radius: Length::new::<kilometer>(69_911.),
            primary: Some(Box::new(Self::sun())),
            eccentricity: 0.048_9,
            semimajor_axis: Length::new::<gigameter>(778.479),
            inclination: Angle::new::<degree>(1.303),
            ascending_node: Angle::new::<degree>(100.464),
            periapsis_argument: Angle::new::<degree>(273.867),
            periapsis_time: Time::new::<day>(2_459_751.897_397_325_840),
            ..Default::default()
        }
    }

    fn mars() -> Self {
        Self {
            mass: Mass::new::<kilogram>(6.417_1e23),
            radius: Length::new::<kilometer>(3_389.5),
            primary: Some(Box::new(Self::sun())),
            eccentricity: 0.093_4,
            semimajor_axis: Length::new::<kilometer>(227_939_366.),
            inclination: Angle::new::<degree>(1.850),
            ascending_node: Angle::new::<degree>(49.578_54),
            periapsis_argument: Angle::new::<degree>(286.5),
            periapsis_time: Time::new::<day>(2_459_751.897_397_325_840),
            ..Default::default()
        }
    }

    fn mercury() -> Self {
        Self {
            mass: Mass::new::<kilogram>(3.301_1e23),
            radius: Length::new::<kilometer>(2_439.7),
            primary: Some(Box::new(Self::sun())),
            eccentricity: 0.205_630,
            semimajor_axis: Length::new::<kilometer>(57_909_050.),
            inclination: Angle::new::<degree>(7.005),
            ascending_node: Angle::new::<degree>(48.331),
            periapsis_argument: Angle::new::<degree>(29.124),
            periapsis_time: Time::new::<day>(2_459_947.345_508_896_280),
            ..Default::default()
        }
    }

    fn neptune() -> Self {
        Self {
            mass: Mass::new::<kilogram>(1.024_13e26),
            radius: Length::new::<kilometer>(24_622.),
            primary: Some(Box::new(Self::sun())),
            eccentricity: 0.008_678,
            semimajor_axis: Length::new::<kilometer>(4.50e9),
            inclination: Angle::new::<degree>(1.770),
            ascending_node: Angle::new::<degree>(131.783),
            periapsis_argument: Angle::new::<degree>(273.187),
            periapsis_time: Time::new::<day>(2_464_955.570_929_014_124),
            ..Default::default()
        }
    }

    fn saturn() -> Self {
        Self {
            mass: Mass::new::<kilogram>(5.683_4e26),
            radius: Length::new::<kilometer>(58_232.),
            primary: Some(Box::new(Self::sun())),
            eccentricity: 0.056_5,
            semimajor_axis: Length::new::<kilometer>(1_433.53e6),
            inclination: Angle::new::<degree>(2.485),
            ascending_node: Angle::new::<degree>(113.665),
            periapsis_argument: Angle::new::<degree>(339.392),
            periapsis_time: Time::new::<day>(2_459_751.897_397_325_840),
            ..Default::default()
        }
    }

    fn uranus() -> Self {
        Self {
            mass: Mass::new::<kilogram>(8.681_0e25),
            radius: Length::new::<kilometer>(25_362.),
            primary: Some(Box::new(Self::sun())),
            eccentricity: 0.047_17,
            semimajor_axis: Length::new::<gigameter>(2_870.972),
            inclination: Angle::new::<degree>(0.773),
            ascending_node: Angle::new::<degree>(74.006),
            periapsis_argument: Angle::new::<degree>(96.998_857),
            periapsis_time: Time::new::<day>(2_469_819.223_219_580_948),
            ..Default::default()
        }
    }

    fn venus() -> Self {
        Self {
            mass: Mass::new::<kilogram>(4.867_5e24),
            radius: Length::new::<kilometer>(6_051.8),
            primary: Some(Box::new(Self::sun())),
            eccentricity: 0.006_772,
            semimajor_axis: Length::new::<kilometer>(108_208_000.),
            inclination: Angle::new::<degree>(3.394_58),
            ascending_node: Angle::new::<degree>(76.680),
            periapsis_argument: Angle::new::<degree>(54.884),
            periapsis_time: Time::new::<day>(2_460_051.982_227_367_815),
            ..Default::default()
        }
    }

    pub fn properties_for(body: Body) -> Self {
        match body {
            Body::Sun => Self::sun(),
            Body::Mercury => Self::mercury(),
            Body::Venus => Self::venus(),
            Body::Earth => Self::earth(),
            Body::Moon => Self::moon(),
            Body::Mars => Self::mars(),
            Body::Jupiter => Self::jupiter(),
            Body::Saturn => Self::saturn(),
            Body::Uranus => Self::uranus(),
            Body::Neptune => Self::neptune(),
        }
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
            Some(_) => kepler::apsis(self.eccentricity, self.semimajor_axis),
        }
    }

    pub fn eccentric_anomaly(&self, jd: Time) -> Angle {
        match &self.primary {
            None => Angle::new::<radian>(f64::NAN),
            Some(primary) => {
                let t = kepler::period(primary.mass, self.mass, self.semimajor_axis);
                kepler::eccentric_anomaly(
                    self.eccentricity, kepler::mean_anomaly(t, self.periapsis_time, jd),
                )
            }
        }
    }

    // See https://en.wikipedia.org/wiki/Orbital_elements
    pub fn orbit_to_ecliptic(&self, vector: &Vector3<f64>) -> Vector3<f64> {
        match &self.primary {
            None => *vector,
            Some(_) => {
                let lon_rot = Rotation3::from_axis_angle(
                    &Vector3::z_axis(), -self.ascending_node.get::<radian>()
                );
                let inc_rot = Rotation3::from_axis_angle(
                    &Vector3::x_axis(), -self.inclination.get::<radian>()
                );
                let orb_rot = Rotation3::from_axis_angle(
                    &Vector3::z_axis(), -self.periapsis_argument.get::<radian>());
                (lon_rot * inc_rot * orb_rot).transform_vector(vector)
            },
        }
    }

    pub fn true_anomaly(&self, jd: Time) -> Angle {
        match &self.primary {
            None => Angle::new::<radian>(f64::NAN),
            Some(_) => kepler::true_anomaly(self.eccentricity, self.eccentric_anomaly(jd)),
        }
    }

    pub fn orbital_distance(&self, jd: Time) -> Length {
        match &self.primary {
            None => Length::new::<meter>(0.),
            Some(_) => kepler::radial_distance(
                self.semimajor_axis, self.eccentricity, self.eccentric_anomaly(jd),
            ),
        }
    }

    pub fn orbital_position(&self, jd: Time) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(_) => {
                let pos_2 = kepler::position_m(self.orbital_distance(jd), self.true_anomaly(jd));
                Vector3::new(pos_2[0], pos_2[1], 0.)
            }
        }
    }

    pub fn primary_ecliptic_position(&self, jd: Time) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(_) => self.orbit_to_ecliptic(&self.orbital_position(jd)),
        }
    }

    pub fn sun_ecliptic_position(&self, jd: Time) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(primary) => primary.sun_ecliptic_position(jd) + self.primary_ecliptic_position(jd),
        }
    }

    pub fn orbital_velocity(&self, jd: Time) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(primary) => {
                let orbital_distance = self.orbital_distance(jd);
                let speed = kepler::speed(
                    primary.mass,
                    self.mass,
                    self.semimajor_axis,
                    orbital_distance,
                );
                let pos_2 = kepler::position_m(orbital_distance, self.true_anomaly(jd));
                let vel_2 = kepler::velocity_mps(self.eccentricity, speed, &pos_2);
                Vector3::new(vel_2[0], vel_2[1], 0.)
            }
        }
    }

    pub fn primary_ecliptic_velocity(&self, jd: Time) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(_) => self.orbit_to_ecliptic(&self.orbital_velocity(jd)),
        }
    }

    pub fn sun_ecliptic_velocity(&self, jd: Time) -> Vector3<f64> {
        match &self.primary {
            None => Vector3::zeros(),
            Some(primary) => primary.sun_ecliptic_velocity(jd) + self.primary_ecliptic_velocity(jd),
        }
    }
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
        }
    }
}


struct OrbitalState {
    mass: Mass,

    // TODO: Replace Vector3<Length> once uom PRs accepted.
    position: Vector3<f64>,

    // TODO: Replace Vector3<Length> once uom PRs accepted.
    velocity: Vector3<f64>,
}

impl OrbitalState {
    fn new(mass: Mass, position: &Vector3<f64>, velocity: &Vector3<f64>) -> Self {
        Self {
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

    fn apply_force(&mut self, force: &Vector3<f64>, dt: Time) {
        self.velocity += force * dt.get::<second>() / self.mass.get::<kilogram>();
        self.position += self.velocity * dt.get::<second>();
    }
}

impl Debug for OrbitalState {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{:?}\t{:?}", self.position, self.velocity)
    }
}

pub struct SolarSystem {
    body_properties: HashMap<Body, BodyProperties>,
    body_states: HashMap<Body, OrbitalState>,
    epoch_jd: Time,
    elapsed_time: Time,
}

impl SolarSystem {
    pub fn init(start_time: Time) -> Self {
        let mut body_properties = HashMap::new();
        let mut states = HashMap::new();

        for body in Body::VARIANTS {
            let props = BodyProperties::properties_for(*body);
            let state =  OrbitalState::new(
                props.mass,
                &props.sun_ecliptic_position(start_time),
                &props.sun_ecliptic_velocity(start_time),
            );
            states.insert(*body, state);
            body_properties.insert(*body, props);
        }

        Self {
            body_properties,
            body_states: states,
            epoch_jd: start_time,
            elapsed_time: Time::new::<second>(0.),
        }
    }

    pub fn advance_time(&mut self, dt: Time) {
        self.elapsed_time += dt;

        let mut net_forces = vec![Vector3::zeros(); Body::VARIANTS.len()];

        for i in 0..Body::VARIANTS.len() {
            let body_i = &self.body_states[&Body::VARIANTS[i]];
            for j in (i + 1)..Body::VARIANTS.len() {
                let body_j = &self.body_states[&Body::VARIANTS[j]];
                let gmm = G * body_i.mass * body_j.mass;
                let r = body_i.position - body_j.position;
                let force = -gmm.value / f64::powi(r.magnitude(), 3) * r;

                net_forces[i] += force;
                net_forces[j] -= force;
            }
        }

        for i in 0..self.body_states.len() {
            self.body_states.get_mut(&Body::VARIANTS[i]).unwrap().apply_force(&net_forces[i], dt);
        }
    }

    pub fn bodies(&self) -> HashSet<Body> {
        self.body_properties.keys().cloned().collect()
    }

    pub fn position_of(&self, body: Body) -> &Vector3<f64> {
        self.body_states.get(&body).unwrap().position()
    }

    pub fn velocity_of(&self, body: Body) -> &Vector3<f64> {
        self.body_states.get(&body).unwrap().velocity()
    }

    // Return the properties for a requested body
    pub fn properties_of(&self, body: Body) -> &BodyProperties {
        self.body_properties.get(&body).unwrap()
    }
}

impl Debug for SolarSystem {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let t_fmt = Time::format_args(second, Description);
        write!(f, "{}:\n", t_fmt.with(self.epoch_jd + self.elapsed_time))?;
        for body in Body::VARIANTS.iter() {
            write!(f, "\t{:?}\t{:?}\n", body, self.body_states[&body])?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::uom_wrapper::si::angle::revolution;
    use crate::test::assert_rel_eq;
    use super::*;

    fn epoch() -> Time {
        Time::new::<day>(2_459_945.5)
    }

    fn primary() -> BodyProperties {
        BodyProperties { ..Default::default() }
    }

    fn mk_body_props(
        inclination: Angle, ascending_node: Angle, periapsis_argument: Angle
    ) -> BodyProperties {
        BodyProperties {
            primary: Some(Box::new(primary())),
            inclination,
            ascending_node,
            periapsis_argument,
            ..Default::default()
        }
    }

    #[test]
    fn test_body_properties_eccentric_anomaly_sun() {
        assert!(BodyProperties::sun().eccentric_anomaly(epoch()).is_nan())
    }

    #[test]
    fn test_body_properties_eccentric_anomaly_not_sun() {
        assert!(!BodyProperties::earth().eccentric_anomaly(epoch()).is_nan())
    }

    #[test]
    fn test_body_properties_orbital_distance_sun() {
        assert_eq!(BodyProperties::sun().orbital_distance(epoch()).get::<meter>(), 0.)
    }

    #[test]
    fn test_body_properties_orbital_distance_not_sun() {
        assert_ne!(BodyProperties::earth().orbital_distance(epoch()).get::<meter>(), 0.)
    }

    #[test]
    fn test_body_properties_orbit_to_ecliptic_periapsis_argument() {
        let props = mk_body_props(
            Angle::new::<revolution>(0.),
            Angle::new::<revolution>(0.),
            Angle::new::<revolution>(0.25),
        );
        let act = props.orbit_to_ecliptic(&Vector3::new(1f64, 0., 0.));
        assert_rel_eq!(act, Vector3::new(0f64, -1., 0.))
    }

    #[test]
    fn test_body_properties_orbit_to_ecliptic_inclination() {
        let props = mk_body_props(
            Angle::new::<revolution>(0.5),
            Angle::new::<revolution>(0.),
            Angle::new::<revolution>(0.),
        );
        let act = props.orbit_to_ecliptic(&Vector3::new(0f64, -1., 0.));
        assert_rel_eq!(act, Vector3::new(0f64, 1., 0.))
    }

    #[test]
    fn test_body_properties_orbit_to_ecliptic_ascending_node() {
        let props = mk_body_props(
            Angle::new::<revolution>(0.),
            Angle::new::<revolution>(0.25),
            Angle::new::<revolution>(0.),
        );
        let act = props.orbit_to_ecliptic(&Vector3::new(0f64, 1., 0.));
        assert_rel_eq!(act, Vector3::new(1f64, 0., 0.))
    }

    #[test]
    fn test_body_properties_orbit_to_ecliptic_combined() {
        let props = mk_body_props(
            Angle::new::<revolution>(0.5),
            Angle::new::<revolution>(0.25),
            Angle::new::<revolution>(0.25),
        );
        let orbit = Vector3::new(1f64, 0., 0.);
        let act = props.orbit_to_ecliptic(&orbit);
        assert_rel_eq!(act, orbit)
    }

    #[test]
    fn test_body_properties_orbit_to_ecliptic_sun() {
        let vector = Vector3::new(1., 2., 3.);
        assert_eq!(BodyProperties::sun().orbit_to_ecliptic(&vector), vector)
    }

    #[test]
    fn test_body_properties_orbit_to_ecliptic_not_sun() {
        let vector = Vector3::new(1., 2., 3.);
        assert_ne!(BodyProperties::earth().orbit_to_ecliptic(&vector), vector)
    }

    #[test]
    fn test_body_properties_orbital_position_sun() {
        assert_eq!(BodyProperties::sun().orbital_position(epoch()), Vector3::zeros());
    }

    #[test]
    fn test_body_properties_orbital_position_not_sun() {
        assert_ne!(BodyProperties::earth().orbital_position(epoch()), Vector3::zeros())
    }

    #[test]
    fn test_body_properties_primary_ecliptic_position_sun() {
        assert_eq!(BodyProperties::sun().primary_ecliptic_position(epoch()), Vector3::zeros())
    }

    #[test]
    fn test_body_properties_primary_ecliptic_position_not_sun() {
        assert_ne!(BodyProperties::earth().primary_ecliptic_position(epoch()), Vector3::zeros())
    }
}
