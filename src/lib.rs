extern crate kiss3d;
use kiss3d::camera::{Camera, FirstPerson};
use kiss3d::conrod::{color, Color};
use kiss3d::light::Light;
use kiss3d::nalgebra::{Vector3, Vector2, Translation3, Point3, Reflection};
use kiss3d::planar_camera::PlanarCamera;
use kiss3d::post_processing::post_processing_effect::PostProcessingEffect;
use kiss3d::renderer::Renderer;
use kiss3d::scene::SceneNode;
use kiss3d::text::Font;
use kiss3d::window::{State, Window};

pub mod clock;
use clock::TimeKeeper;

mod measures;
use measures::{Angle, Displacement, Time};

mod simulation;
use simulation::{Body, BodyProperties, SolarSystem};

use std::collections::HashMap;
use std::time::Duration;
use std::{f32, f64};

// kiss3d::camera::FirstPerson default FOV
const FOV_DEFAULT: f32 = f32::consts::FRAC_PI_4;

// The relative amount the view frustrum needs to extend beyond the view space.
const FRUSTRUM_EASEMENT: f64 = 1.1;

// The relative amount the light has to be above the surface of the Sun to give
// the illusion that the Sun's avatar is glowing.
const LIGHT_EASEMENT: f64 = 2.;

const EPOCH_JD: f64 = 2_459_945.5; // Julian Date in days for 2023-01-01T00:00:00 UTC

const DT: f64 = 3600.;
const FRAME_RATE: f64 = 60.;

const PT_PER_IN: f64 = 72.;
const STANDARD_DPI: f64 = 96.;
const TEXT_SIZE_PT: f64 = 6.;
const LOGICAL_TEXT_SIZE_PX: f64 = TEXT_SIZE_PT * STANDARD_DPI / PT_PER_IN;

// The minimum angular resolution of the human eye, doubled
const EYE_ANGULAR_RES: f64 = 3e-4 * 2.; // rad
const MIN_OBSERVER_DIST: f64 = 0.3; // m

fn min_angular_res(window: &Window) -> Angle {
    let pixel_size = Displacement::from_in(1. / (STANDARD_DPI * window.scale_factor()));
    let obs_dist = Displacement::from_m(MIN_OBSERVER_DIST);
    let screen_angular_res = 2. * Angle::atan(pixel_size / (2. * obs_dist));
    Angle::from_rad(EYE_ANGULAR_RES).max(screen_angular_res)
}

fn avatar_color(body: Body) -> Color {
    match body {
        Body::Earth => color::BLUE,
        Body::Jupiter => color::rgb(1., 0.9, 0.7),
        Body::Mars => color::RED,
        Body::Mercury => color::grayscale(0.5),
        Body::Moon => color::WHITE,
        Body::Neptune => color::rgb(0.5, 0.5, 1.),
        Body::Saturn => color::rgb(0.8, 0.7, 0.5),
        Body::Sun => color::YELLOW,
        Body::Uranus => color::rgb(0.9, 0.9, 1.),
        Body::Venus => color::grayscale(0.8),
    }
}

fn avatar_label(body: Body) -> &'static str {
    match body {
        Body::Earth => "Earth",
        Body::Jupiter => "Jupiter",
        Body::Mars => "Mars",
        Body::Mercury => "Mercury",
        Body::Moon => "Moon",
        Body::Neptune => "Neptune",
        Body::Saturn => "Saturn",
        Body::Sun => "Sun",
        Body::Uranus => "Uranus",
        Body::Venus => "Venus",
    }
}

struct BodyAvatar {
    node: SceneNode,
    radius: Displacement,
    label: String,
}

impl BodyAvatar {
    fn new(
        properties: &BodyProperties,
        color: &Color,
        label: &str,
        min_radius: Displacement,
        window: &mut Window,
    ) -> BodyAvatar {
        let radius = properties.radius().max(min_radius);
        let mut node = window.add_sphere(radius.to_au() as f32);
        node.set_color(color.red(), color.green(), color.blue());

        BodyAvatar {
            node,
            radius,
            label: String::from(label),
        }
    }

    fn mk_avatars(window: &mut Window, solar_system: &SolarSystem) -> HashMap<Body, BodyAvatar> {
        let min_angle = min_angular_res(window);
        let mut avatars: HashMap<Body, BodyAvatar> = HashMap::new();
        for body in [
            Body::Sun,
            Body::Moon,
            Body::Mercury,
            Body::Venus,
            Body::Mars,
            Body::Jupiter,
            Body::Saturn,
            Body::Uranus,
            Body::Neptune,
        ] {
            let properties = solar_system.properties_of(body);
            let max_dist = match body {
                Body::Moon => properties.apsis(),
                _ => {
                    let earth = solar_system.properties_of(Body::Earth);
                    properties.apsis() + earth.apsis()
                }
            };
            let avatar = BodyAvatar::new(
                solar_system.properties_of(body),
                &avatar_color(body),
                avatar_label(body),
                max_dist * min_angle.tan() / 2.,
                window,
            );
            avatars.insert(body, avatar);
        }

        avatars
    }

    fn label(&self) -> &str {
        self.label.as_str()
    }

    fn radius(&self) -> Displacement {
        self.radius
    }

    fn update_position(&mut self, position: &Vector3<f64>) {
        self.node.set_local_translation(Translation3::from(
            (position / Displacement::M_PER_AU).cast::<f32>(),
        ));
    }
}

// The view space distance units are in AU.
pub struct Simulator<T: TimeKeeper> {
    clock: T,
    solar_system: SolarSystem,
    body_avatars: HashMap<Body, BodyAvatar>,
    camera: FirstPerson,
}

impl<T: TimeKeeper> Simulator<T> {
    pub fn new(clock: T, window: &mut Window) -> Self {
        let solar_system = SolarSystem::init(Time::from_day(EPOCH_JD));

        let earth_apsis = solar_system.properties_of(Body::Earth).apsis();
        let neptune_apsis = solar_system.properties_of(Body::Neptune).apsis();
        let view_depth = earth_apsis + neptune_apsis;

        let body_avatars = BodyAvatar::mk_avatars(window, &solar_system);

        let znear = solar_system.properties_of(Body::Moon).periapsis()
            * f64::cos(FOV_DEFAULT as f64 / 2.)
            / FRUSTRUM_EASEMENT;
        let zfar = view_depth * FRUSTRUM_EASEMENT;

        let mut camera = FirstPerson::new_with_frustrum(
            FOV_DEFAULT,
            znear.to_au() as f32,
            zfar.to_au() as f32,
            Point3::new(1., 0., 0.),
            Point3::new(0., 0., 0.),
        );
        camera.set_up_axis_dir(Vector3::z_axis());

        let mut sim = Self {
            clock,
            solar_system,
            body_avatars,
            camera,
        };

        sim.render(window);

        sim
    }

    fn render(&mut self, window: &mut Window) {
        self.position_avatars();
        self.orient_camera();
        self.position_light(window);
        self.label_bodies(window);
    }

    fn label_bodies(&mut self, window: &mut Window) {
        let screen_dims = Vector2::new(window.width() as f32, window.height() as f32);
        let screen_reflect = Reflection::new(Vector2::y_axis(), screen_dims.y / 2.);

        for body in self.body_avatars.keys() {
            let world_pos = Point3::from(
                (self.solar_system.position_of(*body) / Displacement::M_PER_AU).cast::<f32>(),
            );
            let view_pos = self.camera.view_transform().transform_point(&world_pos);

            if view_pos.z < 0. {
                let mut screen_pos: Vector2<f32> = self.camera.project(&world_pos, &screen_dims);
                screen_reflect.reflect(&mut screen_pos);

                let color = avatar_color(*body);
                // XXX - have to double the width and height of the screen. See
                //    https://github.com/sebcrozet/kiss3d/issues/98. This is fixed in PR
                //    https://github.com/sebcrozet/kiss3d/pull/319/.
                // window.draw_text(
                //     self.body_avatars.get(body).unwrap().label(),
                //     &screen_pos.into(),
                //     (LOGICAL_TEXT_SIZE_PX * window.scale_factor()) as f32,
                //     &Font::default(),
                //     &Point3::new(color.red(), color.green(), color.blue()),
                // );
                window.draw_text(
                    self.body_avatars.get(body).unwrap().label(),
                    &(2. * screen_pos).into(),
                    (2. * LOGICAL_TEXT_SIZE_PX * window.scale_factor()) as f32,
                    &Font::default(),
                    &Point3::new(color.red(), color.green(), color.blue()),
                );
                // XXX - ^^^
            }
        }
    }

    fn orient_camera(&mut self) {
        let camera_pos = *self.solar_system.position_of(Body::Earth);

        let eye_loc = Point3::from((camera_pos / Displacement::M_PER_AU).cast::<f32>());
        self.camera.look_at(eye_loc, Point3::origin());
    }

    fn position_avatars(&mut self) {
        for avatar in self.body_avatars.iter_mut() {
            avatar
                .1
                .update_position(self.solar_system.position_of(*avatar.0));
        }
    }

    fn position_light(&self, window: &mut Window) {
        // For the Sun to appear to glow, a light needs to be placed exterior to
        // the Sun in the direction of the camera. If the light is too close to
        // the Sun, it looks like it is reflecting light. As the distance of the
        // light from the Sun increases, the likelihood of a planet passing
        // between the light and the Sun increases, ruining the glowing
        // illusion.
        let camera_pos = self.camera.eye().cast::<f64>();
        let sun_pos =
            Point3::from(self.solar_system.position_of(Body::Sun) / Displacement::M_PER_AU);
        let sun_rad_offset = self.body_avatars.get(&Body::Sun).unwrap().radius() * LIGHT_EASEMENT;
        let sun_disp = (1. - sun_rad_offset.to_au()) * (sun_pos - camera_pos);
        let light_pos = camera_pos + sun_disp;
        window.set_light(Light::Absolute(light_pos.cast::<f32>()));
    }
}

impl<T: TimeKeeper + Clone + 'static> State for Simulator<T> {
    fn cameras_and_effect_and_renderer(
        &mut self,
    ) -> (
        Option<&mut dyn Camera>,
        Option<&mut dyn PlanarCamera>,
        Option<&mut dyn Renderer>,
        Option<&mut dyn PostProcessingEffect>,
    ) {
        (Some(&mut self.camera), None, None, None)
    }

    fn step(&mut self, window: &mut Window) {
        self.clock.advance(Duration::from_secs_f64(1. / FRAME_RATE));
        self.solar_system.advance_time(DT);
        self.render(window);
    }
}
