extern crate kiss3d;
use kiss3d::camera::{Camera, FirstPerson};
use kiss3d::light::Light;
use kiss3d::planar_camera::PlanarCamera;
use kiss3d::post_processing::post_processing_effect::PostProcessingEffect;
use kiss3d::renderer::Renderer;
use kiss3d::scene::SceneNode;
use kiss3d::window::{State, Window};

extern crate nalgebra;
use nalgebra::{Point3, Translation3, Vector3};

pub mod clock;
use clock::TimeKeeper;

mod color;

mod measures;
use measures::{Displacement, Time};

mod simulation;
use simulation::{Body, SolarSystem};

use std::collections::HashMap;
use std::time::Duration;
use std::{f32, f64};

// The default frustrum max Z is 1024 for ArcBall camera
const VIEW_DEPTH: f64 = 1e3;
const VIEW_DIAMETER: f64 = 5e12; // Just past the orbit of Neptune

// The minimum angular size of a
const MIN_ANGULAR_RES: f64 = 5e-4; // rad

fn choose_radius(body: Body, solar_system: &SolarSystem, scale_factor: f64) -> Displacement {
    let obj_props = solar_system.properties_of(body);
    let rad_guess = obj_props.radius();

    let max_dist = match body {
        Body::Earth => Displacement::from_m(0.),
        Body::Moon => obj_props.semimajor_axis(),
        _ => {
            let earth = solar_system.properties_of(Body::Earth);
            obj_props.semimajor_axis() + earth.semimajor_axis()
        }
    };
    let rad_min = max_dist * MIN_ANGULAR_RES.tan();
    rad_guess.max(rad_min) * scale_factor
}

struct ViewController {
    body_views: HashMap<Body, SceneNode>,
    camera: FirstPerson,
    scale_factor: f64,
    sun_radius: Displacement,
}

impl ViewController {
    fn new(window: &mut Window, solar_system: &SolarSystem) -> Self {
        let scale_factor = VIEW_DEPTH / VIEW_DIAMETER;
        let mut sun_radius = Displacement::from_m(0.);

        let mut body_views: HashMap<Body, SceneNode> = HashMap::new();
        for state in solar_system.body_states() {
            let radius = choose_radius(state.body(), solar_system, scale_factor);
            match state.body() {
                Body::Sun => sun_radius = radius,
                _ => {}
            }
            let mut body_view = window.add_sphere(radius.to_m() as f32);
            body_view.set_visible(false);
            let color = solar_system.properties_of(state.body()).color();
            body_view.set_color(
                color.red() as f32,
                color.green() as f32,
                color.blue() as f32,
            );
            body_views.insert(state.body(), body_view);
        }

        window.show();

        let mut camera = FirstPerson::new_with_frustrum(
            f32::consts::FRAC_PI_4,
            1e-2,
            1e4,
            Point3::new(1., 0., 0.),
            Point3::new(0., 0., 0.),
        );
        camera.set_up_axis_dir(Vector3::z_axis());

        Self {
            body_views,
            camera,
            scale_factor,
            sun_radius,
        }
    }

    fn update_view(&mut self, window: &mut Window, solar_system: &SolarSystem) {
        for state in solar_system.body_states() {
            if state.body() != Body::Earth {
                let obj_pos = state.position();
                let position = Translation3::new(
                    self.scale(obj_pos[0]),
                    self.scale(obj_pos[1]),
                    self.scale(obj_pos[2]),
                );
                self.body_views
                    .get_mut(&state.body())
                    .unwrap()
                    .set_local_translation(position);
            }
            let camera_pos = solar_system.position_of(Body::Earth);
            let sun_pos = solar_system.position_of(Body::Sun);

            // Light has to be outside of sun sphere by some distance for sun to
            // glow, so multiply radius by 2.
            let sun_rad_offset = self.sun_radius * 2. / self.scale_factor;

            let light_pos = sun_pos + (camera_pos - sun_pos).normalize() * sun_rad_offset.to_m();
            self.orient_camera(&camera_pos);
            self.set_light(window, &light_pos);
            self.body_views
                .get_mut(&state.body())
                .unwrap()
                .set_visible(true);
        }
    }

    fn orient_camera(&mut self, pos: &Vector3<f64>) {
        self.camera.look_at(
            Point3::new(self.scale(pos[0]), self.scale(pos[1]), self.scale(pos[2])),
            Point3::new(0., 0., 0.),
        )
    }

    fn set_light(&self, window: &mut Window, pos: &Vector3<f64>) {
        let light_pos = Light::Absolute(Point3::new(
            self.scale(pos[0]),
            self.scale(pos[1]),
            self.scale(pos[2]),
        ));
        window.set_light(light_pos);
    }

    fn scale(&self, length: f64) -> f32 {
        (self.scale_factor * length) as f32
    }
}

const DT: f64 = 3600.;
const FRAME_RATE: f64 = 60.;

pub struct Simulator<T: TimeKeeper> {
    solar_system: SolarSystem,
    controller: ViewController,
    clock: T,
}

impl<T: TimeKeeper + Clone> Simulator<T> {
    pub fn new(clock: T, window: &mut Window) -> Self {
        let start_time = Time::from_day(2_459_945.5); // 2023-01-01T00:00:00 UTC
        let solar_system = SolarSystem::init(start_time);
        let controller = ViewController::new(window, &solar_system);

        Self {
            solar_system,
            controller,
            clock,
        }
    }
}

impl<T: TimeKeeper + 'static> State for Simulator<T> {
    fn step(&mut self, window: &mut Window) {
        self.clock
            .advance(&Duration::from_secs_f64(1. / FRAME_RATE));
        self.solar_system.advance_time(DT);
        self.controller.update_view(window, &self.solar_system);
    }

    fn cameras_and_effect_and_renderer(
        &mut self,
    ) -> (
        Option<&mut dyn Camera>,
        Option<&mut dyn PlanarCamera>,
        Option<&mut dyn Renderer>,
        Option<&mut dyn PostProcessingEffect>,
    ) {
        (Some(&mut self.controller.camera), None, None, None)
    }
}
