extern crate bevy;

use bevy::{
    core_pipeline::{bloom::BloomSettings, tonemapping::Tonemapping},
    ecs::component::{ComponentHooks, StorageType},
    prelude::*,
    utils::HashMap,
    window::{PrimaryWindow, WindowMode},
};

extern crate bevy_mod_billboard;
use bevy_mod_billboard::prelude::*;

extern crate bevy_framepace;
use bevy_framepace::FramepacePlugin;

extern crate uom;
use uom::si::{
    angle::radian,
    f32::{Angle, Length},
    f64,
    length::{astronomical_unit, inch, meter},
    time::day,
};

mod simulation;
use simulation::{Body, SolarSystem};

// TODO: Create astronomical_unit_per_day velocity unit.
const MPS_TO_AUPD: f64 = 86400. / 1.495_979_E11;

// The radius of the rendering volume in AU.
const WORLD_RADIUS_AU: f32 = 100.;

// Twice the minimum angular resolution in radians of the human eye
const EYE_ANG_RES_RAD: f32 = 3e-4 * 2.;

// The minimum distance in meters someone typically sits away from a laptop display.
const MIN_OBSERVER_DIST_M: f32 = 0.3;

// The standard DPI of a monitor
const STANDARD_DPI: f32 = 96.;

// The minimum distance in AU away from the camera for an object to be rendered
const ZNEAR_AU: f32 = 0.001;

// The maximum distance in AU away from the camera for an object to be rendered
const ZFAR_AU: f32 = 100.;

// The scaling to prevent the Sun's light from saturating the camera and causing distortions
const LUMINOSITY_SCALE: f32 = 1e-26;

// The offset of the label below the body in normalized device units
const LABEL_OFFSET: f32 = 0.03;

// The scaling applied to the labels to get the to an appropriate size.
const LABEL_SCALE: f32 = 0.0003;

// Manages the visual display properties of a body
struct BodyVisual {
    name: String,
    color: Color,
}

impl BodyVisual {
    pub fn new(name: &str, color: &Color) -> Self {
        Self {
            name: name.to_string(),
            color: *color,
        }
    }

    pub fn name(&self) -> &String {
        &self.name
    }

    pub fn color(&self) -> &Color {
        &self.color
    }
}

#[derive(Resource)]
struct Simulation {
    solar_system: SolarSystem,
    body_visuals: HashMap<Body, BodyVisual>,
}

// This provides an interface to the solar system model. It ensures all the data types match those
// expected by bevy, and it ensures that the coordinate system and units are consistent with the
// World.
impl Simulation {
    // The simulation time step size in seconds
    const DT: f64 = 1800.; // half an hour

    // The Julian Date when the simulation begins (2023-01-01T00:00:00 UTC)
    const EPOCH_JD: f64 = 2_459_945.5;

    pub fn init() -> Self {
        let mut visuals = HashMap::new();

        // Color has be scaled by 10 to take advantage of HDR and bloom effects
        let sun_color = Color::srgb(9.922, 9.843, 8.275);

        let mercury_color = Color::srgb_u8(0x1a, 0x1a, 0x1a);
        let venus_color = Color::srgb_u8(0xe6, 0xe6, 0xe6);
        let earth_color = Color::srgb_u8(0x2f, 0x6a, 0x69);
        let moon_color = Color::srgb_u8(96, 86, 74);
        let mars_color = Color::srgb_u8(0x99, 0x3d, 0x00);
        let jupiter_color = Color::srgb_u8(0xb0, 0x7f, 0x35);
        let saturn_color = Color::srgb_u8(0xb0, 0x8f, 0x36);
        let uranus_color = Color::srgb_u8(0x55, 0x80, 0xaa);
        let neptune_color = Color::srgb_u8(0x36, 0x68, 0x96);
        visuals.insert(Body::Sun, BodyVisual::new("Sun", &sun_color));
        visuals.insert(Body::Mercury, BodyVisual::new("Mercury", &mercury_color));
        visuals.insert(Body::Venus, BodyVisual::new("Venus", &venus_color));
        visuals.insert(Body::Earth, BodyVisual::new("Earth", &earth_color));
        visuals.insert(Body::Moon, BodyVisual::new("Moon", &moon_color));
        visuals.insert(Body::Mars, BodyVisual::new("Mars", &mars_color));
        visuals.insert(Body::Jupiter, BodyVisual::new("Jupiter", &jupiter_color));
        visuals.insert(Body::Saturn, BodyVisual::new("Saturn", &saturn_color));
        visuals.insert(Body::Uranus, BodyVisual::new("Uranus", &uranus_color));
        visuals.insert(Body::Neptune, BodyVisual::new("Neptune", &neptune_color));
        Self {
            solar_system: SolarSystem::init(f64::Time::new::<day>(Self::EPOCH_JD)),
            body_visuals: visuals,
        }
    }

    pub fn advance(&mut self) {
        self.solar_system.advance_time(Self::DT);
    }

    pub fn apsis_of(&self, body: Body) -> f32 {
        self.solar_system
            .properties_of(body)
            .apsis()
            .get::<astronomical_unit>() as f32
    }

    pub fn color_of(&self, body: Body) -> &Color {
        self.body_visuals.get(&body).unwrap().color()
    }

    pub fn luminosity_of(&self, body: Body) -> f32 {
        self.solar_system.properties_of(body).luminosity().value as f32
    }

    pub fn name_of(&self, body: Body) -> &String {
        self.body_visuals.get(&body).unwrap().name()
    }

    pub fn position_of(&self, body: Body) -> Vec3 {
        let pos = self.solar_system.position_of(body);
        Vec3::new(
            f64::Length::new::<meter>(pos.x).get::<astronomical_unit>() as f32,
            f64::Length::new::<meter>(pos.y).get::<astronomical_unit>() as f32,
            f64::Length::new::<meter>(pos.z).get::<astronomical_unit>() as f32,
        )
    }

    pub fn radius_of(&self, body: Body) -> f32 {
        self.solar_system
            .properties_of(body)
            .radius()
            .get::<astronomical_unit>() as f32
    }

    pub fn velocity_of(&self, body: Body) -> Vec3 {
        let world_vel = (self.solar_system.velocity_of(body) * MPS_TO_AUPD).cast::<f32>();
        Vec3::new(world_vel.x, world_vel.y, world_vel.z)
    }
}

// This function advance the time by one step in the solar system model.
fn advance_sim_time(mut sim: ResMut<Simulation>) {
    sim.advance();
}

impl Component for Body {
    const STORAGE_TYPE: StorageType = StorageType::Table;

    fn register_component_hooks(_hooks: &mut ComponentHooks) {}
}

// This is the view model of a celestial body.
#[derive(Component, Default)]
struct BodyModel {
    position: Vec3,
    avatar: Option<Entity>,
    label: Option<Entity>,
}

impl BodyModel {
    pub fn new(body: Body, sim: &Simulation) -> Self {
        Self {
            position: sim.position_of(body),
            ..default()
        }
    }

    pub fn avatar(&self) -> Option<Entity> {
        self.avatar
    }

    pub fn set_avatar(&mut self, avatar: Entity) {
        self.avatar = Some(avatar);
    }

    pub fn label(&self) -> Option<Entity> {
        self.label
    }

    pub fn set_label(&mut self, label: Entity) {
        self.label = Some(label);
    }

    pub fn position(&self) -> &Vec3 {
        &self.position
    }

    pub fn update_position(&mut self, position: &Vec3) {
        self.position = *position;
    }
}

// This adds the celestial bodies being watched to the bevy World.
fn create_body_models(sim: Res<Simulation>, mut commands: Commands) {
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
        commands.spawn((body, BodyModel::new(body, &sim)));
    }
}

// This aligns the position of the bodies in the World with their positions in the solar system
// model.
fn update_bodies(sim: Res<Simulation>, mut bodies: Query<(&Body, &mut BodyModel)>) {
    for (body, mut model) in &mut bodies {
        model.update_position(&sim.position_of(*body));
    }
}

fn min_ang_res(win: &Window) -> f32 {
    let pix_size = Length::new::<inch>(1. / (STANDARD_DPI * win.scale_factor()));
    let obs_dist = Length::new::<meter>(MIN_OBSERVER_DIST_M);
    let disp_ang_res = 2. * (pix_size / (2. * obs_dist)).atan();
    Angle::new::<radian>(EYE_ANG_RES_RAD)
        .max(disp_ang_res)
        .get::<radian>()
}

fn create_avatars(
    sim: Res<Simulation>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    window: Query<&Window, With<PrimaryWindow>>,
    mut bodies: Query<(&Body, &mut BodyModel)>,
) {
    let min_ang = min_ang_res(window.single());

    for (body, mut model) in &mut bodies {
        // Camera is on Earth
        let max_dist = match body {
            Body::Moon => sim.apsis_of(*body),
            _ => sim.apsis_of(*body) + sim.apsis_of(Body::Earth),
        };

        let min_radius = max_dist * min_ang.tan() / 2.;
        let avatar_radius = sim.radius_of(*body).max(min_radius);
        let avatar_color = sim.color_of(*body);
        let avatar_lum = sim.luminosity_of(*body) * LUMINOSITY_SCALE;
        let mut avatar = commands.spawn(PbrBundle {
            mesh: meshes.add(Sphere::new(avatar_radius)),
            material: materials.add(if avatar_lum > 0. {
                StandardMaterial {
                    emissive: (*avatar_color).into(),
                    ..default()
                }
            } else {
                StandardMaterial {
                    base_color: *avatar_color,
                    ..default()
                }
            }),
            transform: Transform::from_translation(*model.position()),
            ..default()
        });
        if avatar_lum > 0. {
            avatar.with_children(|parent| {
                parent.spawn(PointLightBundle {
                    point_light: PointLight {
                        color: *avatar_color,
                        intensity: avatar_lum,
                        range: WORLD_RADIUS_AU,
                        radius: avatar_radius,
                        shadows_enabled: true,
                        ..default()
                    },
                    ..default()
                });
            });
        }
        model.set_avatar(avatar.id());
    }
}

fn update_avatars(bodies: Query<&BodyModel, With<Body>>, mut transforms: Query<&mut Transform>) {
    for model in &bodies {
        let avatar = model.avatar().unwrap();
        *transforms.get_mut(avatar).unwrap() = Transform::from_translation(*model.position());
    }
}

fn mk_lbl_transform(
    sim: &Simulation,
    model: &BodyModel,
    cam: &Camera,
    cam_trans: &GlobalTransform,
) -> Transform {
    let avatar_ndc = cam.world_to_ndc(cam_trans, *model.position()).unwrap();
    let lbl_ndc = avatar_ndc + Vec3::new(0., -LABEL_OFFSET, 0.);
    let lbl_pos = cam.ndc_to_world(cam_trans, lbl_ndc).unwrap();
    let lbl_scale = LABEL_SCALE * model.position().distance(sim.position_of(Body::Earth));
    Transform::from_translation(lbl_pos).with_scale(Vec3::splat(lbl_scale))
}

fn create_labels(
    sim: Res<Simulation>,
    mut commands: Commands,
    mut bodies: Query<(&Body, &mut BodyModel)>,
    cam: Query<(&Camera, &GlobalTransform)>,
) {
    let (cam, cam_trans) = cam.single();
    for (body, mut model) in &mut bodies {
        let lbl = commands.spawn(BillboardTextBundle {
            text: Text::from_section(
                sim.name_of(*body),
                TextStyle {
                    color: sim.color_of(*body).with_luminance(1.),
                    ..default()
                },
            ),
            transform: mk_lbl_transform(&sim, &model, cam, cam_trans),
            ..default()
        });
        model.set_label(lbl.id());
    }
}

fn update_labels(
    sim: Res<Simulation>,
    bodies: Query<&BodyModel, With<Body>>,
    cam: Query<(&Camera, &GlobalTransform)>,
    mut transforms: Query<&mut Transform>,
) {
    let (cam, cam_trans) = cam.single();
    for model in &bodies {
        *transforms.get_mut(model.label().unwrap()).unwrap() =
            mk_lbl_transform(&sim, model, cam, cam_trans);
    }
}

fn mk_cam_transform(sim: &Simulation) -> Transform {
    let cam_pos = sim.position_of(Body::Earth);
    let focus_pos = sim.position_of(Body::Sun);
    let z = focus_pos - cam_pos;
    let x = -sim.velocity_of(Body::Earth);
    let y = z.cross(x);
    Transform::from_translation(cam_pos).looking_at(focus_pos, y)
}

fn create_camera(sim: Res<Simulation>, mut commands: Commands) {
    commands.spawn((
        Camera3dBundle {
            camera: Camera {
                hdr: true,
                ..default()
            },
            projection: Projection::Perspective(PerspectiveProjection {
                near: ZNEAR_AU,
                far: ZFAR_AU,
                ..default()
            }),
            tonemapping: Tonemapping::TonyMcMapface,
            transform: mk_cam_transform(&sim),
            ..default()
        },
        BloomSettings::NATURAL,
    ));
}

fn update_camera(sim: Res<Simulation>, mut cam: Query<&mut Transform, With<Camera>>) {
    *cam.single_mut() = mk_cam_transform(&sim);
}

fn quit(input: Res<ButtonInput<KeyCode>>, mut app_exit_events: ResMut<Events<AppExit>>) {
    if input.any_pressed([KeyCode::ControlLeft, KeyCode::ControlRight])
        && input.just_released(KeyCode::KeyC)
    {
        app_exit_events.send(AppExit::Success);
    }
}

/// Run the simulation.
pub fn run() {
    let win_plug = WindowPlugin {
        primary_window: Some(Window {
            mode: WindowMode::BorderlessFullscreen,
            ..default()
        }),
        ..default()
    };
    App::new()
        .add_plugins((
            DefaultPlugins.set(win_plug),
            BillboardPlugin,
            FramepacePlugin,
        ))
        .insert_resource(Simulation::init())
        .insert_resource(ClearColor(Color::BLACK))
        .add_systems(
            Startup,
            (
                (create_body_models, create_camera),
                (create_avatars, create_labels),
            )
                .chain(),
        )
        .add_systems(
            FixedUpdate,
            (
                quit,
                (
                    advance_sim_time,
                    (update_bodies, update_camera),
                    (update_avatars, update_labels),
                )
                    .chain(),
            ),
        )
        .run();
}
