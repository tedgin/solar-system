extern crate bevy;
use bevy::{prelude::*, window::{Window, WindowPlugin}};
#[cfg(not(target_family = "wasm"))]
use bevy::window::WindowMode;

use solar_system;

#[cfg(target_family = "wasm")]
fn mk_window() -> Window {
    Window {
        canvas: Some("#canvas".into()),
        fit_canvas_to_parent: true,
        ..default()
    }
}
#[cfg(not(target_family = "wasm"))]
fn mk_window() -> Window {
    Window {
        mode: WindowMode::BorderlessFullscreen,
        ..default()
    }
}


#[cfg(not(target_family = "wasm"))]
fn quit(input: Res<ButtonInput<KeyCode>>, mut app_exit_events: ResMut<Events<AppExit>>) {
    if input.any_pressed([KeyCode::ControlLeft, KeyCode::ControlRight])
        && input.just_released(KeyCode::KeyC)
    {
        app_exit_events.send(AppExit::Success);
    }
}


pub fn main() {
    let win_plug = WindowPlugin {
        primary_window: Some(mk_window()),
        ..default()
    };

    let mut app = App::new();
    app.add_plugins(DefaultPlugins.set(win_plug));

    #[cfg(not(target_family = "wasm"))]
    app.add_systems(FixedUpdate, quit);

    solar_system::setup(&mut app).run();
}
