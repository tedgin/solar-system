extern crate kiss3d;
use kiss3d::window::Window;

use solar_system::clock::TimeKeeper;
use solar_system::Simulator;

#[cfg(feature = "wasm")]
extern crate console_error_panic_hook;

#[cfg(feature = "wasm")]
mod wasm_clock;
#[cfg(feature = "wasm")]
use wasm_clock::Clock;

#[cfg(not(feature = "wasm"))]
mod native_clock;
#[cfg(not(feature = "wasm"))]
use native_clock::Clock;

pub fn main() {
    #[cfg(feature = "wasm")]
    console_error_panic_hook::set_once();

    let mut window = Window::new_hidden("Solar System from Earth");
    let app = Simulator::new(Clock::init(), &mut window);
    window.show();
    window.render_loop(app);
}
