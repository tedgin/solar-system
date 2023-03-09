extern crate kiss3d;
use kiss3d::window::Window;

use solar_system::Simulator;

#[cfg(feature = "wasm")]
extern crate console_error_panic_hook;

const FRAME_LIMIT: u64 = 60;

pub fn main() {
    #[cfg(feature = "wasm")]
    console_error_panic_hook::set_once();

    let mut window = Window::new_hidden("Solar System from Earth");
    let app = Simulator::new(&mut window);
    window.set_framerate_limit(Some(FRAME_LIMIT));
    window.show();
    window.render_loop(app);
}
