//#[cfg(feature = "wasm")]
//extern crate console_error_panic_hook;

use solar_system;

pub fn main() {
    //    #[cfg(feature = "wasm")]
    //    console_error_panic_hook::set_once();

    solar_system::run();
}
