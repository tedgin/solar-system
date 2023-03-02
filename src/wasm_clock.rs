use std::time::Duration;

use solar_system::clock::TimeKeeper;

#[derive(Clone)]
pub struct Clock {}

impl TimeKeeper for Clock {
    fn init() -> Self {
        Self {}
    }

    fn advance(&mut self, _: Duration) {}
}
