use std::thread;
use std::time::{Duration, Instant};

use solar_system::clock::TimeKeeper;

#[derive(Clone)]
pub struct Clock {
    step_start_time: Instant,
}

impl TimeKeeper for Clock {
    fn init() -> Self {
        Self {
            step_start_time: Instant::now(),
        }
    }

    fn advance(&mut self, step: Duration) {
        if let Some(remaining_time) = step.checked_sub(self.step_start_time.elapsed()) {
            thread::sleep(remaining_time);
        }
        self.step_start_time = Instant::now();
    }
}
