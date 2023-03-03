use std::time::Duration;

pub trait TimeKeeper {
    fn init() -> Self;
    fn advance(&mut self, step: Duration);
}
