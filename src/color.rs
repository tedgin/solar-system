fn verify_channel_normalized(channel: f64) {
    if channel.is_nan() || channel < 0. || channel > 1. {
        panic!("Color channel value {channel} is not in the range [0, 1]");
    }
}

pub struct Color {
    red: f64,
    green: f64,
    blue: f64,
}

impl Color {
    pub fn blue(&self) -> f64 {
        self.blue
    }

    pub fn green(&self) -> f64 {
        self.green
    }

    pub fn red(&self) -> f64 {
        self.red
    }
}

pub fn mk_color(red: f64, green: f64, blue: f64) -> Color {
    verify_channel_normalized(red);
    verify_channel_normalized(green);
    verify_channel_normalized(blue);
    Color { red, green, blue }
}

pub fn mk_gray(shade: f64) -> Color {
    verify_channel_normalized(shade);
    Color {
        red: shade,
        green: shade,
        blue: shade,
    }
}

pub const BLUE: Color = Color {
    red: 0.,
    green: 0.,
    blue: 1.,
};

pub const RED: Color = Color {
    red: 1.,
    green: 0.,
    blue: 0.,
};

pub const YELLOW: Color = Color {
    red: 1.,
    green: 1.,
    blue: 0.,
};

pub const WHITE: Color = Color {
    red: 1.,
    green: 1.,
    blue: 1.,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "Color channel value NaN is not in the range [0, 1]")]
    fn test_verify_channel_normalized_nan() {
        verify_channel_normalized(f64::NAN)
    }

    #[test]
    #[should_panic(expected = "")]
    fn test_verify_channel_normalized_neg() {
        verify_channel_normalized(-f64::from_bits(0x1))
    }

    #[test]
    fn test_verify_channel_normalized_good() {
        verify_channel_normalized(0.);
        verify_channel_normalized(1.);
    }

    #[test]
    #[should_panic(expected = "")]
    fn test_verify_channel_normalized_too_large() {
        verify_channel_normalized(1. + f64::EPSILON)
    }

    #[test]
    fn test_color_new() {
        let c = mk_color(0.1, 0.2, 0.3);
        assert_eq!(c.red, 0.1);
        assert_eq!(c.green, 0.2);
        assert_eq!(c.blue, 0.3);
    }

    #[test]
    #[should_panic(expected = "")]
    fn test_color_new_bad_red() {
        mk_color(2., 0., 0.);
    }

    #[test]
    #[should_panic(expected = "")]
    fn test_color_new_bad_green() {
        mk_color(0., 2., 0.);
    }

    #[test]
    #[should_panic(expected = "")]
    fn test_color_new_bad_blue() {
        mk_color(0., 0., 2.);
    }

    #[test]
    fn test_color_gray() {
        let gray = mk_gray(0.5);
        assert_eq!(gray.red, 0.5);
        assert_eq!(gray.green, 0.5);
        assert_eq!(gray.blue, 0.5);
    }
}
