pub const LAPSE_RATE_CELSIUS_PER_METRE: f64 = 0.0065;
pub const VAPOR_SCALE_HEIGHT_METRES: f64 = 2000.0;
pub const SECONDS_PER_DAY: f64 = 86_400.0;
pub const SECONDS_PER_HOUR: f64 = 3_600.0;
pub const DAYS_PER_YEAR: f64 = 365.25;
pub const FREEZING_DEPRESSION_CELSIUS_PER_GRAM_PER_KG: f64 = 0.054;
pub const GRAVITY_METRES_PER_SECOND_SQUARED: f64 = 9.80665;
pub const SEA_LEVEL_AIR_KG_PER_CUBIC_METRE: f64 = 1.225;
pub const SURFACE_MOISTURE_TRANSFER_COEFFICIENT: f64 = 1.2e-3;
pub const STEFAN_BOLTZMANN_W_PER_M2_PER_K4: f64 = 5.670_374_419e-8;
pub const KELVIN_AT_ZERO_CELSIUS: f64 = 273.15;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Scale {
    pub cell_metres: f64,
    pub step_seconds: f64,
    pub wind_metres_per_second: f64,
    pub wind_carry: f64,
    pub mixing_share: f64,
}

impl Scale {
    pub fn new(
        equator_km: f64,
        width: u32,
        eddy_diffusivity_m2_per_s: f64,
        wind_metres_per_second: f64,
    ) -> Self {
        let cell_metres = equator_km * 1000.0 / f64::from(width);
        let wind_rate = wind_metres_per_second / cell_metres;
        let step_seconds = 1.0 / (wind_rate + mixing_rate(eddy_diffusivity_m2_per_s, cell_metres));
        let wind_carry = wind_rate * step_seconds;
        Self {
            cell_metres,
            step_seconds,
            wind_metres_per_second,
            wind_carry,
            mixing_share: 1.0 - wind_carry,
        }
    }

    pub fn millimetres(amount: f64) -> f64 {
        amount * VAPOR_SCALE_HEIGHT_METRES
    }

    pub fn amount(millimetres: f64) -> f64 {
        millimetres / VAPOR_SCALE_HEIGHT_METRES
    }

    pub fn steps_per_day(&self) -> f64 {
        SECONDS_PER_DAY / self.step_seconds
    }

    pub fn per_year(&self, per_step: f64) -> f64 {
        per_step * self.steps_per_day() * DAYS_PER_YEAR
    }

    pub fn per_step(&self, per_day: f64) -> f64 {
        per_day / self.steps_per_day()
    }

    pub fn days(&self, steps: u32) -> f64 {
        f64::from(steps) / self.steps_per_day()
    }

    pub fn steps(&self, days: f64) -> u32 {
        libm::ceil(days * self.steps_per_day()) as u32
    }
}

pub const EIGHT_NEIGHBOR_SPREAD_PER_DIFFUSIVITY: f64 = 8.0 / 3.0;

pub fn mixing_rate(eddy_diffusivity_m2_per_s: f64, cell_metres: f64) -> f64 {
    EIGHT_NEIGHBOR_SPREAD_PER_DIFFUSIVITY * eddy_diffusivity_m2_per_s / (cell_metres * cell_metres)
}

pub fn salinity_grams_per_kg(salt_amount: f64, depth_metres: f64) -> f64 {
    Scale::millimetres(salt_amount) / depth_metres
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_step_is_the_longest_in_which_wind_and_mixing_empty_a_cell_at_most_once() {
        let wind = 10.0;
        let scale = Scale::new(40_000.0, 400, 1.0e6, wind);
        assert_eq!(scale.cell_metres, 100_000.0);
        let close = |a: f64, b: f64| (a - b).abs() <= b.abs() * 1e-12;
        let mixing = 8.0 / 3.0 * 1.0e-4;
        assert!(close(scale.step_seconds, 1.0 / (1.0e-4 + mixing)));
        assert!(close(
            scale.wind_carry,
            wind * scale.step_seconds / scale.cell_metres
        ));
        assert!(close(scale.mixing_share, mixing * scale.step_seconds));
        assert!(scale.wind_carry + scale.mixing_share <= 1.0);
    }

    #[test]
    fn equal_shares_to_eight_neighbors_diffuse_at_the_stated_diffusivity() {
        let (diffusivity, cell, seconds) = (1.0e6, 100_000.0, 1000.0);
        let per_neighbor = mixing_rate(diffusivity, cell) * seconds / 8.0;
        let field = |dx: f64, _dy: f64| (dx * cell) * (dx * cell);
        let gained: f64 = [-1.0, 0.0, 1.0]
            .iter()
            .flat_map(|&dy| [-1.0, 0.0, 1.0].map(|dx| (dx, dy)))
            .filter(|&(dx, dy)| (dx, dy) != (0.0, 0.0))
            .map(|(dx, dy)| per_neighbor * (field(dx, dy) - field(0.0, 0.0)))
            .sum();
        let laplacian_of_x_squared = 2.0;
        let expected = diffusivity * seconds * laplacian_of_x_squared;
        assert!(
            (gained - expected).abs() <= expected * 1e-12,
            "{gained} {expected}"
        );
    }

    #[test]
    fn a_finer_grid_takes_shorter_steps() {
        let coarse = Scale::new(40_000.0, 128, 1.0e6, 10.0);
        let fine = Scale::new(40_000.0, 512, 1.0e6, 10.0);
        assert!(fine.step_seconds < coarse.step_seconds / 4.0);
        assert!(fine.wind_carry < coarse.wind_carry);
    }

    #[test]
    fn amounts_convert_to_water_depth_through_the_vapor_column() {
        assert_eq!(Scale::millimetres(0.5), 1000.0);
        assert_eq!(Scale::amount(1000.0), 0.5);
    }

    #[test]
    fn salinity_is_salt_mass_over_water_mass() {
        let salt_for_35_grams_in_1000_metres = Scale::amount(35.0 * 1000.0);
        assert_eq!(
            salinity_grams_per_kg(salt_for_35_grams_in_1000_metres, 1000.0),
            35.0
        );
    }
}
