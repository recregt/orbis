// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use core::f64::consts::{FRAC_PI_2, PI};

use crate::material::Material;
use crate::scale::{
    GRAVITY_METRES_PER_SECOND_SQUARED, LAPSE_RATE_CELSIUS_PER_METRE,
    SEA_LEVEL_AIR_KG_PER_CUBIC_METRE,
};

pub const SURFACE_DRAG_COEFFICIENT: f64 = 1.3e-3;

pub fn eastward_stress(latitude: f64, belts: &Belts, wind_metres_per_second: f64) -> f64 {
    let flow = layer_flow(0, latitude, belts);
    let speed = wind_metres_per_second * flow.south.abs();
    let balance = if flow.east > 0.0 {
        belts.westerly_balance()
    } else {
        1.0
    };
    flow.east.signum()
        * balance
        * SEA_LEVEL_AIR_KG_PER_CUBIC_METRE
        * SURFACE_DRAG_COEFFICIENT
        * speed
        * speed
}

fn belt_torque(from_radians: f64, to_radians: f64) -> f64 {
    let width = to_radians - from_radians;
    if width <= 0.0 {
        return 0.0;
    }
    let spread = libm::sin(2.0 * to_radians) - libm::sin(2.0 * from_radians);
    let wave = 2.0 * PI / width;
    (width + spread / 2.0 + 2.0 * spread / (wave * wave - 4.0)) / 4.0
}

pub fn ground_friction_w_per_m2(north_south_work_w_per_m2: f64, air: &Material) -> f64 {
    north_south_work_w_per_m2 * (1.0 - air.free_air_dissipation_share)
}

pub fn surface_wind_metres_per_second(ground_friction_w_per_m2: f64) -> f64 {
    libm::cbrt(
        ground_friction_w_per_m2.max(0.0)
            / (SEA_LEVEL_AIR_KG_PER_CUBIC_METRE * SURFACE_DRAG_COEFFICIENT),
    )
}

pub fn equatorial_width_metres(planet: &Planet, air: &Material, surface_kelvin: f64) -> f64 {
    let dry_lapse = GRAVITY_METRES_PER_SECOND_SQUARED / air.specific_heat;
    let buoyancy_squared = GRAVITY_METRES_PER_SECOND_SQUARED / surface_kelvin
        * (dry_lapse - LAPSE_RATE_CELSIUS_PER_METRE);
    let wave_speed = libm::sqrt(buoyancy_squared.max(0.0)) * planet.tropopause_metres / PI;
    let beta = 2.0 * planet.rotation_radians_per_second / planet.radius_metres;
    libm::sqrt(wave_speed / beta)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Winds {
    pub planet: Planet,
    pub equatorial_width_metres: f64,
    pub speed_metres_per_second: f64,
}

#[cfg(test)]
pub(crate) const EARTH_WINDS: Winds = Winds {
    planet: Planet {
        radius_metres: 6.378e6,
        rotation_radians_per_second: 7.292e-5,
        tropopause_metres: 12_000.0,
        radiative_contrast: 0.33,
    },
    equatorial_width_metres: 1.3e6,
    speed_metres_per_second: 7.4,
};

pub fn boundary_south_fraction(
    grid_latitude: f64,
    belt_latitude: f64,
    belts: &Belts,
    winds: &Winds,
) -> f64 {
    let shifted = grid_latitude - belt_latitude;
    let toward = layer_flow(0, shifted, belts).south;
    if toward == 0.0 || winds.speed_metres_per_second <= 0.0 {
        return 0.0;
    }
    let speed = ekman_speed(
        eastward_stress(shifted, belts, winds.speed_metres_per_second),
        -grid_latitude * FRAC_PI_2,
        &winds.planet,
        winds.equatorial_width_metres,
        crate::atmosphere::boundary_air_kg_per_square_metre(),
    );
    (speed / winds.speed_metres_per_second).min(1.0) * toward.signum()
}

pub fn ekman_speed(
    stress: f64,
    north_radians: f64,
    planet: &Planet,
    equatorial_width_metres: f64,
    boundary_air_kg_per_square_metre: f64,
) -> f64 {
    let omega = planet.rotation_radians_per_second;
    let coriolis = 2.0 * omega * libm::sin(north_radians);
    let equatorial = 2.0 * omega / planet.radius_metres * equatorial_width_metres;
    stress.abs() * coriolis.abs()
        / ((coriolis * coriolis + equatorial * equatorial) * boundary_air_kg_per_square_metre)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Direction {
    pub east: i64,
    pub south: i64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Flow {
    pub east: f64,
    pub south: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Belts {
    pub hadley_edge: f64,
    pub polar_front: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Planet {
    pub radius_metres: f64,
    pub rotation_radians_per_second: f64,
    pub tropopause_metres: f64,
    pub radiative_contrast: f64,
}

impl Belts {
    pub const EARTHLIKE: Belts = Belts {
        hadley_edge: 1.0 / 3.0,
        polar_front: 2.0 / 3.0,
    };

    pub fn of(planet: &Planet) -> Self {
        let spin = planet.rotation_radians_per_second * planet.radius_metres;
        let edge_radians = libm::sqrt(
            5.0 * GRAVITY_METRES_PER_SECOND_SQUARED
                * planet.tropopause_metres
                * planet.radiative_contrast.max(0.0)
                / (3.0 * spin * spin),
        );
        let hadley_edge = (edge_radians / FRAC_PI_2).min(1.0);
        Self {
            hadley_edge,
            polar_front: (hadley_edge + 1.0) / 2.0,
        }
    }

    fn westerly_balance(&self) -> f64 {
        let (edge, front) = (self.hadley_edge * FRAC_PI_2, self.polar_front * FRAC_PI_2);
        let westerlies = belt_torque(edge, front);
        if westerlies <= 0.0 {
            return 1.0;
        }
        (belt_torque(0.0, edge) + belt_torque(front, FRAC_PI_2)) / westerlies
    }

    fn cell(&self, latitude: f64) -> (Cell, f64) {
        let away = libm::fabs(latitude).min(1.0);
        let within = |from: f64, to: f64| {
            if to > from {
                (away - from) / (to - from)
            } else {
                0.0
            }
        };
        if away < self.hadley_edge {
            (Cell::Hadley, within(0.0, self.hadley_edge))
        } else if away < self.polar_front {
            (Cell::Ferrel, within(self.hadley_edge, self.polar_front))
        } else {
            (Cell::Polar, within(self.polar_front, 1.0))
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Cell {
    Hadley,
    Ferrel,
    Polar,
}

fn toward_equator(latitude: f64) -> f64 {
    if latitude < 0.0 { 1.0 } else { -1.0 }
}

pub fn surface_meridional(latitude: f64, belts: &Belts) -> f64 {
    let (cell, within) = belts.cell(latitude);
    let strength = libm::sin(PI * within);
    match cell {
        Cell::Hadley | Cell::Polar => toward_equator(latitude) * strength,
        Cell::Ferrel => -toward_equator(latitude) * strength,
    }
}

pub fn layer_flow(layer: usize, latitude: f64, belts: &Belts) -> Flow {
    let tropics = belts.cell(latitude).0 == Cell::Hadley;
    match layer {
        0 => Flow {
            east: prevailing(latitude, belts).east as f64,
            south: surface_meridional(latitude, belts),
        },
        1 => Flow {
            east: if tropics { -1.0 } else { 1.0 },
            south: 0.0,
        },
        _ => Flow {
            east: 1.0,
            south: -surface_meridional(latitude, belts),
        },
    }
}

pub fn prevailing(latitude: f64, belts: &Belts) -> Direction {
    let toward = toward_equator(latitude) as i64;
    match belts.cell(latitude).0 {
        Cell::Hadley | Cell::Polar => Direction {
            east: -1,
            south: toward,
        },
        Cell::Ferrel => Direction {
            east: 1,
            south: -toward,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EARTH: Belts = Belts::EARTHLIKE;

    #[test]
    fn surface_stress_puts_back_as_much_spin_as_it_takes_out() {
        for belts in [EARTH, Belts::of(&EARTH_WINDS.planet)] {
            let strips = 20_000;
            let (mut torque, mut easterly) = (0.0, 0.0);
            for k in 0..strips {
                let north = (k as f64 + 0.5) / strips as f64 * FRAC_PI_2;
                let cos = libm::cos(north);
                let push = eastward_stress(-north / FRAC_PI_2, &belts, 7.0) * cos * cos;
                torque += push;
                easterly += push.min(0.0);
            }
            assert!(
                torque.abs() < 1e-6 * easterly.abs(),
                "{torque} of {easterly}"
            );
        }
        assert!((1.5..2.5).contains(&EARTH.westerly_balance()));
    }

    #[test]
    fn the_ground_spends_what_the_free_air_leaves_of_the_north_south_work() {
        let air = crate::material::air();
        let work = 5.0;
        let ground = ground_friction_w_per_m2(work, &air);
        assert_eq!(ground, work * (1.0 - air.free_air_dissipation_share));
        assert!(ground < work);
    }

    #[test]
    fn surface_wind_spends_the_ground_friction_through_drag() {
        let friction = 0.6;
        let speed = surface_wind_metres_per_second(friction);
        let spent = SEA_LEVEL_AIR_KG_PER_CUBIC_METRE * SURFACE_DRAG_COEFFICIENT * speed.powi(3);
        assert!((spent - friction).abs() < friction * 1e-12, "{spent}");
        assert_eq!(surface_wind_metres_per_second(-0.1), 0.0);
    }

    #[test]
    fn surface_air_flows_toward_the_equator_and_sixty_degrees() {
        let north = |degrees: f64| -degrees / 90.0;
        assert!(surface_meridional(north(15.0), &EARTH) > 0.0);
        assert!(surface_meridional(north(45.0), &EARTH) < 0.0);
        assert!(surface_meridional(north(75.0), &EARTH) > 0.0);
        assert!(surface_meridional(-north(15.0), &EARTH) < 0.0);
        assert!(surface_meridional(0.0, &EARTH).abs() < 1e-12);
        assert!(surface_meridional(north(30.0), &EARTH).abs() < 1e-12);
    }

    #[test]
    fn the_upper_layer_returns_what_the_surface_brings() {
        for latitude in [-0.8, -0.4, -0.1, 0.2, 0.5, 0.9] {
            assert_eq!(
                layer_flow(2, latitude, &EARTH).south,
                -layer_flow(0, latitude, &EARTH).south
            );
        }
    }

    #[test]
    fn trade_winds_blow_west_toward_the_equator() {
        assert_eq!(prevailing(-0.1, &EARTH), Direction { east: -1, south: 1 });
        assert_eq!(
            prevailing(0.1, &EARTH),
            Direction {
                east: -1,
                south: -1
            }
        );
    }

    #[test]
    fn westerlies_blow_east_toward_the_poles() {
        assert_eq!(prevailing(-0.5, &EARTH), Direction { east: 1, south: -1 });
        assert_eq!(prevailing(0.5, &EARTH), Direction { east: 1, south: 1 });
    }

    #[test]
    fn polar_easterlies_blow_west_toward_the_equator() {
        assert_eq!(prevailing(-0.9, &EARTH), Direction { east: -1, south: 1 });
        assert_eq!(
            prevailing(0.9, &EARTH),
            Direction {
                east: -1,
                south: -1
            }
        );
    }

    fn earth(rotation_scale: f64) -> Planet {
        Planet {
            radius_metres: 6.378e6,
            rotation_radians_per_second: 7.292e-5 * rotation_scale,
            tropopause_metres: 12_000.0,
            radiative_contrast: 0.34,
        }
    }

    #[test]
    fn earths_tropical_cell_reaches_the_subtropics() {
        let degrees = Belts::of(&earth(1.0)).hadley_edge * 90.0;
        assert!((22.0..34.0).contains(&degrees), "{degrees}");
    }

    #[test]
    fn slower_spin_widens_the_tropical_cell_until_it_reaches_the_pole() {
        let fast = Belts::of(&earth(2.0));
        let earthlike = Belts::of(&earth(1.0));
        let slow = Belts::of(&earth(0.1));
        assert!(fast.hadley_edge < earthlike.hadley_edge);
        assert!(earthlike.hadley_edge < slow.hadley_edge);
        assert_eq!(slow.hadley_edge, 1.0);
        assert_eq!(slow.polar_front, 1.0);
    }

    #[test]
    fn a_stronger_equator_to_pole_contrast_widens_the_tropical_cell() {
        let mild = Belts::of(&Planet {
            radiative_contrast: 0.1,
            ..earth(1.0)
        });
        let harsh = Belts::of(&Planet {
            radiative_contrast: 0.3,
            ..earth(1.0)
        });
        assert!(mild.hadley_edge < harsh.hadley_edge);
    }
}
