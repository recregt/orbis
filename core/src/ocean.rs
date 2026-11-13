// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec;
use alloc::vec::Vec;
use core::f64::consts::{PI, TAU};

use crate::cycle::Grid;
use crate::energy::OCEAN_MIXED_LAYER_METRES;
use crate::material::Material;
use crate::scale::GRAVITY_METRES_PER_SECOND_SQUARED;
use crate::wind::{self, Belts};

pub const SVERDRUP_CUBIC_METRES_PER_SECOND: f64 = 1.0e6;

pub fn equatorial_width_metres(
    water: &Material,
    equator_to_pole_kelvin: f64,
    radius_metres: f64,
    rotation_radians_per_second: f64,
) -> f64 {
    let reduced_gravity =
        GRAVITY_METRES_PER_SECOND_SQUARED * water.expands_per_kelvin * equator_to_pole_kelvin;
    let wave_speed = libm::sqrt(reduced_gravity.max(0.0) * OCEAN_MIXED_LAYER_METRES);
    let beta = 2.0 * rotation_radians_per_second / radius_metres;
    libm::sqrt(wave_speed / beta)
}

pub struct Ocean<'a> {
    pub grid: &'a Grid<'a>,
    pub belts: Belts,
    pub radius_metres: f64,
    pub rotation_radians_per_second: f64,
    pub water_kg_per_cubic_metre: f64,
    pub equatorial_width_metres: f64,
    pub wind_metres_per_second: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Currents {
    pub gyre: Vec<f64>,
    pub overturning: Vec<f64>,
}

impl Currents {
    pub fn still(cells: usize) -> Self {
        Self {
            gyre: vec![0.0; cells],
            overturning: vec![0.0; cells],
        }
    }

    pub fn strongest_gyre_sverdrups(&self) -> f64 {
        self.gyre.iter().fold(0.0, |m: f64, &v| m.max(v.abs())) / SVERDRUP_CUBIC_METRES_PER_SECOND
    }

    pub fn overturning_sverdrups(&self, width: usize) -> Vec<f64> {
        self.overturning
            .chunks(width)
            .map(|row| row.iter().sum::<f64>() / SVERDRUP_CUBIC_METRES_PER_SECOND)
            .collect()
    }
}

impl Ocean<'_> {
    fn latitude_north(&self, row: usize) -> f64 {
        (0.5 - (row as f64 + 0.5) / f64::from(self.grid.height)) * PI
    }

    pub fn eastward_stress(&self, row: usize) -> f64 {
        wind::eastward_stress(
            self.grid.latitude(row as u32),
            &self.belts,
            self.wind_metres_per_second,
        )
    }

    fn stress_curl(&self, row: usize) -> f64 {
        let height = self.grid.height as usize;
        let d_phi = PI / height as f64;
        let weighted = |r: usize| self.eastward_stress(r) * libm::cos(self.latitude_north(r));
        let north = if row == 0 {
            weighted(row)
        } else {
            weighted(row - 1)
        };
        let south = if row + 1 == height {
            weighted(row)
        } else {
            weighted(row + 1)
        };
        let span = d_phi
            * (if row == 0 || row + 1 == height {
                1.0
            } else {
                2.0
            });
        -(north - south) / span / (self.radius_metres * libm::cos(self.latitude_north(row)))
    }

    pub fn currents(&self) -> Currents {
        Currents {
            gyre: self.gyres(),
            overturning: self.overturning(),
        }
    }

    fn overturning(&self) -> Vec<f64> {
        let grid = self.grid;
        let (width, height) = (grid.width as usize, grid.height as usize);
        let omega = self.rotation_radians_per_second;
        let equatorial_coriolis = 2.0 * omega / self.radius_metres * self.equatorial_width_metres;
        let mut northward = vec![0.0; grid.len()];
        for row in 0..height.saturating_sub(1) {
            let latitude = (0.5 - (row as f64 + 1.0) / height as f64) * PI;
            let coriolis = 2.0 * omega * libm::sin(latitude);
            let stress = (self.eastward_stress(row) + self.eastward_stress(row + 1)) / 2.0;
            let per_metre = -stress * coriolis
                / (self.water_kg_per_cubic_metre
                    * (coriolis * coriolis + equatorial_coriolis * equatorial_coriolis));
            let face_metres = self.radius_metres * libm::cos(latitude) * TAU / width as f64;
            for x in 0..width {
                let (here, south) = (row * width + x, (row + 1) * width + x);
                if grid.ocean[here] && grid.ocean[south] {
                    northward[south] = per_metre * face_metres;
                }
            }
        }
        northward
    }

    fn gyres(&self) -> Vec<f64> {
        let grid = self.grid;
        let (width, height) = (grid.width as usize, grid.height as usize);
        let a = self.radius_metres;
        let omega = self.rotation_radians_per_second;
        let d_lambda = TAU / width as f64;
        let mut transport = vec![0.0; grid.len()];
        for row in 0..height {
            let c = libm::cos(self.latitude_north(row));
            let cell_metres = a * c * d_lambda;
            let beta = 2.0 * omega * c / a;
            let sverdrup = -self.stress_curl(row) / (self.water_kg_per_cubic_metre * beta);
            let boundary_cells = 1.0 / (c * c);
            let ocean = |x: usize| grid.ocean[row * width + x % width];
            let Some(shore) = (0..width).find(|&x| !ocean(x)) else {
                continue;
            };
            let mut x = shore;
            while x < shore + width {
                if !ocean(x) {
                    x += 1;
                    continue;
                }
                let west = x;
                while x < shore + width && ocean(x) {
                    x += 1;
                }
                let east = x;
                for cell in west..east {
                    let to_east_coast = (east - cell) as f64 - 0.5;
                    let from_west_coast = (cell - west) as f64 + 0.5;
                    let interior = sverdrup * to_east_coast * cell_metres;
                    let closure = 1.0 - libm::exp(-from_west_coast / boundary_cells);
                    transport[row * width + cell % width] = interior * closure;
                }
            }
        }
        transport
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Basin {
        width: u32,
        height: u32,
        zeros: Vec<f64>,
        ocean: Vec<bool>,
    }

    impl Basin {
        fn new(width: u32, height: u32, land_columns: u32) -> Self {
            let cells = (width * height) as usize;
            Self {
                width,
                height,
                zeros: vec![0.0; cells],
                ocean: (0..cells as u32)
                    .map(|i| i % width >= land_columns)
                    .collect(),
            }
        }

        fn currents(&self, rotation_scale: f64) -> (Currents, Vec<f64>) {
            let grid = Grid {
                width: self.width,
                height: self.height,
                elevation: &self.zeros,
                celsius: &self.zeros,
                condensation_celsius: &self.zeros,
                ocean: &self.ocean,
                sun_latitude: 0.0,
            };
            let ocean = Ocean {
                grid: &grid,
                belts: Belts::EARTHLIKE,
                radius_metres: 6.378e6,
                rotation_radians_per_second: 7.292e-5 * rotation_scale,
                water_kg_per_cubic_metre: 1025.0,
                equatorial_width_metres: EARTH_EQUATORIAL_WIDTH_METRES,
                wind_metres_per_second: crate::wind::EARTH_WINDS.speed_metres_per_second,
            };
            let stress = (0..self.height as usize)
                .map(|r| ocean.eastward_stress(r))
                .collect();
            (ocean.currents(), stress)
        }
    }

    const EARTH_EQUATORIAL_WIDTH_METRES: f64 = 3.0e5;

    fn at(currents: &Currents, width: u32, x: u32, row: u32) -> f64 {
        currents.gyre[(row * width + x) as usize]
    }

    fn zonal_overturning(basin: &Basin, rotation_scale: f64) -> Vec<f64> {
        let (currents, _) = basin.currents(rotation_scale);
        currents.overturning_sverdrups(basin.width as usize)
    }

    #[test]
    fn the_trades_push_surface_water_poleward_and_the_westerlies_push_it_back() {
        let (width, height) = (72, 36);
        let basin = Basin::new(width, height, 12);
        let northward = zonal_overturning(&basin, 1.0);
        let face_at = |north_degrees: f64| ((90.0 - north_degrees) / 5.0) as usize;
        assert!(
            northward[face_at(15.0)] > 0.0,
            "{}",
            northward[face_at(15.0)]
        );
        assert!(
            northward[face_at(-15.0)] < 0.0,
            "{}",
            northward[face_at(-15.0)]
        );
        assert!(
            northward[face_at(45.0)] < 0.0,
            "{}",
            northward[face_at(45.0)]
        );
        assert!(
            northward[face_at(-45.0)] > 0.0,
            "{}",
            northward[face_at(-45.0)]
        );
    }

    #[test]
    fn no_water_crosses_the_equator_where_the_coriolis_force_vanishes() {
        let height = 64;
        let basin = Basin::new(2 * height, height, height / 3);
        let northward = zonal_overturning(&basin, 1.0);
        assert!(northward.iter().all(|v| v.is_finite()), "{northward:?}");
        assert_eq!(northward[height as usize / 2], 0.0);
    }

    #[test]
    fn trades_and_westerlies_have_earthlike_strength() {
        let basin = Basin::new(72, 36, 12);
        let (_, stress) = basin.currents(1.0);
        let row_at = |north_degrees: f64| ((90.0 - north_degrees) / 5.0) as usize;
        assert!(stress[row_at(15.0)] < -0.05 && stress[row_at(15.0)] > -0.25);
        assert!(stress[row_at(45.0)] > 0.05 && stress[row_at(45.0)] < 0.25);
    }

    #[test]
    fn the_subtropical_gyre_turns_clockwise_in_the_north_with_earthlike_strength() {
        let (width, height) = (72, 36);
        let basin = Basin::new(width, height, 12);
        let (currents, _) = basin.currents(1.0);
        let subtropics = 12;
        let middle = at(&currents, width, 42, subtropics) / SVERDRUP_CUBIC_METRES_PER_SECOND;
        assert!(middle > 0.0, "clockwise in the north: {middle}");
        assert!(
            (10.0..200.0).contains(&currents.strongest_gyre_sverdrups()),
            "{} Sv",
            currents.strongest_gyre_sverdrups()
        );
        let southern = at(&currents, width, 42, height - 1 - subtropics);
        assert!(southern < 0.0, "anticlockwise in the south: {southern}");
    }

    #[test]
    fn the_return_flow_crowds_against_the_western_shore() {
        let (width, height) = (72, 36);
        let basin = Basin::new(width, height, 12);
        let (currents, _) = basin.currents(1.0);
        let row = 12;
        let strongest = (12..width)
            .max_by(|&a, &b| at(&currents, width, a, row).total_cmp(&at(&currents, width, b, row)))
            .unwrap_or(12);
        assert!(
            (13..17).contains(&strongest),
            "the gyre peaks a few cells off the western shore, at {strongest}"
        );
        let at_the_wall = at(&currents, width, 12, row);
        let peak = at(&currents, width, strongest, row);
        assert!(at_the_wall < 0.6 * peak, "wall {at_the_wall} peak {peak}");
        let off_eastern_shore = at(&currents, width, width - 1, row);
        assert!(
            off_eastern_shore < 0.1 * peak,
            "east {off_eastern_shore} peak {peak}"
        );
    }
}
