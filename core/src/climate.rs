// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec;
use alloc::vec::Vec;

use crate::atmosphere::{Atmosphere, Ground, LAYERS, Seasons, evaporation_share};
use crate::cycle::{Circulation, Grid, Weather};
use crate::energy::SeasonalCycle;
use crate::material::Material;
use crate::scale::{SECONDS_PER_DAY, Scale};
use crate::wind::{Belts, Winds};

pub struct ClimateRun<'a> {
    pub fine: &'a Grid<'a>,
    pub downscale: u32,
    pub seasons: &'a SeasonalCycle,
    pub belts: Belts,
    pub equator_km: f64,
    pub eddy_diffusivity_m2_per_s: f64,
    pub weather_pattern: &'a [f64],
    pub weather_celsius: f64,
    pub water: &'a Material,
    pub air: &'a Material,
    pub layer_humidity: &'a [[f64; LAYERS]],
    pub soil_capacity_mm: f64,
    pub winds: Winds,
}

pub struct Climate {
    pub scale: Scale,
    pub steps: u32,
    pub rain_mm_per_year: Vec<f64>,
    pub coarse_rain_mm_per_year: Vec<f64>,
    pub evaporation_mm_per_year: Vec<f64>,
    pub column_water_mm: Vec<f64>,
    pub boundary_humidity: Vec<f64>,
    pub coarse_boundary_humidity: Vec<f64>,
    pub coarse_layer_humidity: Vec<[f64; LAYERS]>,
    pub coarse_layer_rain_mm: Vec<[f64; LAYERS]>,
    pub coarse_wetness: Vec<f64>,
    pub moved_per_step: Vec<f64>,
}

pub(crate) struct Blocks {
    fine_width: usize,
    fine_height: usize,
    factor: usize,
    pub(crate) width: usize,
    pub(crate) height: usize,
}

impl Blocks {
    pub(crate) fn new(fine_width: u32, fine_height: u32, factor: u32) -> Self {
        let factor = factor.max(1) as usize;
        let (fine_width, fine_height) = (fine_width as usize, fine_height as usize);
        Self {
            fine_width,
            fine_height,
            factor,
            width: fine_width.div_ceil(factor),
            height: fine_height.div_ceil(factor),
        }
    }

    fn block_of(&self, fine: usize) -> usize {
        let (x, y) = (fine % self.fine_width, fine / self.fine_width);
        (y / self.factor) * self.width + x / self.factor
    }

    pub(crate) fn mean(&self, values: &[f64]) -> Vec<f64> {
        let mut sums = vec![0.0; self.width * self.height];
        let mut counts = vec![0.0; self.width * self.height];
        for (i, &value) in values.iter().enumerate() {
            let block = self.block_of(i);
            sums[block] += value;
            counts[block] += 1.0;
        }
        sums.iter().zip(&counts).map(|(s, c)| s / c).collect()
    }

    pub(crate) fn majority(&self, flags: &[bool]) -> Vec<bool> {
        let shares = self.mean(
            &flags
                .iter()
                .map(|&f| if f { 1.0 } else { 0.0 })
                .collect::<Vec<_>>(),
        );
        shares.iter().map(|&s| s > 0.5).collect()
    }

    fn corners(&self, fine: usize) -> [(usize, f64); 4] {
        let factor = self.factor as f64;
        let last_row = self.height - 1;
        let (x, y) = (fine % self.fine_width, fine / self.fine_width);
        let cy = ((y as f64 + 0.5) / factor - 0.5).clamp(0.0, last_row as f64);
        let y0 = libm::floor(cy);
        let ty = cy - y0;
        let (r0, r1) = (y0 as usize, (y0 as usize + 1).min(last_row));
        let cx = (x as f64 + 0.5) / factor - 0.5;
        let x0 = libm::floor(cx);
        let tx = cx - x0;
        let width = self.width as i64;
        let c0 = (x0 as i64).rem_euclid(width) as usize;
        let c1 = (x0 as i64 + 1).rem_euclid(width) as usize;
        let at = |r: usize, c: usize| r * self.width + c;
        [
            (at(r0, c0), (1.0 - tx) * (1.0 - ty)),
            (at(r0, c1), tx * (1.0 - ty)),
            (at(r1, c0), (1.0 - tx) * ty),
            (at(r1, c1), tx * ty),
        ]
    }

    pub(crate) fn upsample(&self, coarse: &[f64]) -> Vec<f64> {
        (0..self.fine_width * self.fine_height)
            .map(|fine| self.corners(fine).iter().map(|&(c, w)| coarse[c] * w).sum())
            .collect()
    }

    pub(crate) fn upsample_by_surface(
        &self,
        coarse: &[f64],
        coarse_ocean: &[bool],
        fine_ocean: &[bool],
    ) -> Vec<f64> {
        (0..self.fine_width * self.fine_height)
            .map(|fine| {
                let corners = self.corners(fine);
                let same: f64 = corners
                    .iter()
                    .filter(|&&(c, _)| coarse_ocean[c] == fine_ocean[fine])
                    .map(|&(_, w)| w)
                    .sum();
                corners
                    .iter()
                    .filter(|&&(c, _)| same == 0.0 || coarse_ocean[c] == fine_ocean[fine])
                    .map(|&(c, w)| coarse[c] * w / if same == 0.0 { 1.0 } else { same })
                    .sum()
            })
            .collect()
    }
}

impl ClimateRun<'_> {
    pub fn run(&self) -> Climate {
        let fine = self.fine;
        let blocks = Blocks::new(fine.width, fine.height, self.downscale);
        let elevation = blocks.mean(fine.elevation);
        let celsius = blocks.mean(fine.celsius);
        let condensation = blocks.mean(fine.condensation_celsius);
        let ocean = blocks.majority(fine.ocean);
        let pattern = if self.weather_pattern.is_empty() {
            Vec::new()
        } else {
            blocks.mean(self.weather_pattern)
        };
        let coarse = Grid {
            width: blocks.width as u32,
            height: blocks.height as u32,
            elevation: &elevation,
            celsius: &celsius,
            condensation_celsius: &condensation,
            ocean: &ocean,
            sun_latitude: 0.0,
        };

        let scale = Scale::new(
            self.equator_km,
            coarse.width,
            self.eddy_diffusivity_m2_per_s,
            self.winds.speed_metres_per_second,
        );
        let steps = scale.steps(self.seasons.year_seconds / SECONDS_PER_DAY);
        let weather = Weather {
            pattern: &pattern,
            amplitude_celsius: self.weather_celsius,
        };
        let outcome = Atmosphere {
            grid: &coarse,
            water: self.water,
            circulation: Circulation {
                gas_spread: scale.mixing_share,
                wind_carry: scale.wind_carry,
                age_steps: steps,
                belts: self.belts,
                winds: self.winds,
            },
            weather: &weather,
            seasons: Some(Seasons {
                cycle: self.seasons,
                step_seconds: scale.step_seconds,
                air: self.air,
                layer_humidity: self.layer_humidity,
            }),
            ground: Ground {
                soil_capacity_mm: self.soil_capacity_mm,
                evaporation_share: evaporation_share(&scale),
            },
        }
        .run(steps);

        let rain: Vec<f64> = outcome
            .per_step
            .iter()
            .map(|e| scale.per_year(e.precipitation_mm))
            .collect();
        let evaporation: Vec<f64> = outcome
            .per_step
            .iter()
            .map(|e| scale.per_year(e.evaporation_mm))
            .collect();
        let humidity: Vec<f64> = outcome
            .per_step
            .iter()
            .map(|e| e.boundary_humidity)
            .collect();
        let column: Vec<f64> = (0..coarse.len())
            .map(|i| outcome.air.water_mm.iter().map(|layer| layer[i]).sum())
            .collect();
        Climate {
            scale,
            steps,
            rain_mm_per_year: blocks.upsample(&rain),
            coarse_rain_mm_per_year: rain,
            evaporation_mm_per_year: blocks.upsample_by_surface(&evaporation, &ocean, fine.ocean),
            column_water_mm: blocks.upsample(&column),
            boundary_humidity: blocks.upsample(&humidity),
            coarse_boundary_humidity: humidity,
            coarse_layer_humidity: outcome.per_step.iter().map(|e| e.layer_humidity).collect(),
            coarse_layer_rain_mm: outcome.per_step.iter().map(|e| e.layer_rain_mm).collect(),
            coarse_wetness: outcome.per_step.iter().map(|e| e.wetness).collect(),
            moved_per_step: outcome.moved_per_step,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coarsening_averages_each_block() {
        let blocks = Blocks::new(4, 2, 2);
        assert_eq!(
            blocks.mean(&[1.0, 3.0, 10.0, 10.0, 5.0, 7.0, 0.0, 20.0]),
            [4.0, 10.0]
        );
    }

    #[test]
    fn a_block_is_sea_when_most_of_it_is_sea() {
        let blocks = Blocks::new(4, 2, 2);
        let ocean = [true, true, true, false, true, false, false, false];
        assert_eq!(blocks.majority(&ocean), [true, false]);
    }

    #[test]
    fn upsampling_keeps_a_uniform_field_uniform() {
        let blocks = Blocks::new(8, 4, 2);
        assert!(blocks.upsample(&[3.0; 8]).iter().all(|&v| v == 3.0));
    }

    #[test]
    fn land_takes_its_evaporation_only_from_land_blocks() {
        let blocks = Blocks::new(4, 2, 2);
        let coarse_ocean = [false, true];
        let fine_ocean = [false, false, true, true, false, false, true, true];
        let fine = blocks.upsample_by_surface(&[400.0, 1200.0], &coarse_ocean, &fine_ocean);
        assert_eq!(
            fine,
            [400.0, 400.0, 1200.0, 1200.0, 400.0, 400.0, 1200.0, 1200.0]
        );
    }

    #[test]
    fn upsampling_blends_between_neighbouring_blocks_and_wraps() {
        let blocks = Blocks::new(4, 2, 2);
        let fine = blocks.upsample(&[0.0, 8.0]);
        assert_eq!(fine[..4], [2.0, 2.0, 6.0, 6.0]);
    }
}
