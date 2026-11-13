// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::collections::BinaryHeap;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;
use serde::{Deserialize, Serialize};

use crate::Seed;
use crate::hash::{Purpose, lattice, unit};

const RHO_LAYER: u64 = u64::MAX;

pub const AREA_EXPONENT: f64 = 0.5;
const FLOOD_RISE_PER_METRE: f64 = 1e-6;
pub const FLATTEST_SLOPE: f64 = FLOOD_RISE_PER_METRE;
const CREEP_SWEEPS: usize = 4;
const SETTLED_SHARE_OF_UPLIFT: f64 = 0.01;
const STEPS_WITHOUT_PROGRESS: usize = 8;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Laws {
    pub erodibility: f64,
    pub creep_m2_per_year: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Landscape {
    pub columns: usize,
    pub rows: usize,
    pub wraps: bool,
    pub cell_metres: f64,
    pub height: Vec<f64>,
    pub fixed: Vec<bool>,
    pub routing: Routing,
}

#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct Routing {
    pub seed: u64,
    pub origin_cells: (i64, i64),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Settling {
    pub steps: usize,
    pub settled: bool,
    pub change_m_per_year: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Drainage {
    pub receiver: Vec<usize>,
    pub surface: Vec<f64>,
    pub area_m2: Vec<f64>,
    pub flow_m3_per_year: Vec<f64>,
}

impl Drainage {
    pub fn lake_depth(&self, landscape: &Landscape, i: usize) -> f64 {
        self.surface[i] - landscape.height[i]
    }
}

#[derive(PartialEq)]
struct Lowest(f64, usize);

impl Eq for Lowest {}

impl PartialOrd for Lowest {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Lowest {
    fn cmp(&self, other: &Self) -> Ordering {
        other.0.total_cmp(&self.0).then(other.1.cmp(&self.1))
    }
}

const NEIGHBOURS: [(i64, i64); 8] = [
    (-1, -1),
    (0, -1),
    (1, -1),
    (-1, 0),
    (1, 0),
    (-1, 1),
    (0, 1),
    (1, 1),
];

const SIDES: [(i64, i64); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];

impl Landscape {
    fn at(&self, x: i64, y: i64) -> Option<usize> {
        let (columns, rows) = (self.columns as i64, self.rows as i64);
        let x = if self.wraps { x.rem_euclid(columns) } else { x };
        ((0..columns).contains(&x) && (0..rows).contains(&y)).then(|| (y * columns + x) as usize)
    }

    fn neighbours(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (x, y) = ((i % self.columns) as i64, (i / self.columns) as i64);
        NEIGHBOURS.iter().filter_map(move |&(dx, dy)| {
            let distance = if dx != 0 && dy != 0 {
                core::f64::consts::SQRT_2
            } else {
                1.0
            };
            self.at(x + dx, y + dy)
                .map(|j| (j, distance * self.cell_metres))
        })
    }

    pub fn fill_hollows(&mut self) {
        self.height = self.flooded();
    }

    fn flooded(&self) -> Vec<f64> {
        let mut surface = self.height.clone();
        let mut done = vec![false; surface.len()];
        let mut open = BinaryHeap::new();
        for (i, &fixed) in self.fixed.iter().enumerate() {
            if fixed {
                done[i] = true;
                open.push(Lowest(surface[i], i));
            }
        }
        let rise = FLOOD_RISE_PER_METRE * self.cell_metres;
        while let Some(Lowest(level, i)) = open.pop() {
            for (j, _) in self.neighbours(i) {
                if !done[j] {
                    done[j] = true;
                    surface[j] = surface[j].max(level + rise);
                    open.push(Lowest(surface[j], j));
                }
            }
        }
        surface
    }

    fn rho(&self, i: usize) -> f64 {
        let (x, y) = ((i % self.columns) as i64, (i / self.columns) as i64);
        let (x, y) = (
            self.routing.origin_cells.0 + x,
            self.routing.origin_cells.1 + y,
        );
        unit(lattice(
            Seed(self.routing.seed),
            Purpose::Relief,
            RHO_LAYER,
            x,
            y,
            0,
        ))
    }

    pub fn drainage(&self, rain_m_per_year: &[f64]) -> Drainage {
        let surface = self.flooded();
        let cells = surface.len();
        let mut receiver: Vec<usize> = (0..cells).collect();
        for i in 0..cells {
            if self.fixed[i] {
                continue;
            }
            let diagonal = (2.0 - self.rho(i)) * self.cell_metres;
            let mut steepest = 0.0;
            for (j, distance) in self.neighbours(i) {
                let run = if distance > self.cell_metres {
                    diagonal
                } else {
                    distance
                };
                let slope = (surface[i] - surface[j]) / run;
                if slope > steepest {
                    steepest = slope;
                    receiver[i] = j;
                }
            }
        }
        let mut order: Vec<usize> = (0..cells).collect();
        order.sort_by(|&a, &b| surface[b].total_cmp(&surface[a]).then(a.cmp(&b)));
        let cell_area = self.cell_metres * self.cell_metres;
        let mut area_m2 = vec![cell_area; cells];
        let mut flow_m3_per_year: Vec<f64> = rain_m_per_year
            .iter()
            .map(|rain| rain * cell_area)
            .collect();
        for &i in &order {
            let r = receiver[i];
            if r != i {
                area_m2[r] += area_m2[i];
                flow_m3_per_year[r] += flow_m3_per_year[i];
            }
        }
        Drainage {
            receiver,
            surface,
            area_m2,
            flow_m3_per_year,
        }
    }

    pub fn evolve(&mut self, uplift_m_per_year: &[f64], laws: Laws, years: f64, step_years: f64) {
        let mut elapsed = 0.0;
        while elapsed < years {
            let dt = step_years.min(years - elapsed);
            self.step(uplift_m_per_year, laws, dt);
            elapsed += dt;
        }
    }

    pub fn settle(&mut self, uplift_m_per_year: &[f64], laws: Laws, step_years: f64) -> Settling {
        let free = self.fixed.iter().filter(|fixed| !**fixed).count().max(1) as f64;
        let mean_uplift = uplift_m_per_year
            .iter()
            .zip(&self.fixed)
            .filter(|(_, fixed)| !**fixed)
            .map(|(u, _)| u.abs())
            .sum::<f64>()
            / free;
        let settled = SETTLED_SHARE_OF_UPLIFT * mean_uplift;
        let mut best = f64::INFINITY;
        let mut since_best = 0;
        let mut steps = 0;
        let mut previous = f64::INFINITY;
        loop {
            let before = self.height.clone();
            self.step(uplift_m_per_year, laws, step_years);
            steps += 1;
            let change = self
                .height
                .iter()
                .zip(&before)
                .map(|(after, before)| (after - before).abs() / step_years)
                .sum::<f64>()
                / free;
            let ratio = change / previous;
            previous = change;
            let still_to_come = if ratio < 1.0 {
                change * ratio / (1.0 - ratio)
            } else {
                f64::INFINITY
            };
            let settling = |settled| Settling {
                steps,
                settled,
                change_m_per_year: change,
            };
            if change <= settled && still_to_come <= settled {
                return settling(true);
            }
            if change < best {
                best = change;
                since_best = 0;
            } else {
                since_best += 1;
                if since_best >= STEPS_WITHOUT_PROGRESS {
                    return settling(false);
                }
            }
        }
    }

    fn step(&mut self, uplift_m_per_year: &[f64], laws: Laws, dt: f64) {
        let before = self.height.clone();
        for ((height, &fixed), &uplift) in self
            .height
            .iter_mut()
            .zip(&self.fixed)
            .zip(uplift_m_per_year)
        {
            if !fixed {
                *height += uplift * dt;
            }
        }
        let rain = vec![1.0; self.height.len()];
        self.incise(&self.drainage(&rain), laws, dt);
        self.creep(laws, dt, before);
    }

    fn incise(&mut self, drainage: &Drainage, laws: Laws, dt: f64) {
        let mut order: Vec<usize> = (0..self.height.len()).collect();
        order.sort_by(|&a, &b| {
            drainage.surface[a]
                .total_cmp(&drainage.surface[b])
                .then(a.cmp(&b))
        });
        let lake: Vec<bool> = drainage
            .surface
            .iter()
            .zip(&self.height)
            .map(|(surface, height)| surface > height)
            .collect();
        for &i in &order {
            let r = drainage.receiver[i];
            if self.fixed[i] || r == i || lake[i] {
                continue;
            }
            let distance = self
                .neighbours(i)
                .find(|&(j, _)| j == r)
                .map_or(self.cell_metres, |(_, d)| d);
            let base = if lake[r] {
                drainage.surface[r]
            } else {
                self.height[r]
            };
            let factor =
                laws.erodibility * dt * libm::pow(drainage.area_m2[i], AREA_EXPONENT) / distance;
            let lowered = (self.height[i] + factor * base) / (1.0 + factor);
            self.height[i] = lowered.max(base);
        }
    }

    fn creep(&mut self, laws: Laws, dt: f64, guess: Vec<f64>) {
        let rate = laws.creep_m2_per_year * dt / (self.cell_metres * self.cell_metres);
        let start = core::mem::replace(&mut self.height, guess);
        for _ in 0..CREEP_SWEEPS {
            for (i, &own) in start.iter().enumerate() {
                if self.fixed[i] {
                    continue;
                }
                let (x, y) = ((i % self.columns) as i64, (i / self.columns) as i64);
                let mut sum = 0.0;
                let mut count = 0.0;
                for (dx, dy) in SIDES {
                    if let Some(j) = self.at(x + dx, y + dy) {
                        sum += self.height[j];
                        count += 1.0;
                    }
                }
                self.height[i] = (own + rate * sum) / (1.0 + rate * count);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ramp(columns: usize, rows: usize, cell_metres: f64) -> Landscape {
        let mut height = vec![0.0; columns * rows];
        let mut fixed = vec![false; columns * rows];
        for y in 0..rows {
            for x in 0..columns {
                height[y * columns + x] = x as f64 * 0.01 + ((x * 7 + y * 13) % 5) as f64 * 1e-4;
                fixed[y * columns + x] = x == 0;
            }
        }
        Landscape {
            columns,
            rows,
            wraps: false,
            cell_metres,
            height,
            fixed,
            routing: Routing::default(),
        }
    }

    fn channel_slopes_match_stream_power(land: &Landscape, uplift: f64, laws: Laws) -> usize {
        let drainage = land.drainage(&vec![1.0; land.height.len()]);
        let cell_area = land.cell_metres * land.cell_metres;
        let mut checked = 0;
        for i in 0..land.height.len() {
            let r = drainage.receiver[i];
            if land.fixed[i] || r == i || drainage.area_m2[i] < 4.0 * cell_area {
                continue;
            }
            let distance = land
                .neighbours(i)
                .find(|&(j, _)| j == r)
                .map_or(land.cell_metres, |(_, d)| d);
            let slope = (land.height[i] - land.height[r]) / distance;
            let expected =
                uplift / laws.erodibility * libm::pow(drainage.area_m2[i], -AREA_EXPONENT);
            assert!(
                (slope - expected).abs() < 0.05 * expected,
                "{slope} vs {expected}"
            );
            checked += 1;
        }
        checked
    }

    #[test]
    fn channels_settle_on_the_stream_power_slope() {
        let laws = Laws {
            erodibility: 1e-5,
            creep_m2_per_year: 1e-12,
        };
        let uplift = 1e-3;
        let mut land = ramp(24, 24, 100.0);
        let rates = vec![uplift; land.height.len()];
        land.settle(&rates, laws, 2e4);
        assert!(channel_slopes_match_stream_power(&land, uplift, laws) > 10);
    }

    #[test]
    fn a_wrapping_strip_drains_across_its_seam() {
        let mut land = ramp(12, 4, 100.0);
        for y in 0..4 {
            for x in 0..12 {
                land.height[y * 12 + x] = if x == 1 { 5.0 } else { 0.1 * x as f64 };
                land.fixed[y * 12 + x] = false;
            }
        }
        land.fixed[0] = true;
        land.height[0] = 0.0;
        land.wraps = true;
        let drainage = land.drainage(&vec![1.0; land.height.len()]);
        assert_eq!(drainage.receiver[11], 0);
    }

    #[test]
    fn water_on_a_tilted_plane_follows_the_true_slope_not_the_grid() {
        let size = 200;
        let heading = 30f64.to_radians();
        let (gx, gy) = (libm::cos(heading), libm::sin(heading));
        let mut height = vec![0.0; size * size];
        for y in 0..size {
            for x in 0..size {
                height[y * size + x] = -(gx * x as f64 + gy * y as f64);
            }
        }
        let land = Landscape {
            columns: size,
            rows: size,
            wraps: false,
            cell_metres: 1.0,
            height,
            fixed: (0..size * size)
                .map(|i| i % size == size - 1 || i / size == size - 1)
                .collect(),
            routing: Routing {
                seed: 11,
                origin_cells: (0, 0),
            },
        };
        let drainage = land.drainage(&vec![1.0; size * size]);
        let mut errors = Vec::new();
        for start in [(5, 5), (5, 60), (60, 5), (30, 30)] {
            let mut at = start.1 * size + start.0;
            for _ in 0..size * 2 {
                let next = drainage.receiver[at];
                if next == at {
                    break;
                }
                at = next;
            }
            let (dx, dy) = (
                (at % size) as f64 - start.0 as f64,
                (at / size) as f64 - start.1 as f64,
            );
            errors.push((libm::atan2(dy, dx) - heading).to_degrees().abs());
        }
        let mean = errors.iter().sum::<f64>() / errors.len() as f64;
        assert!(
            mean < 5.0,
            "flow strays {mean:.1} degrees from the slope: {errors:?}"
        );
    }

    #[test]
    fn every_drop_of_rain_leaves_through_the_fixed_edge() {
        let land = ramp(16, 16, 50.0);
        let rain = vec![0.8; land.height.len()];
        let drainage = land.drainage(&rain);
        let fallen: f64 = rain.iter().map(|r| r * 50.0 * 50.0).sum();
        let left: f64 = (0..land.height.len())
            .filter(|&i| land.fixed[i])
            .map(|i| drainage.flow_m3_per_year[i])
            .sum();
        assert!((left - fallen).abs() < 1e-6 * fallen, "{left} vs {fallen}");
    }

    #[test]
    fn a_hollow_fills_into_a_lake_that_spills_toward_the_outlet() {
        let size = 9;
        let mut land = Landscape {
            columns: size,
            rows: size,
            wraps: false,
            cell_metres: 10.0,
            height: vec![5.0; size * size],
            fixed: vec![false; size * size],
            routing: Routing::default(),
        };
        land.fixed[4] = true;
        land.height[4] = 0.0;
        land.height[4 * size + 4] = 1.0;
        let drainage = land.drainage(&vec![1.0; size * size]);
        assert!(drainage.lake_depth(&land, 4 * size + 4) > 3.9);
        let mut at = 4 * size + 4;
        for _ in 0..size * size {
            if land.fixed[at] {
                break;
            }
            at = drainage.receiver[at];
        }
        assert_eq!(at, 4);
    }

    #[test]
    fn creep_settles_a_hill_into_the_uplift_parabola() {
        let size = 21;
        let cell = 1.0;
        let mut land = Landscape {
            columns: size,
            rows: 1,
            wraps: false,
            cell_metres: cell,
            height: vec![0.0; size],
            fixed: (0..size).map(|x| x == 0 || x == size - 1).collect(),
            routing: Routing::default(),
        };
        let laws = Laws {
            erodibility: 0.0,
            creep_m2_per_year: 0.01,
        };
        let uplift = 1e-4;
        land.settle(&vec![uplift; size], laws, 1e3);
        let half = (size - 1) as f64 / 2.0 * cell;
        let crest = uplift * half * half / (2.0 * laws.creep_m2_per_year);
        let middle = land.height[size / 2];
        assert!((middle - crest).abs() < 0.02 * crest, "{middle} vs {crest}");
    }
}
