// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec;
use alloc::vec::Vec;
use core::f64::consts::TAU;
use serde::{Deserialize, Serialize};

use crate::Seed;
use crate::atlas::Atlas;
use crate::erosion::{Drainage, FLATTEST_SLOPE, Landscape, Laws, Routing, Settling};
use crate::hash::{Purpose, lattice, unit};

pub const PLANET_REFINEMENT: usize = 16;
pub const REGION_CELLS: usize = 256;
pub const REGION_SPACING_CELLS: usize = REGION_CELLS / 2;
const BLUR_PASSES: usize = 2;
const SECONDS_PER_YEAR: f64 = 365.25 * 86_400.0;
const RIVER_WIDTH_PER_ROOT_DISCHARGE: f64 = 7.2;
const MEANDER_WAVELENGTH_PER_WIDTH: f64 = 10.9;
const MEANDER_AMPLITUDE_PER_WIDTH: f64 = 2.7;
const MEANDER_SPREAD: f64 = 0.5;
const MEANDER_STREAM: i64 = 1;
const ROUTING_STREAM: i64 = 2;
const COURSE_SMOOTHING_CELLS: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Rock {
    pub laws: Laws,
    pub relief_erosion_per_year: f64,
}

impl Rock {
    fn uplift(&self, metres_above_sea: f64) -> f64 {
        self.relief_erosion_per_year * metres_above_sea.max(0.0)
    }

    fn step_years(&self) -> f64 {
        1.0 / self.laws.erodibility
    }

    fn channel_head_m2(&self, cell_metres: f64) -> f64 {
        let root = self.laws.creep_m2_per_year / (self.laws.erodibility * cell_metres);
        root * root
    }

    fn cell_response_years(&self, cell_metres: f64) -> f64 {
        let creep = cell_metres * cell_metres / self.laws.creep_m2_per_year;
        self.step_years().min(creep)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Grid {
    columns: usize,
    rows: usize,
    wraps: bool,
    cell_metres: f64,
    origin_metres: (f64, f64),
}

impl Grid {
    fn centre(&self, column: usize, row: usize) -> (f64, f64) {
        (
            self.origin_metres.0 + (column as f64 + 0.5) * self.cell_metres,
            self.origin_metres.1 + (row as f64 + 0.5) * self.cell_metres,
        )
    }

    fn edge(&self, column: usize, row: usize) -> bool {
        !self.wraps
            && (column == 0 || column + 1 == self.columns || row == 0 || row + 1 == self.rows)
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Field {
    pub land: Landscape,
    pub origin_metres: (f64, f64),
    pub settling: Settling,
}

impl Field {
    pub fn metres(&self, x_metres: f64, y_metres: f64) -> f64 {
        let land = &self.land;
        let fx = (x_metres - self.origin_metres.0) / land.cell_metres - 0.5;
        let fy = (y_metres - self.origin_metres.1) / land.cell_metres - 0.5;
        let (column, row) = (libm::floor(fx), libm::floor(fy));
        let (tx, ty) = (fx - column, fy - row);
        let (column, row) = (column as i64, row as i64);
        let (columns, rows) = (land.columns as i64, land.rows as i64);
        let at = |x: i64, y: i64| {
            let x = if land.wraps {
                x.rem_euclid(columns)
            } else {
                x.clamp(0, columns - 1)
            };
            land.height[(y.clamp(0, rows - 1) * columns + x) as usize]
        };
        let line = |y: i64| {
            catmull_rom(
                [
                    at(column - 1, y),
                    at(column, y),
                    at(column + 1, y),
                    at(column + 2, y),
                ],
                tx,
            )
        };
        catmull_rom([line(row - 1), line(row), line(row + 1), line(row + 2)], ty)
    }

    pub fn index(&self, x_metres: f64, y_metres: f64) -> Option<usize> {
        let land = &self.land;
        let column = libm::floor((x_metres - self.origin_metres.0) / land.cell_metres) as i64;
        let row = libm::floor((y_metres - self.origin_metres.1) / land.cell_metres) as i64;
        ((0..land.columns as i64).contains(&column) && (0..land.rows as i64).contains(&row))
            .then(|| row as usize * land.columns + column as usize)
    }

    pub fn centre(&self, i: usize) -> (f64, f64) {
        let land = &self.land;
        (
            self.origin_metres.0 + ((i % land.columns) as f64 + 0.5) * land.cell_metres,
            self.origin_metres.1 + ((i / land.columns) as f64 + 0.5) * land.cell_metres,
        )
    }
}

fn catmull_rom(p: [f64; 4], t: f64) -> f64 {
    let a = -0.5 * p[0] + 1.5 * p[1] - 1.5 * p[2] + 0.5 * p[3];
    let b = p[0] - 2.5 * p[1] + 2.0 * p[2] - 0.5 * p[3];
    let c = -0.5 * p[0] + 0.5 * p[2];
    ((a * t + b) * t + c) * t + p[1]
}

#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    pub field: Field,
    pub receiver: Vec<usize>,
    pub area_m2: Vec<f64>,
    main_donor: Vec<Option<usize>>,
    wetness: Vec<f64>,
    course: Vec<(f64, f64)>,
    phase: Vec<f64>,
    wander: Vec<f64>,
    reach_cells: i64,
    channel_head_m2: f64,
    narrowest_metres: f64,
}

impl Region {
    fn new(
        field: Field,
        rain_m_per_year: f64,
        rock: Rock,
        narrowest_metres: f64,
        seed: Seed,
        layer: u64,
    ) -> Self {
        let channel_head_m2 = rock.channel_head_m2(field.land.cell_metres);
        let drainage = field.land.drainage(&vec![1.0; field.land.height.len()]);
        let cells = drainage.receiver.len();
        let land = &field.land;
        let first = (
            libm::round(field.origin_metres.0 / land.cell_metres) as i64,
            libm::round(field.origin_metres.1 / land.cell_metres) as i64,
        );
        let wander: Vec<f64> = (0..cells)
            .map(|i| {
                let h = lattice(
                    seed,
                    Purpose::Relief,
                    layer,
                    first.0 + (i % land.columns) as i64,
                    first.1 + (i / land.columns) as i64,
                    MEANDER_STREAM,
                );
                1.0 + MEANDER_SPREAD * (2.0 * unit(h) - 1.0)
            })
            .collect();
        let drawn = |area: f64| river_width_metres(area, rain_m_per_year).max(narrowest_metres);
        let mut main_donor: Vec<Option<usize>> = vec![None; cells];
        for i in 0..cells {
            let r = drainage.receiver[i];
            if r != i && main_donor[r].is_none_or(|d| drainage.area_m2[i] > drainage.area_m2[d]) {
                main_donor[r] = Some(i);
            }
        }
        let mut order: Vec<usize> = (0..cells).collect();
        order.sort_by(|&a, &b| {
            drainage.surface[a]
                .total_cmp(&drainage.surface[b])
                .then(a.cmp(&b))
        });
        let mut phase = vec![0.0; cells];
        for &i in &order {
            let r = drainage.receiver[i];
            if r != i {
                let length = distance(field.centre(i), field.centre(r));
                let wavelength =
                    MEANDER_WAVELENGTH_PER_WIDTH * drawn(drainage.area_m2[i]) * wander[i];
                phase[i] = phase[r] + TAU * length / wavelength;
            }
        }
        let widest = drawn(drainage.area_m2.iter().copied().fold(0.0, f64::max));
        let reach = widest / 2.0 + MEANDER_AMPLITUDE_PER_WIDTH * widest * (1.0 + MEANDER_SPREAD);
        let reach_cells =
            libm::ceil(reach / field.land.cell_metres) as i64 + 1 + COURSE_SMOOTHING_CELLS as i64;
        let course = (0..cells)
            .map(|i| {
                let mut sum = field.centre(i);
                let mut count = 1.0;
                let mut add = |j: usize| {
                    let (x, y) = field.centre(j);
                    sum = (sum.0 + x, sum.1 + y);
                    count += 1.0;
                };
                let (mut down, mut up) = (i, i);
                for _ in 0..COURSE_SMOOTHING_CELLS {
                    let next = drainage.receiver[down];
                    if next != down {
                        down = next;
                        add(down);
                    }
                    if let Some(donor) = main_donor[up] {
                        up = donor;
                        add(up);
                    }
                }
                (sum.0 / count, sum.1 / count)
            })
            .collect();
        let wet = wetness(&field.land, &drainage);
        Self {
            field,
            receiver: drainage.receiver,
            wetness: wet,
            area_m2: drainage.area_m2,
            main_donor,
            course,
            phase,
            wander,
            reach_cells,
            channel_head_m2,
            narrowest_metres,
        }
    }

    fn channel(&self, i: usize) -> bool {
        self.receiver[i] != i
            && self.field.land.height[i] >= 0.0
            && self.area_m2[i] >= self.channel_head_m2
    }

    fn course_distance(&self, i: usize, point: (f64, f64), rain_m_per_year: f64) -> f64 {
        let r = self.receiver[i];
        let (a, b) = (self.course[i], self.course[r]);
        let before =
            self.main_donor[i].map_or((2.0 * a.0 - b.0, 2.0 * a.1 - b.1), |d| self.course[d]);
        let after = if self.receiver[r] == r {
            (2.0 * b.0 - a.0, 2.0 * b.1 - a.1)
        } else {
            self.course[self.receiver[r]]
        };
        let width = river_width_metres(self.area_m2[i], rain_m_per_year).max(self.narrowest_metres);
        let amplitude = MEANDER_AMPLITUDE_PER_WIDTH * width;
        let samples = libm::ceil(2.0 * distance(a, b) / self.narrowest_metres).max(1.0) as usize;
        let mut nearest = f64::INFINITY;
        let mut previous: Option<(f64, f64)> = None;
        for k in 0..=samples {
            let t = k as f64 / samples as f64;
            let at = |t: f64| {
                (
                    catmull_rom([before.0, a.0, b.0, after.0], t),
                    catmull_rom([before.1, a.1, b.1, after.1], t),
                )
            };
            let (x, y) = at(t);
            let ahead = at((t + 1.0 / samples as f64).min(1.0));
            let behind = at((t - 1.0 / samples as f64).max(0.0));
            let (tx, ty) = (ahead.0 - behind.0, ahead.1 - behind.1);
            let length = libm::hypot(tx, ty).max(f64::MIN_POSITIVE);
            let phase = self.phase[i] + t * (self.phase[r] - self.phase[i]);
            let wander = self.wander[i] + t * (self.wander[r] - self.wander[i]);
            let swing = amplitude * wander * libm::sin(phase);
            let bend = (x - ty / length * swing, y + tx / length * swing);
            if let Some(last) = previous {
                nearest = nearest.min(distance_to_segment(point, last, bend));
            }
            previous = Some(bend);
        }
        nearest - width / 2.0
    }

    pub fn wetness_at(&self, x_metres: f64, y_metres: f64) -> f64 {
        let land = &self.field.land;
        let fx = (x_metres - self.field.origin_metres.0) / land.cell_metres - 0.5;
        let fy = (y_metres - self.field.origin_metres.1) / land.cell_metres - 0.5;
        let (column, row) = (libm::floor(fx), libm::floor(fy));
        let (tx, ty) = (fx - column, fy - row);
        let (column, row) = (column as i64, row as i64);
        let at = |x: i64, y: i64| {
            let x = x.clamp(0, land.columns as i64 - 1) as usize;
            let y = y.clamp(0, land.rows as i64 - 1) as usize;
            libm::log(self.wetness[y * land.columns + x])
        };
        libm::exp(
            (1.0 - ty) * ((1.0 - tx) * at(column, row) + tx * at(column + 1, row))
                + ty * ((1.0 - tx) * at(column, row + 1) + tx * at(column + 1, row + 1)),
        )
    }

    pub fn river_at(&self, x_metres: f64, y_metres: f64, rain_m_per_year: f64) -> bool {
        let Some(here) = self.field.index(x_metres, y_metres) else {
            return false;
        };
        let land = &self.field.land;
        let (column, row) = ((here % land.columns) as i64, (here / land.columns) as i64);
        for dy in -self.reach_cells..=self.reach_cells {
            for dx in -self.reach_cells..=self.reach_cells {
                let (c, r) = (column + dx, row + dy);
                if !(0..land.columns as i64).contains(&c) || !(0..land.rows as i64).contains(&r) {
                    continue;
                }
                let i = r as usize * land.columns + c as usize;
                if self.channel(i)
                    && self.course_distance(i, (x_metres, y_metres), rain_m_per_year) < 0.0
                {
                    return true;
                }
            }
        }
        false
    }
}

fn distance(a: (f64, f64), b: (f64, f64)) -> f64 {
    libm::hypot(b.0 - a.0, b.1 - a.1)
}

pub fn river_width_metres(area_m2: f64, rain_m_per_year: f64) -> f64 {
    let discharge = area_m2 * rain_m_per_year.max(0.0) / SECONDS_PER_YEAR;
    RIVER_WIDTH_PER_ROOT_DISCHARGE * libm::sqrt(discharge)
}

fn distance_to_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (abx, aby) = (b.0 - a.0, b.1 - a.1);
    let length = abx * abx + aby * aby;
    let t = if length > 0.0 {
        (((p.0 - a.0) * abx + (p.1 - a.1) * aby) / length).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let (dx, dy) = (p.0 - a.0 - t * abx, p.1 - a.1 - t * aby);
    libm::sqrt(dx * dx + dy * dy)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sample {
    pub metres: f64,
    pub east_slope: f64,
    pub south_slope: f64,
    pub river: bool,
    pub wetness: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Local {
    pub river: bool,
    pub wetness: f64,
}

pub fn wetness(land: &Landscape, drainage: &Drainage) -> Vec<f64> {
    let index: Vec<f64> = (0..land.height.len())
        .map(|i| {
            let r = drainage.receiver[i];
            let run = if r == i {
                land.cell_metres
            } else {
                distance_between(land, i, r)
            };
            let slope = ((drainage.surface[i] - drainage.surface[r]) / run).max(FLATTEST_SLOPE);
            libm::log(drainage.area_m2[i] / land.cell_metres / slope)
        })
        .collect();
    let counted: Vec<f64> = index
        .iter()
        .zip(&land.height)
        .filter(|(_, height)| **height >= 0.0)
        .map(|(value, _)| *value)
        .collect();
    let mean = if counted.is_empty() {
        0.0
    } else {
        counted.iter().sum::<f64>() / counted.len() as f64
    };
    index.iter().map(|value| libm::exp(value - mean)).collect()
}

fn distance_between(land: &Landscape, a: usize, b: usize) -> f64 {
    let (ax, ay) = ((a % land.columns) as f64, (a / land.columns) as f64);
    let (bx, by) = ((b % land.columns) as f64, (b / land.columns) as f64);
    let dx = (ax - bx).abs().min(land.columns as f64 - (ax - bx).abs());
    libm::hypot(dx, ay - by) * land.cell_metres
}

pub trait Terrain {
    fn sample(&self, x_metres: f64, y_metres: f64) -> Sample;
    fn coarse_metres(&self, x_metres: f64, y_metres: f64) -> f64;
    fn prepare(&self, _x_metres: f64, _y_metres: f64) {}
}

impl Terrain for Atlas {
    fn sample(&self, x_metres: f64, y_metres: f64) -> Sample {
        let at = |x: f64, y: f64| self.blend(x, y).metres;
        let h = self.cell_metres / 2.0;
        Sample {
            metres: at(x_metres, y_metres),
            east_slope: (at(x_metres + h, y_metres) - at(x_metres - h, y_metres)) / (2.0 * h),
            south_slope: (at(x_metres, y_metres + h) - at(x_metres, y_metres - h)) / (2.0 * h),
            river: false,
            wetness: 1.0,
        }
    }

    fn coarse_metres(&self, x_metres: f64, y_metres: f64) -> f64 {
        self.blend(x_metres, y_metres).metres
    }
}

pub fn blended_metres(
    x_metres: f64,
    y_metres: f64,
    spacing: f64,
    around_metres: f64,
    region: &dyn Fn((i64, i64)) -> alloc::sync::Arc<Region>,
) -> f64 {
    let x_metres = wrapped(x_metres, around_metres);
    let around = libm::round(around_metres / spacing) as i64;
    covering(x_metres, y_metres, spacing)
        .iter()
        .filter(|(_, weight)| *weight > 0.0)
        .map(|&((i, j), weight)| {
            let wrapped = i.rem_euclid(around);
            let shift = (i - wrapped) as f64 * spacing;
            weight
                * region((wrapped, j))
                    .field
                    .metres(x_metres - shift, y_metres)
        })
        .sum()
}

pub fn local_at(
    x_metres: f64,
    y_metres: f64,
    spacing: f64,
    around_metres: f64,
    rain_m_per_year: f64,
    region: &dyn Fn((i64, i64)) -> alloc::sync::Arc<Region>,
) -> Local {
    let x_metres = wrapped(x_metres, around_metres);
    let around = libm::round(around_metres / spacing) as i64;
    let mut wetness = 0.0;
    let mut strongest = ((0, 0), 0.0, f64::NEG_INFINITY);
    for ((i, j), weight) in covering(x_metres, y_metres, spacing) {
        let wrapped = i.rem_euclid(around);
        let shift = (i - wrapped) as f64 * spacing;
        if weight > 0.0 {
            wetness += weight * region((wrapped, j)).wetness_at(x_metres - shift, y_metres);
        }
        if weight > strongest.2 {
            strongest = ((wrapped, j), shift, weight);
        }
    }
    let (index, shift, _) = strongest;
    Local {
        river: region(index).river_at(x_metres - shift, y_metres, rain_m_per_year),
        wetness,
    }
}

pub fn wrapped(x_metres: f64, around_metres: f64) -> f64 {
    x_metres - libm::floor(x_metres / around_metres) * around_metres
}

pub fn spacing_metres(cell_metres: f64) -> f64 {
    cell_metres * REGION_SPACING_CELLS as f64
}

pub fn covering(x_metres: f64, y_metres: f64, spacing: f64) -> [((i64, i64), f64); 4] {
    let (fx, fy) = (x_metres / spacing, y_metres / spacing);
    let (i, j) = (libm::floor(fx), libm::floor(fy));
    let (tx, ty) = (fx - i, fy - j);
    let (i, j) = (i as i64, j as i64);
    [
        ((i, j), (1.0 - tx) * (1.0 - ty)),
        ((i + 1, j), tx * (1.0 - ty)),
        ((i, j + 1), (1.0 - tx) * ty),
        ((i + 1, j + 1), tx * ty),
    ]
}

fn settled(
    seed: Seed,
    layer: u64,
    grid: Grid,
    refinement: usize,
    parent: impl Fn(f64, f64) -> f64,
    rock: Rock,
) -> Field {
    let cells = grid.columns * grid.rows;
    let response_years = rock.cell_response_years(grid.cell_metres);
    let first = (
        libm::round(grid.origin_metres.0 / grid.cell_metres) as i64,
        libm::round(grid.origin_metres.1 / grid.cell_metres) as i64,
    );
    let mut height = Vec::with_capacity(cells);
    let mut fixed = Vec::with_capacity(cells);
    let mut uplift = Vec::with_capacity(cells);
    let mut jiggle = Vec::with_capacity(cells);
    for row in 0..grid.rows {
        for column in 0..grid.columns {
            let (x, y) = grid.centre(column, row);
            let below = parent(x, y);
            let pinned = grid.edge(column, row) || below < 0.0;
            let h = lattice(
                seed,
                Purpose::Relief,
                layer,
                first.0 + column as i64,
                first.1 + row as i64,
                0,
            );
            let rise = rock.uplift(below);
            height.push(below);
            if !pinned {
                jiggle.push((height.len() - 1, rise * response_years * unit(h)));
            }
            fixed.push(pinned);
            uplift.push(rise);
        }
    }
    let mut land = Landscape {
        columns: grid.columns,
        rows: grid.rows,
        wraps: grid.wraps,
        cell_metres: grid.cell_metres,
        height,
        fixed,
        routing: Routing {
            seed: lattice(seed, Purpose::Relief, layer, 0, 0, ROUTING_STREAM),
            origin_cells: first,
        },
    };
    let parent_heights = land.height.clone();
    for (i, rise) in jiggle {
        land.height[i] += rise;
    }
    let settling = land.settle(&uplift, rock.laws, response_years);
    let smooth = blurred(&land.height, grid, refinement / 2);
    for i in 0..cells {
        if !land.fixed[i] {
            land.height[i] = parent_heights[i] + land.height[i] - smooth[i];
        }
    }
    land.fill_hollows();
    Field {
        land,
        origin_metres: grid.origin_metres,
        settling,
    }
}

pub fn planet(seed: Seed, atlas: &Atlas, rock: Rock) -> Field {
    let grid = Grid {
        columns: atlas.width as usize * PLANET_REFINEMENT,
        rows: atlas.height as usize * PLANET_REFINEMENT,
        wraps: true,
        cell_metres: atlas.cell_metres / PLANET_REFINEMENT as f64,
        origin_metres: (0.0, 0.0),
    };
    settled(
        seed,
        0,
        grid,
        PLANET_REFINEMENT,
        |x, y| atlas.blend(x, y).metres,
        rock,
    )
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Level {
    pub layer: u64,
    pub parent_cell_metres: f64,
    pub refinement: usize,
}

impl Level {
    pub fn cell_metres(&self) -> f64 {
        self.parent_cell_metres / self.refinement as f64
    }

    pub fn spacing_metres(&self) -> f64 {
        spacing_metres(self.cell_metres())
    }
}

pub fn region(
    seed: Seed,
    level: Level,
    index: (i64, i64),
    parent: impl Fn(f64, f64) -> f64,
    rain_m_per_year: f64,
    narrowest_metres: f64,
    rock: Rock,
) -> Region {
    let spacing = level.spacing_metres();
    let grid = Grid {
        columns: REGION_CELLS,
        rows: REGION_CELLS,
        wraps: false,
        cell_metres: level.cell_metres(),
        origin_metres: (
            (index.0 - 1) as f64 * spacing,
            (index.1 - 1) as f64 * spacing,
        ),
    };
    let field = settled(seed, level.layer, grid, level.refinement, parent, rock);
    Region::new(
        field,
        rain_m_per_year,
        rock,
        narrowest_metres,
        seed,
        level.layer,
    )
}

fn blurred(values: &[f64], grid: Grid, radius: usize) -> Vec<f64> {
    let mut smooth = values.to_vec();
    for _ in 0..BLUR_PASSES {
        for (along, across, step_along, step_across, wraps) in [
            (grid.columns, grid.rows, 1, grid.columns, grid.wraps),
            (grid.rows, grid.columns, grid.columns, 1, false),
        ] {
            let before = smooth.clone();
            for line in 0..across {
                for at in 0..along {
                    let mut sum = 0.0;
                    for offset in -(radius as i64)..=radius as i64 {
                        let k = at as i64 + offset;
                        let k = if wraps {
                            k.rem_euclid(along as i64)
                        } else {
                            k.clamp(0, along as i64 - 1)
                        } as usize;
                        sum += before[line * step_across + k * step_along];
                    }
                    smooth[line * step_across + at * step_along] = sum / (2 * radius + 1) as f64;
                }
            }
        }
    }
    smooth
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_valley_floor_is_wetter_than_the_ridge_above_it() {
        let size = 21;
        let middle = size / 2;
        let land = Landscape {
            columns: size,
            rows: size,
            wraps: false,
            cell_metres: 10.0,
            height: (0..size * size)
                .map(|i| {
                    let (x, y) = (i % size, i / size);
                    x.abs_diff(middle) as f64 * 2.0 + (size - y) as f64 * 0.5
                })
                .collect(),
            fixed: (0..size * size).map(|i| i / size == size - 1).collect(),
            routing: Routing::default(),
        };
        let drainage = land.drainage(&vec![1.0; size * size]);
        let wet = wetness(&land, &drainage);
        let row = (size - 4) * size;
        assert!(wet[row + middle] > 10.0 * wet[row]);
        assert!(wet[row + middle] > 1.0 && wet[row] < 1.0);
    }

    #[test]
    fn eroded_detail_keeps_the_parent_height_at_the_parent_scale() {
        let rock = Rock {
            laws: Laws {
                erodibility: 1e-5,
                creep_m2_per_year: 0.01,
            },
            relief_erosion_per_year: 1.535e-7,
        };
        let refinement = 8;
        let grid = Grid {
            columns: 48,
            rows: 48,
            wraps: false,
            cell_metres: 100.0,
            origin_metres: (0.0, 0.0),
        };
        let parent = |x: f64, y: f64| 200.0 + 0.2 * x + 0.1 * y;
        let field = settled(Seed(7), 1, grid, refinement, parent, rock);
        let land = &field.land;
        let block = refinement;
        for by in 1..grid.rows / block - 1 {
            for bx in 1..grid.columns / block - 1 {
                let (mut child, mut wanted) = (0.0, 0.0);
                for y in by * block..(by + 1) * block {
                    for x in bx * block..(bx + 1) * block {
                        let (cx, cy) = grid.centre(x, y);
                        child += land.height[y * grid.columns + x];
                        wanted += parent(cx, cy);
                    }
                }
                let cells = (block * block) as f64;
                let (child, wanted) = (child / cells, wanted / cells);
                assert!(
                    (child - wanted).abs() < 0.02 * wanted,
                    "{child} vs {wanted}"
                );
            }
        }
    }
}
