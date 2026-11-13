// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec;
use alloc::vec::Vec;
use core::f64::consts::FRAC_PI_2;

use crate::material::{Material, Phase};
use crate::wind::{self, Belts, Winds};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circulation {
    pub gas_spread: f64,
    pub wind_carry: f64,
    pub age_steps: u32,
    pub belts: Belts,
    pub winds: Winds,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sea {
    Reservoir,
    SinksThrough,
}

pub struct Substance<'a> {
    pub material: &'a Material,
    pub sea: Sea,
    pub soil_capacity: f64,
    pub evaporation_limit: Option<f64>,
}

pub struct Grid<'a> {
    pub width: u32,
    pub height: u32,
    pub elevation: &'a [f64],
    pub celsius: &'a [f64],
    pub condensation_celsius: &'a [f64],
    pub ocean: &'a [bool],
    pub sun_latitude: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Pools {
    pub vapor: Vec<f64>,
    pub condensed: Vec<f64>,
    pub soil: Vec<f64>,
}

impl Pools {
    pub fn empty(cells: usize) -> Self {
        Self::with_condensed(vec![0.0; cells])
    }

    pub fn humid(grid: &Grid, material: &Material, soil_capacity: f64) -> Self {
        let mut pools = Self::empty(grid.len());
        for i in 0..grid.len() {
            pools.vapor[i] = material.saturation(grid.celsius[i]);
            pools.soil[i] = if grid.ocean[i] { 0.0 } else { soil_capacity };
        }
        pools
    }

    pub fn with_condensed(condensed: Vec<f64>) -> Self {
        let cells = condensed.len();
        Self {
            vapor: vec![0.0; cells],
            condensed,
            soil: vec![0.0; cells],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Flux {
    pub flow: f64,
    pub condensed: f64,
    pub vented: f64,
    pub soaked: f64,
    pub evaporated: f64,
}

pub struct Settled {
    pub pools: Pools,
    pub flux: Vec<Flux>,
    pub moved_per_step: Vec<f64>,
}

impl Grid<'_> {
    pub(crate) fn len(&self) -> usize {
        self.elevation.len()
    }

    pub(crate) fn latitude(&self, row: u32) -> f64 {
        (f64::from(row) + 0.5) / f64::from(self.height) * 2.0 - 1.0
    }

    #[inline]
    pub(crate) fn neighbors(&self, i: usize) -> ([usize; 8], usize) {
        let cell = i as u32;
        self.neighbors_at((cell % self.width) as usize, (cell / self.width) as usize)
    }

    #[inline]
    pub(crate) fn neighbors_at(&self, x: usize, y: usize) -> ([usize; 8], usize) {
        let width = self.width as usize;
        let columns = [
            if x == 0 { width - 1 } else { x - 1 },
            x,
            if x + 1 == width { 0 } else { x + 1 },
        ];
        let mut out = [0; 8];
        let mut count = 0;
        for ny in [y.wrapping_sub(1), y, y + 1] {
            if ny >= self.height as usize {
                continue;
            }
            for (k, &nx) in columns.iter().enumerate() {
                if ny == y && k == 1 {
                    continue;
                }
                out[count] = ny * width + nx;
                count += 1;
            }
        }
        (out, count)
    }

    pub(crate) fn downhill(&self, i: usize) -> Option<usize> {
        let (neighbors, count) = self.neighbors(i);
        neighbors[..count]
            .iter()
            .copied()
            .filter(|&j| self.elevation[j] < self.elevation[i])
            .min_by(|&a, &b| self.elevation[a].total_cmp(&self.elevation[b]))
    }

    fn downwind_row(&self, row: u32, circulation: &Circulation) -> (usize, i64, f64) {
        let belts = &circulation.belts;
        let wind = wind::prevailing(self.latitude(row) - self.sun_latitude, belts);
        let target = i64::from(row) + wind.south;
        if !(0..i64::from(self.height)).contains(&target) {
            return (row as usize, wind.east, 0.0);
        }
        let south = wind::boundary_south_fraction(
            self.latitude(row),
            self.sun_latitude,
            belts,
            &circulation.winds,
        );
        (target as usize, wind.east, south.abs())
    }

    pub(crate) fn downhill_map(&self) -> Vec<Option<usize>> {
        (0..self.len()).map(|i| self.downhill(i)).collect()
    }

    fn vent(&self, i: usize, total: f64, next: &mut [f64]) {
        let (neighbors, count) = self.neighbors(i);
        let share = total / count as f64;
        for &j in &neighbors[..count] {
            next[j] += share * self.area_ratio(i, j);
        }
    }

    pub fn breadth(&self, row: u32) -> f64 {
        row_breadth(row, self.height)
    }

    pub fn area(&self, i: usize) -> f64 {
        self.breadth(i as u32 / self.width)
    }

    pub(crate) fn area_ratio(&self, from: usize, to: usize) -> f64 {
        self.area(from) / self.area(to)
    }

    pub fn total(&self, per_area: &[f64]) -> f64 {
        per_area
            .iter()
            .enumerate()
            .map(|(i, amount)| amount * self.area(i))
            .sum()
    }

    pub(crate) fn row_weights(&self, row: usize) -> [f64; 3] {
        let own = self.breadth(row as u32);
        let toward = |other: usize| (own + self.breadth(other as u32)) / (2.0 * own);
        [
            row.checked_sub(1).map_or(0.0, toward),
            1.0,
            if row + 1 < self.height as usize {
                toward(row + 1)
            } else {
                0.0
            },
        ]
    }
}

pub fn row_breadth(row: u32, height: u32) -> f64 {
    let latitude = (f64::from(row) + 0.5) / f64::from(height) * 2.0 - 1.0;
    libm::cos(latitude * FRAC_PI_2)
}

pub fn deposit(grid: &Grid, weights: &[f64], mean_amount: f64) -> Vec<f64> {
    let weighted: f64 = weights
        .iter()
        .enumerate()
        .map(|(i, w)| w * grid.area(i))
        .sum();
    if weighted <= 0.0 {
        return vec![mean_amount; weights.len()];
    }
    let total = mean_amount * grid.total(&vec![1.0; weights.len()]);
    weights.iter().map(|w| total * w / weighted).collect()
}

fn equilibrate(grid: &Grid, substance: &Substance, pools: &mut Pools, flux: &mut [Flux]) {
    for (i, cell_flux) in flux.iter_mut().enumerate() {
        let saturation = substance.material.saturation(grid.condensation_celsius[i]);
        if substance.sea == Sea::Reservoir && grid.ocean[i] {
            if saturation.is_finite() {
                pools.vapor[i] = match substance.evaporation_limit {
                    Some(limit) => (pools.vapor[i] + limit).min(saturation),
                    None => saturation,
                };
            }
            pools.condensed[i] = 0.0;
            continue;
        }
        let total = pools.vapor[i] + pools.condensed[i];
        let vapor = total.min(saturation);
        let condensed = total - vapor;
        if condensed > pools.condensed[i] {
            cell_flux.condensed += condensed - pools.condensed[i];
        }
        pools.vapor[i] = vapor;
        pools.condensed[i] = condensed;
    }
}

fn soak_and_breathe(grid: &Grid, substance: &Substance, pools: &mut Pools, flux: &mut [Flux]) {
    let capacity = substance.soil_capacity;
    if capacity <= 0.0 {
        return;
    }
    for (i, cell_flux) in flux.iter_mut().enumerate() {
        if grid.ocean[i] {
            continue;
        }
        let condensed = pools.condensed[i];
        if condensed > 0.0
            && substance
                .material
                .condensed_phase(grid.celsius[i], condensed)
                == Phase::Liquid
        {
            let soaked = condensed.min(capacity - pools.soil[i]).max(0.0);
            pools.condensed[i] -= soaked;
            pools.soil[i] += soaked;
            cell_flux.soaked += soaked;
        }
        let deficit = (substance.material.saturation(grid.celsius[i]) - pools.vapor[i]).max(0.0);
        let wetness = pools.soil[i] / capacity;
        let rate = substance
            .evaporation_limit
            .map_or(deficit, |limit| limit * wetness);
        let breathed = pools.soil[i].min(rate).min(deficit);
        pools.soil[i] -= breathed;
        pools.vapor[i] += breathed;
        cell_flux.evaporated += breathed;
    }
}

pub(crate) fn add_row_and_sides(row: &[f64], weight: f64, gathered: &mut [f64]) {
    let width = row.len();
    let at = |x: usize, dx: isize| row[(x as isize + dx).rem_euclid(width as isize) as usize];
    if width >= 3 {
        for (sum, three) in gathered[1..width - 1].iter_mut().zip(row.windows(3)) {
            *sum += weight * (three[0] + three[1] + three[2]);
        }
        for x in [0, width - 1] {
            gathered[x] += weight * (at(x, -1) + row[x] + at(x, 1));
        }
    } else {
        for (x, sum) in gathered.iter_mut().enumerate() {
            *sum += weight * (at(x, -1) + row[x] + at(x, 1));
        }
    }
}

pub(crate) fn gather_mixing(grid: &Grid, field: &[f64], y: usize, gathered: &mut [f64]) -> f64 {
    let width = grid.width as usize;
    let weights = grid.row_weights(y);
    gathered.fill(0.0);
    let mut neighbors = -1.0;
    for (r, weight) in [y.checked_sub(1), Some(y), Some(y + 1)]
        .into_iter()
        .zip(weights)
    {
        let Some(r) = r.filter(|_| weight > 0.0) else {
            continue;
        };
        add_row_and_sides(&field[r * width..(r + 1) * width], weight, gathered);
        neighbors += 3.0 * weight;
    }
    neighbors
}

fn wrap(x: usize, dx: i64, width: usize) -> usize {
    (x as i64 + dx).rem_euclid(width as i64) as usize
}

fn carry_vapor(grid: &Grid, circulation: &Circulation, vapor: &[f64], next: &mut [f64]) -> f64 {
    let width = grid.width as usize;
    let height = grid.height as usize;
    let share = circulation.gas_spread / 8.0;
    let carry = circulation.wind_carry;
    let winds: Vec<(usize, i64, f64)> = (0..grid.height)
        .map(|row| grid.downwind_row(row, circulation))
        .collect();
    let mut moved = 0.0;
    let mut gathered = vec![0.0; width];
    for y in 0..height {
        let neighbors = gather_mixing(grid, vapor, y, &mut gathered);
        let leaving = share * neighbors + carry;
        let row = &vapor[y * width..(y + 1) * width];
        let out_row = &mut next[y * width..(y + 1) * width];
        for (x, out) in out_row.iter_mut().enumerate() {
            let own = row[x];
            let leaves = own * leaving;
            *out = own - leaves + share * (gathered[x] - own);
            moved += leaves;
        }
        let (_, east, south) = winds[y];
        for (x, out) in out_row.iter_mut().enumerate() {
            *out += row[wrap(x, -east, width)] * carry * (1.0 - south);
        }
        for r in [y.checked_sub(1), (y + 1 < height).then_some(y + 1)]
            .into_iter()
            .flatten()
        {
            let (target, _, south) = winds[r];
            if target != y {
                continue;
            }
            let arriving = carry * south * grid.breadth(r as u32) / grid.breadth(y as u32);
            let source = &vapor[r * width..(r + 1) * width];
            for (out, &from) in out_row.iter_mut().zip(source) {
                *out += from * arriving;
            }
        }
    }
    moved
}

fn transport(
    grid: &Grid,
    substance: &Substance,
    circulation: &Circulation,
    downhill: &[Option<usize>],
    pools: &Pools,
    flux: &mut [Flux],
) -> (Pools, f64) {
    let mut next = Pools::empty(grid.len());
    next.soil.copy_from_slice(&pools.soil);
    let vapor: Vec<f64> = pools.vapor.iter().map(|&v| v.max(0.0)).collect();
    let mut moved = carry_vapor(grid, circulation, &vapor, &mut next.vapor);
    for (i, cell_flux) in flux.iter_mut().enumerate() {
        let condensed = pools.condensed[i];
        if condensed <= 0.0 {
            continue;
        }
        match substance
            .material
            .condensed_phase(grid.celsius[i], condensed)
        {
            Phase::Solid | Phase::Gas => next.condensed[i] += condensed,
            Phase::Liquid => match downhill[i] {
                Some(j) => {
                    next.condensed[j] += condensed * grid.area_ratio(i, j);
                    cell_flux.flow += condensed;
                    moved += condensed;
                }
                None if grid.ocean[i] => {
                    grid.vent(i, condensed, &mut next.vapor);
                    cell_flux.vented += condensed;
                    moved += condensed;
                }
                None => next.condensed[i] += condensed,
            },
        }
    }
    (next, moved)
}

pub fn step(
    grid: &Grid,
    substance: &Substance,
    circulation: &Circulation,
    pools: Pools,
    flux: &mut [Flux],
) -> (Pools, f64) {
    step_along(
        grid,
        substance,
        circulation,
        &grid.downhill_map(),
        pools,
        flux,
    )
}

fn step_along(
    grid: &Grid,
    substance: &Substance,
    circulation: &Circulation,
    downhill: &[Option<usize>],
    mut pools: Pools,
    flux: &mut [Flux],
) -> (Pools, f64) {
    equilibrate(grid, substance, &mut pools, flux);
    soak_and_breathe(grid, substance, &mut pools, flux);
    transport(grid, substance, circulation, downhill, &pools, flux)
}

pub struct Weather<'a> {
    pub pattern: &'a [f64],
    pub amplitude_celsius: f64,
}

impl Weather<'_> {
    pub const CALM: Weather<'static> = Weather {
        pattern: &[],
        amplitude_celsius: 0.0,
    };

    pub(crate) fn apply(&self, base: &[f64], width: usize, cells_moved: f64, out: &mut [f64]) {
        if self.pattern.is_empty() {
            out.copy_from_slice(base);
            return;
        }
        let shift = (libm::floor(cells_moved) as usize) % width;
        let rows = out
            .chunks_mut(width)
            .zip(base.chunks(width))
            .zip(self.pattern.chunks(width));
        for ((out_row, base_row), pattern_row) in rows {
            for (x, (out, base)) in out_row.iter_mut().zip(base_row).enumerate() {
                let source = if x >= shift {
                    x - shift
                } else {
                    x + width - shift
                };
                *out = base - pattern_row[source] * self.amplitude_celsius;
            }
        }
    }
}

pub fn settle(
    grid: &Grid,
    substance: &Substance,
    circulation: &Circulation,
    weather: &Weather,
    initial: Pools,
) -> Settled {
    let mut pools = initial;
    let mut flux = vec![Flux::default(); grid.len()];
    let mut moved_per_step = Vec::with_capacity(circulation.age_steps as usize);
    let mut condensation = grid.condensation_celsius.to_vec();
    let downhill = grid.downhill_map();
    for s in 0..circulation.age_steps {
        weather.apply(
            grid.condensation_celsius,
            grid.width as usize,
            f64::from(s) * circulation.wind_carry,
            &mut condensation,
        );
        let weathered = Grid {
            condensation_celsius: &condensation,
            ..*grid
        };
        let (next, moved) = step_along(
            &weathered,
            substance,
            circulation,
            &downhill,
            pools,
            &mut flux,
        );
        pools = next;
        moved_per_step.push(moved);
    }
    equilibrate(
        grid,
        substance,
        &mut pools,
        &mut vec![Flux::default(); grid.len()],
    );
    let steps = f64::from(circulation.age_steps.max(1));
    for f in &mut flux {
        f.flow /= steps;
        f.condensed /= steps;
        f.vented /= steps;
        f.soaked /= steps;
        f.evaporated /= steps;
    }
    Settled {
        pools,
        flux,
        moved_per_step,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::{ipse, water};

    const CIRCULATION: Circulation = Circulation {
        gas_spread: 0.5,
        wind_carry: 0.25,
        age_steps: 60,
        belts: Belts::EARTHLIKE,
        winds: crate::wind::EARTH_WINDS,
    };

    struct World {
        width: u32,
        height: u32,
        elevation: Vec<f64>,
        celsius: Vec<f64>,
        condensation: Vec<f64>,
        ocean: Vec<bool>,
    }

    impl World {
        fn new(width: u32, height: u32, elevation: Vec<f64>, celsius: Vec<f64>) -> Self {
            Self {
                width,
                height,
                ocean: elevation.iter().map(|&e| e < 0.0).collect(),
                condensation: celsius.clone(),
                elevation,
                celsius,
            }
        }

        fn bowl(width: u32, height: u32, celsius: f64) -> Self {
            let n = (width * height) as usize;
            let elevation = (0..n)
                .map(|i| {
                    let x = (i % width as usize) as f64;
                    let y = (i / width as usize) as f64;
                    (x - 3.0).abs() + (y - 2.0).abs()
                })
                .collect();
            Self::new(width, height, elevation, vec![celsius; n])
        }

        fn from_profile(profile: impl Fn(u32) -> f64, celsius: impl Fn(f64) -> f64) -> Self {
            let (width, height) = (16, 12);
            let elevation: Vec<f64> = (0..width * height).map(|i| profile(i % width)).collect();
            let celsius = elevation.iter().map(|&e| celsius(e)).collect();
            Self::new(width, height, elevation, celsius)
        }

        fn ridge() -> Self {
            Self::from_profile(
                |x| match x {
                    0..=4 => -1.0,
                    5..=8 => f64::from(x - 4) / 4.0,
                    _ => 0.1,
                },
                |e| 25.0 - 25.0 * e.max(0.0),
            )
        }

        fn lowland() -> Self {
            Self::from_profile(|x| if x <= 4 { -1.0 } else { 0.1 }, |_| 25.0)
        }

        fn lift_land_air(mut self, cooling: f64) -> Self {
            for (i, condensation) in self.condensation.iter_mut().enumerate() {
                if !self.ocean[i] {
                    *condensation = self.celsius[i] - cooling;
                }
            }
            self
        }

        fn grid(&self) -> Grid<'_> {
            Grid {
                width: self.width,
                height: self.height,
                elevation: &self.elevation,
                celsius: &self.celsius,
                condensation_celsius: &self.condensation,
                ocean: &self.ocean,
                sun_latitude: 0.0,
            }
        }

        fn cells(&self) -> usize {
            self.elevation.len()
        }
    }

    fn inert(material: &Material, sea: Sea) -> Substance<'_> {
        Substance {
            material,
            sea,
            soil_capacity: 0.0,
            evaporation_limit: None,
        }
    }

    fn with_soil(material: &Material, soil_capacity: f64, evaporation_limit: f64) -> Substance<'_> {
        Substance {
            soil_capacity,
            evaporation_limit: Some(evaporation_limit),
            ..inert(material, Sea::Reservoir)
        }
    }

    #[test]
    fn rain_soaks_into_the_ground_before_running_off() {
        let world = World::bowl(8, 5, 20.0);
        let material = water();
        let substance = with_soil(&material, 0.125, 0.0);
        let mut pools = condensed_at(40, 21, 0.5);
        let mut flux = vec![Flux::default(); 40];
        soak_and_breathe(&world.grid(), &substance, &mut pools, &mut flux);
        assert_eq!(pools.soil[21], 0.125);
        assert_eq!(pools.condensed[21], 0.375);
        assert_eq!(flux[21].soaked, 0.125);
    }

    #[test]
    fn soil_water_breathes_back_into_the_air_at_a_limited_rate() {
        let world = World::bowl(8, 5, 20.0);
        let material = water();
        let rate = 1.0 / 256.0;
        let substance = with_soil(&material, 0.125, rate);
        let mut pools = Pools::empty(40);
        pools.soil[21] = 0.125;
        let mut flux = vec![Flux::default(); 40];
        soak_and_breathe(&world.grid(), &substance, &mut pools, &mut flux);
        assert_eq!(pools.vapor[21], rate);
        assert_eq!(pools.soil[21], 0.125 - rate);
        assert_eq!(flux[21].evaporated, rate);
    }

    #[test]
    fn stored_rain_falls_again_further_inland() {
        let world = World::lowland().lift_land_air(6.5);
        let material = water();
        let inland_rain = |substance: &Substance| -> f64 {
            let settled = settle(
                &world.grid(),
                substance,
                &CIRCULATION,
                &Weather::CALM,
                Pools::humid(&world.grid(), substance.material, substance.soil_capacity),
            );
            (0..world.cells())
                .filter(|&i| i % 16 >= 10)
                .map(|i| settled.flux[i].condensed)
                .sum()
        };
        let bare = inland_rain(&with_soil(&material, 0.0, 0.0001));
        let soaked = inland_rain(&with_soil(&material, 0.075, 0.0001));
        assert!(soaked > bare, "bare {bare} soaked {soaked}");
    }

    fn one_step(world: &World, material: &Material, sea: Sea, pools: Pools) -> (Pools, Vec<Flux>) {
        let mut flux = vec![Flux::default(); world.cells()];
        let substance = inert(material, sea);
        let (next, _) = step(&world.grid(), &substance, &CIRCULATION, pools, &mut flux);
        (next, flux)
    }

    fn condensed_at(cells: usize, index: usize, amount: f64) -> Pools {
        let mut pools = Pools::empty(cells);
        pools.condensed[index] = amount;
        pools
    }

    #[test]
    fn a_substance_without_a_sea_is_conserved() {
        let world = World::bowl(8, 5, 20.0);
        let material = ipse();
        let substance = inert(&material, Sea::SinksThrough);
        let grid = world.grid();
        let initial = Pools::with_condensed(deposit(&grid, &world.elevation, 50.0));
        let total = |p: &Pools| grid.total(&p.vapor) + grid.total(&p.condensed);
        let before = total(&initial);
        let settled = settle(
            &world.grid(),
            &substance,
            &CIRCULATION,
            &Weather::CALM,
            initial,
        );
        assert!((total(&settled.pools) - before).abs() <= before * 1e-12);
    }

    #[test]
    fn excess_beyond_saturation_condenses_and_runs_downhill() {
        let world = World::bowl(8, 5, 20.0);
        let (next, flux) = one_step(
            &world,
            &ipse(),
            Sea::SinksThrough,
            condensed_at(40, 21, 60.0),
        );
        assert_eq!(flux[21].condensed, 0.0);
        assert_eq!(next.condensed[20], 20.0);
        assert_eq!(flux[21].flow, 20.0);
    }

    #[test]
    fn liquid_pools_in_land_basins_and_vents_from_ocean_basins() {
        let mut world = World::bowl(8, 5, 20.0);
        let bottom = 2 * 8 + 3;
        let (pooled, _) = one_step(
            &world,
            &ipse(),
            Sea::SinksThrough,
            condensed_at(40, bottom, 60.0),
        );
        assert_eq!(pooled.condensed[bottom], 20.0);

        world.ocean[bottom] = true;
        let (vented, flux) = one_step(
            &world,
            &ipse(),
            Sea::SinksThrough,
            condensed_at(40, bottom, 60.0),
        );
        assert_eq!(vented.condensed[bottom], 0.0);
        assert_eq!(flux[bottom].vented, 20.0);
    }

    fn downwind(grid: &Grid, cell: usize) -> [(usize, f64); 2] {
        let width = grid.width as usize;
        let (row, x) = (cell / width, cell % width);
        let (target, east, south) = grid.downwind_row(row as u32, &CIRCULATION);
        [
            (row * width + wrap(x, east, width), 1.0 - south),
            (target * width + x, south),
        ]
    }

    #[test]
    fn vapor_drifts_downwind() {
        let mut world = World::bowl(8, 6, 20.0);
        world.elevation.fill(0.0);
        let grid = world.grid();
        let source = 2 * 8 + 4;
        let mut pools = Pools::empty(48);
        pools.vapor[source] = 8.0;
        let material = ipse();
        let substance = inert(&material, Sea::SinksThrough);
        let step_with = |circulation: &Circulation| {
            let mut flux = vec![Flux::default(); world.cells()];
            step(&grid, &substance, circulation, pools.clone(), &mut flux).0
        };
        let windy = step_with(&CIRCULATION);
        let still = step_with(&Circulation {
            wind_carry: 0.0,
            ..CIRCULATION
        });
        for (cell, share) in downwind(&grid, source) {
            let drifted = (windy.vapor[cell] - still.vapor[cell]) * grid.area(cell);
            let carried = 8.0 * CIRCULATION.wind_carry * share * grid.area(source);
            assert!(carried > 0.0);
            assert!(
                (drifted - carried).abs() <= carried * 1e-12,
                "{drifted} {carried}"
            );
        }
    }

    #[test]
    fn the_sea_saturates_the_air_and_swallows_rivers() {
        let world = World::ridge();
        let material = water();
        let mut pools = Pools::empty(world.cells());
        pools.condensed[0] = 5.0;
        let mut flux = vec![Flux::default(); world.cells()];
        let substance = inert(&material, Sea::Reservoir);
        equilibrate(&world.grid(), &substance, &mut pools, &mut flux);
        assert_eq!(pools.vapor[0], material.saturation(25.0));
        assert_eq!(pools.condensed[0], 0.0);
    }

    #[test]
    fn rain_falls_on_the_windward_slope() {
        let world = World::ridge();
        let material = water();
        let substance = inert(&material, Sea::Reservoir);
        let settled = settle(
            &world.grid(),
            &substance,
            &CIRCULATION,
            &Weather::CALM,
            Pools::empty(world.cells()),
        );
        let rain = |x: usize| -> f64 { (2..4).map(|y| settled.flux[y * 16 + x].condensed).sum() };
        assert!(
            rain(7) > rain(11) * 2.0,
            "windward {} leeward {}",
            rain(7),
            rain(11)
        );
    }

    #[test]
    fn lifted_air_rains_on_warm_lowlands() {
        let land_rain = |world: World| -> f64 {
            let material = water();
            let substance = inert(&material, Sea::Reservoir);
            let settled = settle(
                &world.grid(),
                &substance,
                &CIRCULATION,
                &Weather::CALM,
                Pools::empty(world.cells()),
            );
            (0..world.cells())
                .filter(|&i| !world.ocean[i])
                .map(|i| settled.flux[i].condensed)
                .sum()
        };
        let still = land_rain(World::lowland());
        let lifted = land_rain(World::lowland().lift_land_air(6.5));
        assert!(lifted > still, "still {still} lifted {lifted}");
    }

    #[test]
    fn passing_weather_brings_rain_to_more_land() {
        let world = World::lowland();
        let material = water();
        let substance = inert(&material, Sea::Reservoir);
        let pattern: Vec<f64> = (0..world.cells())
            .map(|i| if (i % 16 / 4) % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let rained_on = |weather: &Weather| -> usize {
            let settled = settle(
                &world.grid(),
                &substance,
                &CIRCULATION,
                weather,
                Pools::empty(world.cells()),
            );
            (0..world.cells())
                .filter(|&i| !world.ocean[i] && settled.flux[i].condensed > 0.0)
                .count()
        };
        let calm = rained_on(&Weather::CALM);
        let stormy = rained_on(&Weather {
            pattern: &pattern,
            amplitude_celsius: 6.0,
        });
        assert!(stormy > calm, "calm {calm} stormy {stormy}");
    }

    #[test]
    fn wind_belts_follow_the_sun() {
        let world = World::bowl(8, 12, 20.0);
        let mut grid = world.grid();
        let row = 3;
        let cell = row * 8 + 4;
        assert!((grid.latitude(row as u32) * 90.0).abs() > 30.0);
        let westerly = downwind(&grid, cell)[0].0;
        grid.sun_latitude = grid.latitude(row as u32);
        let trade = downwind(&grid, cell)[0].0;
        assert_eq!(westerly % 8, 5);
        assert_eq!(trade % 8, 3);
    }

    #[test]
    fn neighbors_wrap_east_west() {
        let world = World::bowl(8, 5, 20.0);
        let (neighbors, count) = world.grid().neighbors(2 * 8);
        assert!(neighbors[..count].contains(&(2 * 8 + 7)));
    }
}
