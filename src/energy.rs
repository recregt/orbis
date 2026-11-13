use alloc::vec;
use alloc::vec::Vec;
use core::f64::consts::{FRAC_PI_2, PI, TAU};

use crate::atmosphere::{LAYERS, RADIATION_SLABS, SLABS, Slab, TROPOSPHERE};
use crate::cycle::Grid;
use crate::material::Materials;
use crate::ocean::Currents;
use crate::parallel::ForEach;
use crate::radiation::{self, Longwave, glow};
use crate::scale::{
    GRAVITY_METRES_PER_SECOND_SQUARED, KELVIN_AT_ZERO_CELSIUS, SEA_LEVEL_AIR_KG_PER_CUBIC_METRE,
    SECONDS_PER_DAY, STEFAN_BOLTZMANN_W_PER_M2_PER_K4, SURFACE_MOISTURE_TRANSFER_COEFFICIENT,
    mixing_rate,
};
use crate::sky::Sky;
use crate::wind::{ground_friction_w_per_m2, surface_wind_metres_per_second};

pub const OUTGOING_W_PER_M2_AT_ZERO_CELSIUS: f64 = 203.3;
pub const OUTGOING_W_PER_M2_PER_KELVIN: f64 = 2.09;
pub const UNFROZEN_COALBEDO_MEAN: f64 = 0.697;
pub const UNFROZEN_COALBEDO_TOWARD_POLES: f64 = -0.0779;
pub const ALBEDO_FROZEN: f64 = 0.62;
pub const ICE_EDGE_ANNUAL_MEAN_CELSIUS: f64 = -10.0;
pub const OCEAN_MIXED_LAYER_METRES: f64 = 50.0;
pub const SETTLED_KELVIN: f64 = 0.01;
pub const YEARS_AT_MOST: u32 = 200;
pub const OVERTURNING_LAYERS: f64 = 2.0;

const PASCALS_PER_HECTOPASCAL: f64 = 100.0;

pub fn aether_warming_kelvin(aether_w_per_m2: f64, absorptance: f64, celsius: f64) -> f64 {
    let kelvin = celsius + KELVIN_AT_ZERO_CELSIUS;
    absorptance * aether_w_per_m2
        / (4.0 * STEFAN_BOLTZMANN_W_PER_M2_PER_K4 * kelvin * kelvin * kelvin)
}

fn icebound_at(annual_mean_celsius: f64) -> bool {
    annual_mean_celsius < ICE_EDGE_ANNUAL_MEAN_CELSIUS
}

fn absorbed(insolation: f64, icebound: bool, unfrozen_albedo: f64) -> f64 {
    let albedo = if icebound {
        ALBEDO_FROZEN
    } else {
        unfrozen_albedo
    };
    insolation * (1.0 - albedo)
}

pub(crate) fn annual_insolation(grid: &Grid, sky: &Sky) -> Vec<f64> {
    let days = libm::ceil(sky.year_seconds() / SECONDS_PER_DAY) as usize;
    let day_seconds = sky.year_seconds() / days as f64;
    let mut annual = vec![0.0; grid.height as usize];
    for day in 0..days {
        let seconds = (day as f64 + 0.5) * day_seconds;
        for (row, total) in (0..grid.height).zip(&mut annual) {
            *total += sky.daily_insolation(seconds, -grid.latitude(row) * FRAC_PI_2);
        }
    }
    annual.iter().map(|total| total / days as f64).collect()
}

fn return_temperatures(overturning: &[f64], width: usize, surface: &[f64], deep: &mut [f64]) {
    let cells = surface.len();
    let settle = |i: usize, deep: &mut [f64]| {
        let north_face = overturning[i];
        let south_face = overturning.get(i + width).copied().unwrap_or(0.0);
        let from_north = if i >= width { north_face.max(0.0) } else { 0.0 };
        let from_south = (-south_face).max(0.0);
        let sinking = (south_face - north_face).max(0.0);
        let arriving = from_north + from_south + sinking;
        deep[i] = if arriving > 0.0 {
            let north = if from_north > 0.0 {
                deep[i - width]
            } else {
                0.0
            };
            let south = if from_south > 0.0 {
                deep[i + width]
            } else {
                0.0
            };
            (from_north * north + from_south * south + sinking * surface[i]) / arriving
        } else {
            surface[i]
        };
    };
    for i in 0..cells {
        settle(i, deep);
    }
    for i in (0..cells).rev() {
        settle(i, deep);
    }
}

const HARMONICS: usize = 4;

pub struct EnergyBalance<'a> {
    pub grid: &'a Grid<'a>,
    pub cell_metres: f64,
    pub sky: &'a Sky<'a>,
    pub materials: &'a Materials,
    pub atmosphere_diffusivity_m2_per_s: f64,
}

pub struct Exchanged<'a> {
    pub currents: &'a Currents,
    pub layer_humidity: &'a [[f64; LAYERS]],
    pub wetness: &'a [f64],
    pub wind: Wind,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Wind {
    Given(f64),
    FromHeatEngine,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Surface {
    pub to_space_w_per_m2: f64,
    pub down_to_ground_w_per_m2: f64,
    pub air_radiative_cooling_w_per_m2: f64,
    pub sensible_w_per_m2: f64,
    pub latent_w_per_m2: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SeasonalCycle {
    pub mean: Vec<f64>,
    pub warmest: Vec<f64>,
    pub harmonics: Vec<[f64; HARMONICS]>,
    pub ground_mean: Vec<f64>,
    pub ground_harmonics: Vec<[f64; HARMONICS]>,
    pub year_seconds: f64,
    pub years: u32,
    pub settled: bool,
    pub northward: Vec<Transport>,
    pub flux_equator: f64,
    pub flux_equator_harmonics: [f64; HARMONICS],
    pub atmosphere_work_w_per_m2: f64,
    pub wind_metres_per_second: f64,
    pub surface: Surface,
    pub free_air_cooling_w_per_m2: Vec<f64>,
    pub end_of_year: Vec<f64>,
    pub ground_end_of_year: Vec<f64>,
}

fn waves(angle: f64) -> [f64; HARMONICS] {
    [
        libm::cos(angle),
        libm::sin(angle),
        libm::cos(2.0 * angle),
        libm::sin(2.0 * angle),
    ]
}

fn flux_equator(atmosphere_gain_by_row: &[f64]) -> f64 {
    let height = atmosphere_gain_by_row.len();
    let spacing = 2.0 / height as f64;
    let mut northward = 0.0;
    let mut previous: Option<(f64, f64)> = None;
    let mut closest: Option<f64> = None;
    for (row, gain) in atmosphere_gain_by_row[..height.saturating_sub(1)]
        .iter()
        .enumerate()
    {
        northward += gain;
        let latitude = (row as f64 + 1.0) * spacing - 1.0;
        if let Some((above, above_latitude)) = previous
            && above > 0.0
            && northward <= 0.0
        {
            let crossing = above_latitude + spacing * above / (above - northward);
            if closest.is_none_or(|c: f64| crossing.abs() < c.abs()) {
                closest = Some(crossing);
            }
        }
        previous = Some((northward, latitude));
    }
    closest.unwrap_or(0.0)
}

pub fn north_south_work_watts(northward_watts: &[f64], row_kelvin: &[f64]) -> f64 {
    northward_watts
        .iter()
        .zip(row_kelvin.windows(2))
        .map(|(crossing, pair)| {
            let (north, south) = (pair[0], pair[1]);
            crossing * (south - north) / ((north + south) / 2.0)
        })
        .sum()
}

fn row_kelvin(celsius: &[f64], width: usize) -> Vec<f64> {
    celsius
        .chunks(width)
        .map(|row| row.iter().sum::<f64>() / width as f64 + KELVIN_AT_ZERO_CELSIUS)
        .collect()
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Transport {
    pub atmosphere_watts: f64,
    pub ocean_watts: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Radiative {
    pub contrast: f64,
    pub equator_to_pole_kelvin: f64,
    pub mean_kelvin: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Budget {
    ground: f64,
    air: f64,
    ground_damping: f64,
    air_damping: f64,
    atmosphere: f64,
    ocean: f64,
    longwave: Longwave,
    sunlight_in_air: f64,
    sensible: f64,
    latent: f64,
    parcels: [f64; SLABS],
}

impl SeasonalCycle {
    pub fn phase(&self, seconds: f64) -> [f64; HARMONICS] {
        waves(TAU * seconds / self.year_seconds)
    }

    pub fn flux_equator_at(&self, phase: &[f64; HARMONICS]) -> f64 {
        self.flux_equator
            + self
                .flux_equator_harmonics
                .iter()
                .zip(phase)
                .map(|(amplitude, wave)| amplitude * wave)
                .sum::<f64>()
    }

    pub fn anomaly(&self, i: usize, phase: &[f64; HARMONICS]) -> f64 {
        self.harmonics[i]
            .iter()
            .zip(phase)
            .map(|(amplitude, wave)| amplitude * wave)
            .sum()
    }

    pub fn ground_anomaly(&self, i: usize, phase: &[f64; HARMONICS]) -> f64 {
        self.ground_harmonics[i]
            .iter()
            .zip(phase)
            .map(|(amplitude, wave)| amplitude * wave)
            .sum()
    }
}

impl EnergyBalance<'_> {
    fn air_column_capacity(&self) -> f64 {
        self.materials.air.specific_heat * TROPOSPHERE[0].base_hpa * PASCALS_PER_HECTOPASCAL
            / GRAVITY_METRES_PER_SECOND_SQUARED
    }

    fn ground_capacities(&self, year_seconds: f64) -> Vec<f64> {
        let (water, rock) = (&self.materials.water, &self.materials.rock);
        let sea = water.specific_heat * water.density * OCEAN_MIXED_LAYER_METRES;
        let annual_depth = libm::sqrt(rock.thermal_diffusivity() * year_seconds / PI);
        let land = rock.specific_heat * rock.density * annual_depth;
        self.grid
            .ocean
            .iter()
            .map(|&ocean| if ocean { sea } else { land })
            .collect()
    }

    fn latitude_north(&self, row: u32) -> f64 {
        -self.grid.latitude(row) * FRAC_PI_2
    }

    pub fn radiative(&self) -> Radiative {
        let annual = annual_insolation(self.grid, self.sky);
        let (mut mean, mut shape, mut shape_norm, mut area) = (0.0, 0.0, 0.0, 0.0);
        for row in 0..self.grid.height {
            let r = row as usize;
            let equilibrium = (absorbed(annual[r], false, self.unfrozen_albedo(row))
                - OUTGOING_W_PER_M2_AT_ZERO_CELSIUS)
                / OUTGOING_W_PER_M2_PER_KELVIN
                + KELVIN_AT_ZERO_CELSIUS;
            let x = libm::sin(self.latitude_north(row));
            let second_legendre = (3.0 * x * x - 1.0) / 2.0;
            let weight = self.grid.breadth(row);
            mean += equilibrium * weight;
            area += weight;
            shape += equilibrium * second_legendre * weight;
            shape_norm += second_legendre * second_legendre * weight;
        }
        let equator_to_pole_kelvin = -1.5 * (shape / shape_norm);
        Radiative {
            contrast: equator_to_pole_kelvin / (mean / area),
            equator_to_pole_kelvin,
            mean_kelvin: mean / area,
        }
    }

    fn unfrozen_albedo(&self, row: u32) -> f64 {
        let x = libm::sin(self.latitude_north(row));
        let second_legendre = (3.0 * x * x - 1.0) / 2.0;
        1.0 - (UNFROZEN_COALBEDO_MEAN + UNFROZEN_COALBEDO_TOWARD_POLES * second_legendre)
    }

    fn specific_humidity(&self, celsius: f64, relative_humidity: f64) -> f64 {
        relative_humidity * self.materials.water.saturation(celsius) / self.materials.air.density
    }

    fn boundary_humidity(&self, air_celsius: f64, relative_humidity: f64) -> f64 {
        self.specific_humidity(air_celsius, relative_humidity)
    }

    fn moist_enthalpy(&self, air_celsius: f64, relative_humidity: f64) -> f64 {
        let (air, water) = (&self.materials.air, &self.materials.water);
        air_celsius
            + water.latent_heat / air.specific_heat
                * self.boundary_humidity(air_celsius, relative_humidity)
    }

    fn ocean_faces(&self, ocean_transport: &[f64]) -> Vec<[(usize, f64); 4]> {
        let grid = self.grid;
        let (width, height) = (grid.width as usize, grid.height as usize);
        let wet = |x: usize, y: Option<usize>| {
            y.is_some_and(|y| y < height && grid.ocean[y * width + x % width])
        };
        let corner = |cx: usize, cy: usize| -> f64 {
            let around = [
                (cx, cy.checked_sub(1)),
                (cx + 1, cy.checked_sub(1)),
                (cx, Some(cy)),
                (cx + 1, Some(cy)),
            ];
            if around.iter().all(|&(x, y)| wet(x, y)) {
                around
                    .iter()
                    .map(|&(x, y)| ocean_transport[y.unwrap_or(0) * width + x % width])
                    .sum::<f64>()
                    / 4.0
            } else {
                0.0
            }
        };
        (0..grid.len())
            .map(|i| {
                let (x, y) = (i % width, i / width);
                let west_x = x + width - 1;
                let (north_west, north_east) = (corner(west_x, y), corner(x, y));
                let (south_west, south_east) = (corner(west_x, y + 1), corner(x, y + 1));
                let east = (y * width + (x + 1) % width, north_east - south_east);
                let west = (y * width + west_x % width, south_west - north_west);
                let north = (
                    y.checked_sub(1).map_or(i, |n| n * width + x),
                    north_west - north_east,
                );
                let south = (
                    if y + 1 < height {
                        (y + 1) * width + x
                    } else {
                        i
                    },
                    south_east - south_west,
                );
                [east, west, north, south]
            })
            .collect()
    }

    pub fn run(
        &self,
        exchanged: &Exchanged,
        from: Option<&SeasonalCycle>,
        parallel: &dyn ForEach<Budget>,
    ) -> SeasonalCycle {
        let grid = self.grid;
        let cells = grid.len();
        let width = grid.width as usize;
        let year_seconds = self.sky.year_seconds();
        let ground_capacity = self.ground_capacities(year_seconds);
        let air = self.air_column_capacity();
        let (air_material, water) = (&self.materials.air, &self.materials.water);
        let water_heat_per_volume = water.specific_heat * water.density;
        let mixing_rate = mixing_rate(self.atmosphere_diffusivity_m2_per_s, self.cell_metres);
        let transfer_at = |wind_metres_per_second: f64| {
            SEA_LEVEL_AIR_KG_PER_CUBIC_METRE
                * SURFACE_MOISTURE_TRANSFER_COEFFICIENT
                * wind_metres_per_second
        };
        let layer_humidity = exchanged.layer_humidity;
        let breadth: Vec<f64> = (0..grid.height).map(|row| grid.breadth(row)).collect();
        let radius = self.cell_metres * f64::from(grid.width) / TAU;
        let d_phi = PI / f64::from(grid.height);
        let area: Vec<f64> = (0..cells)
            .map(|i| radius * radius * (TAU / width as f64) * d_phi * breadth[i / width])
            .collect();
        let total_area: f64 = area.iter().sum();
        let links: Vec<([(usize, f64); 8], usize)> = (0..cells)
            .map(|i| {
                let (neighbors, count) = grid.neighbors(i);
                let own = breadth[i / width];
                let mut weighted = [(0, 0.0); 8];
                for (slot, &j) in weighted.iter_mut().zip(&neighbors[..count]) {
                    *slot = (j, (own + breadth[j / width]) / (2.0 * own));
                }
                (weighted, count)
            })
            .collect();
        let currents = exchanged.currents;
        let faces = self.ocean_faces(&currents.gyre);
        let overturning = &currents.overturning;
        let overturned_across = |i: usize| overturning.get(i).copied().unwrap_or(0.0);
        let leaving_by_current: Vec<f64> = faces
            .iter()
            .enumerate()
            .map(|(i, faces)| {
                let out: f64 = faces.iter().map(|&(_, f)| (-f).max(0.0)).sum();
                let layers = OVERTURNING_LAYERS
                    * (overturned_across(i).abs() + overturned_across(i + width).abs());
                water_heat_per_volume * (out + layers) / area[i]
            })
            .collect();
        let lifted_north = |ground: &[f64], deep: &[f64], i: usize| -> f64 {
            let flux = overturned_across(i);
            if flux == 0.0 {
                return 0.0;
            }
            let (from, to) = if flux > 0.0 {
                (i, i - width)
            } else {
                (i - width, i)
            };
            water_heat_per_volume * flux * (ground[from] - deep[to])
        };

        let annual_sunlight = annual_insolation(self.grid, self.sky);
        let unfrozen: Vec<f64> = (0..grid.height)
            .map(|row| self.unfrozen_albedo(row))
            .collect();
        let mut air_celsius: Vec<f64> = (0..cells)
            .map(|i| {
                let row = i / width;
                let emitting = absorbed(annual_sunlight[row], false, unfrozen[row])
                    / STEFAN_BOLTZMANN_W_PER_M2_PER_K4;
                libm::sqrt(libm::sqrt(emitting)) - KELVIN_AT_ZERO_CELSIUS
            })
            .collect();
        let mut ground = air_celsius.clone();
        if let Some(earlier) = from {
            air_celsius.clone_from(&earlier.end_of_year);
            ground.clone_from(&earlier.ground_end_of_year);
        }
        let glow_per_kelvin_at = |celsius: f64| {
            let kelvin = celsius + KELVIN_AT_ZERO_CELSIUS;
            4.0 * STEFAN_BOLTZMANN_W_PER_M2_PER_K4 * kelvin * kelvin * kelvin
        };
        let humidity_per_kelvin =
            |celsius: f64| water.saturation_per_kelvin(celsius) / air_material.density;
        let all_air: f64 = RADIATION_SLABS.iter().map(Slab::air_kg_per_m2).sum();
        let mut layer_air = [0.0; LAYERS];
        for slab in RADIATION_SLABS {
            layer_air[slab.layer] += slab.air_kg_per_m2();
        }
        let free_air_share: f64 = layer_air[1..].iter().sum::<f64>() / all_air;
        let mixing_weight: Vec<f64> = links
            .iter()
            .map(|(weighted, count)| weighted[..*count].iter().map(|&(_, w)| w).sum())
            .collect();
        let enthalpies = |air_celsius: &[f64]| -> Vec<f64> {
            air_celsius
                .iter()
                .zip(layer_humidity)
                .map(|(&celsius, layers)| self.moist_enthalpy(celsius, layers[0]))
                .collect()
        };
        let budget = |ground: &[f64],
                      air_celsius: &[f64],
                      enthalpy: &[f64],
                      deep: &[f64],
                      sunlight: &[f64],
                      icebound: &[bool],
                      transfer: f64,
                      previous: Option<&[Budget]>,
                      i: usize|
         -> Budget {
            let row = i / width;
            let (weighted, count) = &links[i];
            let mixed: f64 = weighted[..*count]
                .iter()
                .map(|&(j, w)| w * (enthalpy[j] - enthalpy[i]))
                .sum();
            let carried: f64 = faces[i]
                .iter()
                .map(|&(j, flux)| {
                    let source = if flux > 0.0 { ground[j] } else { ground[i] };
                    water_heat_per_volume * flux * source
                })
                .sum::<f64>();
            let overturned = lifted_north(ground, deep, i + width) - lifted_north(ground, deep, i);
            let column = radiation::column_from(
                previous.map(|earlier| &earlier[i].parcels),
                ground[i],
                air_celsius[i],
                &layer_humidity[i],
                air_material,
                water,
            );
            let longwave = column.longwave;
            let taken = absorbed(sunlight[row], icebound[i], unfrozen[row]);
            let sunlight_in_air = (air_material.sunlight_absorptance * sunlight[row]).min(taken);
            let sensible = transfer * air_material.specific_heat * (ground[i] - air_celsius[i]);
            let latent = transfer
                * water.latent_heat
                * exchanged.wetness[i]
                * (self.specific_humidity(ground[i], 1.0)
                    - self.boundary_humidity(air_celsius[i], layer_humidity[i][0]))
                .max(0.0);
            let ground_glow = glow(ground[i] + KELVIN_AT_ZERO_CELSIUS);
            let atmosphere = air * mixing_rate / 8.0 * mixed;
            let ocean = (carried + overturned) / area[i];
            let relative = layer_humidity[i][0];
            let boundary = air_celsius[i];
            let sensible_per_kelvin = transfer * air_material.specific_heat;
            let latent_per_kelvin = transfer * water.latent_heat * exchanged.wetness[i];
            let enthalpy_per_kelvin = 1.0
                + water.latent_heat / air_material.specific_heat
                    * relative
                    * humidity_per_kelvin(boundary);
            Budget {
                ground_damping: glow_per_kelvin_at(ground[i])
                    + sensible_per_kelvin
                    + latent_per_kelvin * humidity_per_kelvin(ground[i])
                    + leaving_by_current[i],
                air_damping: 2.0 * glow_per_kelvin_at(air_celsius[i])
                    + sensible_per_kelvin
                    + latent_per_kelvin * relative * humidity_per_kelvin(boundary)
                    + air * mixing_rate / 8.0 * mixing_weight[i] * enthalpy_per_kelvin,
                ground: taken - sunlight_in_air + longwave.down_to_ground_w_per_m2
                    - ground_glow
                    - sensible
                    - latent
                    + ocean,
                air: sunlight_in_air + ground_glow
                    - longwave.down_to_ground_w_per_m2
                    - longwave.to_space_w_per_m2
                    + sensible
                    + latent
                    + atmosphere,
                atmosphere,
                ocean,
                longwave,
                sunlight_in_air,
                sensible,
                latent,
                parcels: column.parcels,
            }
        };

        let fill = |out: &mut [Budget],
                    ground: &[f64],
                    air_celsius: &[f64],
                    deep: &[f64],
                    sunlight: &[f64],
                    icebound: &[bool],
                    transfer: f64,
                    previous: Option<&[Budget]>| {
            let enthalpy = enthalpies(air_celsius);
            parallel.each(out, &|i, slot| {
                *slot = budget(
                    ground,
                    air_celsius,
                    &enthalpy,
                    deep,
                    sunlight,
                    icebound,
                    transfer,
                    previous,
                    i,
                );
            });
        };
        let wind_for = |work_w_per_m2: f64| match exchanged.wind {
            Wind::Given(speed) => speed,
            Wind::FromHeatEngine => surface_wind_metres_per_second(ground_friction_w_per_m2(
                work_w_per_m2,
                air_material,
            )),
        };
        let work_now = |ground: &[f64],
                        air_celsius: &[f64],
                        deep: &[f64],
                        icebound: &[bool],
                        transfer: f64| {
            let mut out = vec![Budget::default(); cells];
            fill(
                &mut out,
                ground,
                air_celsius,
                deep,
                &annual_sunlight,
                icebound,
                transfer,
                None,
            );
            let mut crossing = 0.0;
            let northward: Vec<f64> = out
                .chunks(width)
                .zip(area.chunks(width))
                .map(|(row, row_area)| {
                    crossing += row
                        .iter()
                        .zip(row_area)
                        .map(|(budget, a)| budget.atmosphere * a)
                        .sum::<f64>();
                    crossing
                })
                .collect();
            north_south_work_watts(&northward, &row_kelvin(air_celsius, width)) / total_area
        };
        let mut deep = ground.clone();
        let mut icebound: Vec<bool> = air_celsius.iter().map(|&t| icebound_at(t)).collect();
        let mut wind = match (exchanged.wind, from) {
            (Wind::Given(speed), _) => speed,
            (Wind::FromHeatEngine, Some(earlier)) => wind_for(earlier.atmosphere_work_w_per_m2),
            (Wind::FromHeatEngine, None) => {
                wind_for(work_now(&ground, &air_celsius, &deep, &icebound, 0.0))
            }
        };
        let steps = libm::ceil(year_seconds / SECONDS_PER_DAY) as usize;
        let step_seconds = year_seconds / steps as f64;
        let sunlight: Vec<Vec<f64>> = (0..steps)
            .map(|s| {
                let seconds = (s as f64 + 0.5) * step_seconds;
                (0..grid.height)
                    .map(|row| self.sky.daily_insolation(seconds, self.latitude_north(row)))
                    .collect()
            })
            .collect();

        let (mut next_ground, mut next_air) = (vec![0.0; cells], vec![0.0; cells]);
        let mut heats = vec![Budget::default(); cells];
        let mut last = vec![Budget::default(); cells];
        let mut have_last = false;
        let spin_up_years = if let Some(earlier) = from {
            for (frozen, &t) in icebound.iter_mut().zip(&earlier.mean) {
                *frozen = icebound_at(t);
            }
            0
        } else {
            YEARS_AT_MOST
        };
        for _ in 0..spin_up_years {
            let (ground_before, air_before) = (ground.clone(), air_celsius.clone());
            let transfer = transfer_at(wind);
            {
                return_temperatures(overturning, width, &ground, &mut deep);
                fill(
                    &mut heats,
                    &ground,
                    &air_celsius,
                    &deep,
                    &annual_sunlight,
                    &icebound,
                    transfer,
                    have_last.then_some(last.as_slice()),
                );
                for (i, heat) in heats.iter().enumerate() {
                    next_ground[i] = ground[i]
                        + heat.ground * year_seconds
                            / (ground_capacity[i] + year_seconds * heat.ground_damping);
                    next_air[i] = air_celsius[i]
                        + heat.air * year_seconds / (air + year_seconds * heat.air_damping);
                }
                core::mem::swap(&mut ground, &mut next_ground);
                core::mem::swap(&mut air_celsius, &mut next_air);
                core::mem::swap(&mut heats, &mut last);
                have_last = true;
            }
            for (frozen, &t) in icebound.iter_mut().zip(&air_celsius) {
                *frozen = icebound_at(t);
            }
            wind = wind_for(work_now(
                &ground,
                &air_celsius,
                &deep,
                &icebound,
                transfer_at(wind),
            ));
            let unchanged = |before: &[f64], after: &[f64]| {
                before
                    .iter()
                    .zip(after)
                    .all(|(a, b)| (a - b).abs() < SETTLED_KELVIN)
            };
            if unchanged(&ground_before, &ground) && unchanged(&air_before, &air_celsius) {
                break;
            }
        }
        let mut previous_means: Option<(Vec<f64>, Vec<f64>)> = None;
        let mut years = 0;
        loop {
            years += 1;
            let mut mean = vec![0.0; cells];
            let mut ground_mean = vec![0.0; cells];
            let mut warmest = vec![f64::NEG_INFINITY; cells];
            let mut harmonics = vec![[0.0; HARMONICS]; cells];
            let mut ground_harmonics = vec![[0.0; HARMONICS]; cells];
            let mut free_air_cooling = vec![0.0; cells];
            let mut surface = Surface::default();
            let mut gained = vec![Transport::default(); grid.height as usize];
            let mut equator = (0.0, [0.0; HARMONICS]);
            let mut gained_now = vec![0.0; grid.height as usize];
            let transfer = transfer_at(wind);
            for (s, rows) in sunlight.iter().enumerate() {
                gained_now.fill(0.0);
                return_temperatures(overturning, width, &ground, &mut deep);
                fill(
                    &mut heats,
                    &ground,
                    &air_celsius,
                    &deep,
                    rows,
                    &icebound,
                    transfer,
                    have_last.then_some(last.as_slice()),
                );
                for (i, heat) in heats.iter().enumerate() {
                    let row = &mut gained[i / width];
                    row.atmosphere_watts += heat.atmosphere * area[i];
                    row.ocean_watts += heat.ocean * area[i];
                    gained_now[i / width] += heat.atmosphere * area[i];
                    let cooling = heat.longwave.air_cooling_w_per_m2() - heat.sunlight_in_air;
                    free_air_cooling[i] += -heat.longwave.layer_heating_w_per_m2()[1..]
                        .iter()
                        .sum::<f64>()
                        - heat.sunlight_in_air * free_air_share;
                    let share = area[i] / total_area;
                    surface.to_space_w_per_m2 += heat.longwave.to_space_w_per_m2 * share;
                    surface.down_to_ground_w_per_m2 +=
                        heat.longwave.down_to_ground_w_per_m2 * share;
                    surface.air_radiative_cooling_w_per_m2 += cooling * share;
                    surface.sensible_w_per_m2 += heat.sensible * share;
                    surface.latent_w_per_m2 += heat.latent * share;
                    next_ground[i] = ground[i]
                        + heat.ground * step_seconds
                            / (ground_capacity[i] + step_seconds * heat.ground_damping);
                    next_air[i] = air_celsius[i]
                        + heat.air * step_seconds / (air + step_seconds * heat.air_damping);
                }
                core::mem::swap(&mut ground, &mut next_ground);
                core::mem::swap(&mut air_celsius, &mut next_air);
                core::mem::swap(&mut heats, &mut last);
                have_last = true;
                let crossing = flux_equator(&gained_now);
                equator.0 += crossing / steps as f64;
                let middle = waves(TAU * (s as f64 + 0.5) / steps as f64);
                for (h, wave) in equator.1.iter_mut().zip(&middle) {
                    *h += crossing * wave * 2.0 / steps as f64;
                }
                let after = waves(TAU * (s as f64 + 1.0) / steps as f64);
                for i in 0..cells {
                    mean[i] += air_celsius[i];
                    ground_mean[i] += ground[i];
                    warmest[i] = warmest[i].max(air_celsius[i]);
                    for (h, wave) in harmonics[i].iter_mut().zip(&after) {
                        *h += air_celsius[i] * wave;
                    }
                    for (h, wave) in ground_harmonics[i].iter_mut().zip(&after) {
                        *h += ground[i] * wave;
                    }
                }
            }
            let per_step = 1.0 / steps as f64;
            for i in 0..cells {
                mean[i] *= per_step;
                ground_mean[i] *= per_step;
                free_air_cooling[i] *= per_step;
                for h in harmonics[i].iter_mut().chain(&mut ground_harmonics[i]) {
                    *h *= 2.0 * per_step;
                }
                icebound[i] = icebound_at(mean[i]);
            }
            for flux in [
                &mut surface.to_space_w_per_m2,
                &mut surface.down_to_ground_w_per_m2,
                &mut surface.air_radiative_cooling_w_per_m2,
                &mut surface.sensible_w_per_m2,
                &mut surface.latent_w_per_m2,
            ] {
                *flux *= per_step;
            }
            let settled = previous_means
                .as_ref()
                .is_some_and(|(air_before, ground_before)| {
                    air_before
                        .iter()
                        .zip(&mean)
                        .chain(ground_before.iter().zip(&ground_mean))
                        .all(|(a, b)| (a - b).abs() < SETTLED_KELVIN)
                });
            let mut cap = Transport::default();
            let mut northward = Vec::with_capacity(gained.len());
            for row in &gained {
                cap.atmosphere_watts += row.atmosphere_watts * per_step;
                cap.ocean_watts += row.ocean_watts * per_step;
                northward.push(cap);
            }
            let crossings: Vec<f64> = northward.iter().map(|t| t.atmosphere_watts).collect();
            let atmosphere_work_w_per_m2 =
                north_south_work_watts(&crossings, &row_kelvin(&mean, width)) / total_area;
            if settled || years >= YEARS_AT_MOST {
                return SeasonalCycle {
                    wind_metres_per_second: wind,
                    mean,
                    warmest,
                    harmonics,
                    ground_mean,
                    ground_harmonics,
                    year_seconds,
                    years,
                    settled,
                    northward,
                    flux_equator: equator.0,
                    flux_equator_harmonics: equator.1,
                    atmosphere_work_w_per_m2,
                    surface,
                    free_air_cooling_w_per_m2: free_air_cooling,
                    end_of_year: air_celsius,
                    ground_end_of_year: ground,
                };
            }
            wind = wind_for(atmosphere_work_w_per_m2);
            previous_means = Some((mean, ground_mean));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::test_materials;
    use crate::ocean::{Ocean, equatorial_width_metres};
    use crate::sky::{Companion, EARTH_ORBIT, Heavens, sun_alone};
    use crate::wind::{Belts, EARTH_WINDS};

    struct Planet {
        width: u32,
        height: u32,
        zeros: Vec<f64>,
        ocean: Vec<bool>,
        heavens: Heavens,
        layer_humidity: [f64; LAYERS],
    }

    impl Planet {
        fn new(width: u32, height: u32, land: impl Fn(u32, u32) -> bool) -> Self {
            let cells = (width * height) as usize;
            Self {
                width,
                height,
                zeros: vec![0.0; cells],
                ocean: (0..cells as u32)
                    .map(|i| !land(i % width, i / width))
                    .collect(),
                heavens: sun_alone(),
                layer_humidity: EARTH_LAYER_HUMIDITY,
            }
        }

        fn cycle(&self) -> SeasonalCycle {
            self.run(Flow::Full)
        }

        fn run(&self, flow: Flow) -> SeasonalCycle {
            let grid = Grid {
                width: self.width,
                height: self.height,
                elevation: &self.zeros,
                celsius: &self.zeros,
                condensation_celsius: &self.zeros,
                ocean: &self.ocean,
                sun_latitude: 0.0,
            };
            let sky = self.heavens.sky(EARTH_ORBIT);
            let materials = test_materials();
            let balance = EnergyBalance {
                grid: &grid,
                cell_metres: TAU * EARTH_RADIUS_METRES / f64::from(self.width),
                sky: &sky,
                materials: &materials,
                atmosphere_diffusivity_m2_per_s: EARTH_EDDY_DIFFUSIVITY_M2_PER_S,
            };
            let ocean = Ocean {
                grid: &grid,
                belts: Belts::EARTHLIKE,
                radius_metres: EARTH_RADIUS_METRES,
                rotation_radians_per_second: EARTH_ROTATION_RADIANS_PER_SECOND,
                water_kg_per_cubic_metre: materials.water.density,
                equatorial_width_metres: equatorial_width_metres(
                    &materials.water,
                    balance.radiative().equator_to_pole_kelvin,
                    EARTH_RADIUS_METRES,
                    EARTH_ROTATION_RADIANS_PER_SECOND,
                ),
                wind_metres_per_second: EARTH_WINDS.speed_metres_per_second,
            };
            let still = Currents::still(self.ocean.len());
            let currents = match flow {
                Flow::Still => still,
                Flow::GyresOnly => Currents {
                    gyre: ocean.currents().gyre,
                    ..still
                },
                Flow::Full => ocean.currents(),
            };
            balance.run(
                &Exchanged {
                    currents: &currents,
                    layer_humidity: &vec![self.layer_humidity; self.ocean.len()],
                    wetness: &vec![1.0; self.ocean.len()],
                    wind: Wind::Given(EARTH_WINDS.speed_metres_per_second),
                },
                None,
                &crate::parallel::Sequential,
            )
        }
    }

    #[derive(Clone, Copy)]
    enum Flow {
        Still,
        GyresOnly,
        Full,
    }

    const EARTH_RADIUS_METRES: f64 = 6.378e6;
    const EARTH_ROTATION_RADIANS_PER_SECOND: f64 = 7.292e-5;
    const EARTH_EDDY_DIFFUSIVITY_M2_PER_S: f64 = 1.0e6;
    const EARTH_LAYER_HUMIDITY: [f64; LAYERS] = [0.8, 0.5, 0.4];

    #[test]
    fn winds_and_currents_only_move_heat_so_nothing_crosses_the_south_pole() {
        let (width, height) = (36, 18);
        let planet = Planet::new(width, height, |x, _| {
            (width / 2..width * 2 / 3).contains(&x)
        });
        let cycle = planet.cycle();
        let strongest = cycle
            .northward
            .iter()
            .map(|t| t.atmosphere_watts.abs().max(t.ocean_watts.abs()))
            .fold(0.0, f64::max);
        let Some(pole) = cycle.northward.last() else {
            panic!("no transport reported");
        };
        assert!(strongest > 0.0);
        assert!(
            pole.atmosphere_watts.abs() < 1.0e-9 * strongest,
            "{pole:?} of {strongest}"
        );
        assert!(
            pole.ocean_watts.abs() < 1.0e-9 * strongest,
            "{pole:?} of {strongest}"
        );
    }

    #[test]
    fn gyres_carry_heat_from_the_subtropics_toward_the_poles_in_both_hemispheres() {
        let (width, height) = (36, 18);
        let planet = Planet::new(width, height, |x, _| x < width / 6);
        let (moving, still) = (planet.run(Flow::GyresOnly), planet.run(Flow::Still));
        let change =
            |row: u32| row_mean(width, row, &moving.mean) - row_mean(width, row, &still.mean);
        let (north_subtropics, south_subtropics) = (height / 2 - 3, height / 2 + 2);
        let (north_high, south_high) = (height / 6, height - 1 - height / 6);
        assert!(
            change(north_subtropics) < 0.0,
            "{}",
            change(north_subtropics)
        );
        assert!(
            change(south_subtropics) < 0.0,
            "{}",
            change(south_subtropics)
        );
        assert!(change(north_high) > 0.0, "{}", change(north_high));
        assert!(change(south_high) > 0.0, "{}", change(south_high));
    }

    #[test]
    fn an_aether_absorber_warms_until_its_extra_glow_matches_what_it_absorbs() {
        let (aether, celsius) = (50.0, 0.0);
        let kelvin = celsius + KELVIN_AT_ZERO_CELSIUS;
        let glow_per_kelvin = 4.0 * STEFAN_BOLTZMANN_W_PER_M2_PER_K4 * kelvin * kelvin * kelvin;
        for absorptance in [0.0, 0.5, 1.0] {
            let warming = aether_warming_kelvin(aether, absorptance, celsius);
            assert!(
                (warming * glow_per_kelvin - absorptance * aether).abs() < 1e-12,
                "absorptance {absorptance}: {warming} K"
            );
        }
    }

    #[test]
    fn a_companion_below_the_orbit_warms_the_south_more_than_the_north() {
        let (width, height) = (36, 18);
        let alone = Planet::new(width, height, |x, _| x < width / 3);
        let mut lit = Planet::new(width, height, |x, _| x < width / 3);
        let sun_watts = lit.heavens.star.luminosity_watts;
        lit.heavens.companions.push(Companion {
            name: "ember".into(),
            luminosity_watts: 1.0e3 * sun_watts,
            distance_metres: 400.0 * EARTH_ORBIT.semi_major_axis_metres,
            ecliptic_latitude_degrees: -60.0,
            ecliptic_longitude_degrees: 0.0,
            aether_watts: 0.0,
        });
        let (before, after) = (alone.cycle(), lit.cycle());
        let warming =
            |row: u32| row_mean(width, row, &after.mean) - row_mean(width, row, &before.mean);
        let (north, south) = (warming(0), warming(height - 1));
        assert!(south > north && north >= 0.0, "north {north} south {south}");
    }

    #[test]
    fn heat_carried_down_a_temperature_step_does_work_and_up_it_takes_work() {
        let (north, south) = (280.0, 300.0);
        let northward = north_south_work_watts(&[4.0, 0.0], &[north, south]);
        let middle = (north + south) / 2.0;
        assert!(
            (northward - 4.0 * (south - north) / middle).abs() < 1e-12,
            "{northward}"
        );
        assert_eq!(
            north_south_work_watts(&[-4.0, 0.0], &[north, south]),
            -northward
        );
    }

    #[test]
    fn moist_air_carries_its_vapors_latent_heat_between_cells() {
        let planet = Planet::new(4, 2, |_, _| false);
        let grid = Grid {
            width: planet.width,
            height: planet.height,
            elevation: &planet.zeros,
            celsius: &planet.zeros,
            condensation_celsius: &planet.zeros,
            ocean: &planet.ocean,
            sun_latitude: 0.0,
        };
        let sky = planet.heavens.sky(EARTH_ORBIT);
        let materials = test_materials();
        let balance = EnergyBalance {
            grid: &grid,
            cell_metres: EARTH_RADIUS_METRES,
            sky: &sky,
            materials: &materials,
            atmosphere_diffusivity_m2_per_s: EARTH_EDDY_DIFFUSIVITY_M2_PER_S,
        };
        let contrast = |relative_humidity: f64| {
            balance.moist_enthalpy(25.0, relative_humidity)
                - balance.moist_enthalpy(5.0, relative_humidity)
        };
        assert_eq!(contrast(0.0), 20.0);
        assert!(
            contrast(0.4) < contrast(0.8),
            "{} {}",
            contrast(0.4),
            contrast(0.8)
        );
    }

    #[test]
    fn the_ground_hands_the_air_earths_share_of_heat_by_evaporation_and_warm_contact() {
        let cycle = Planet::new(36, 18, |x, _| x * 10 < 36 * 3).cycle();
        let surface = cycle.surface;
        assert!(
            (60.0..110.0).contains(&surface.latent_w_per_m2),
            "{surface:?}"
        );
        assert!(
            (5.0..40.0).contains(&surface.sensible_w_per_m2),
            "{surface:?}"
        );
        assert!(
            (220.0..260.0).contains(&surface.to_space_w_per_m2),
            "{surface:?}"
        );
        let carried = surface.latent_w_per_m2 + surface.sensible_w_per_m2;
        assert!(
            (carried - surface.air_radiative_cooling_w_per_m2).abs()
                < 0.05 * surface.air_radiative_cooling_w_per_m2,
            "{surface:?}"
        );
    }

    #[test]
    fn the_flux_equator_is_where_the_atmosphere_stops_carrying_heat_north() {
        let gains = [1.0, 1.0, -1.0, -3.0, 1.0, 1.0];
        let crossing = flux_equator(&gains);
        let row_spacing = 2.0 / 6.0;
        let a_third_of_the_way_south_of_the_equator = row_spacing / 3.0;
        assert!(
            (crossing - a_third_of_the_way_south_of_the_equator).abs() < 1e-12,
            "{crossing}"
        );
    }

    #[test]
    fn the_flux_equator_follows_the_sun_into_the_summer_hemisphere() {
        let (width, height) = (36, 18);
        let planet = Planet::new(width, height, |x, _| x < width / 3);
        let cycle = planet.cycle();
        let at = |fraction: f64| cycle.flux_equator_at(&cycle.phase(fraction * cycle.year_seconds));
        let (june, december) = (0.47, 0.97);
        assert!(
            at(june) < at(december),
            "june {} december {}",
            at(june),
            at(december)
        );
        let swing_degrees = (at(december) - at(june)) * 90.0;
        assert!(swing_degrees < 2.0 * 23.44, "{swing_degrees}");
    }

    #[test]
    fn water_rising_at_the_equator_carries_the_warmth_of_the_subtropics_where_it_sank() {
        let flux = 1.0e7;
        let overturning = [0.0, 0.0, flux, 0.0, -flux, 0.0];
        let surface = [0.0, 10.0, 30.0, 30.0, 12.0, 0.0];
        let mut deep = [-100.0; 6];
        return_temperatures(&overturning, 1, &surface, &mut deep);
        assert_eq!(deep[2], surface[1], "{deep:?}");
        assert_eq!(deep[3], surface[4], "{deep:?}");
    }

    #[test]
    fn the_trades_cool_the_equator_by_lifting_water_that_sank_in_the_subtropics() {
        let (width, height) = (36, 18);
        let planet = Planet::new(width, height, |x, _| x < width / 6);
        let (full, gyres) = (planet.run(Flow::Full), planet.run(Flow::GyresOnly));
        for row in [height / 2 - 1, height / 2] {
            let change = row_mean(width, row, &full.mean) - row_mean(width, row, &gyres.mean);
            assert!(change < 0.0, "row {row}: {change}");
        }
    }

    fn row_mean(width: u32, row: u32, field: &[f64]) -> f64 {
        let start = (row * width) as usize;
        field[start..start + width as usize].iter().sum::<f64>() / f64::from(width)
    }

    #[test]
    fn an_earthlike_ocean_world_settles_warm_at_the_equator_and_below_freezing_at_the_poles() {
        let (width, height) = (36, 18);
        let planet = Planet::new(width, height, |_, _| false);
        let cycle = planet.cycle();
        assert!(cycle.settled, "settled after {} years", cycle.years);
        let equator = row_mean(width, height / 2, &cycle.mean);
        let pole = row_mean(width, 0, &cycle.mean);
        assert!((20.0..32.0).contains(&equator), "equator {equator}");
        assert!(pole < 0.0, "pole {pole}");
    }

    #[test]
    fn an_earthlike_planet_settles_near_earths_mean_temperature() {
        let (width, height) = (36, 18);
        let planet = Planet::new(width, height, |x, _| x * 10 < width * 3);
        let cycle = planet.cycle();
        assert!(cycle.settled, "settled after {} years", cycle.years);
        let (mut weighted, mut area) = (0.0, 0.0);
        for row in 0..height {
            let north = (0.5 - (f64::from(row) + 0.5) / f64::from(height)) * PI;
            let breadth = libm::cos(north);
            weighted += breadth * row_mean(width, row, &cycle.mean);
            area += breadth;
        }
        let global = weighted / area;
        assert!((12.0..18.0).contains(&global), "global mean {global}");
        let equator = row_mean(width, height / 2, &cycle.mean);
        assert!(equator > 20.0, "equator {equator}");
    }

    #[test]
    fn land_swings_with_the_seasons_more_than_the_sea() {
        let (width, height) = (36, 18);
        let planet = Planet::new(width, height, |x, _| x < width / 2);
        let cycle = planet.cycle();
        let row = height / 4;
        let swing = |x: u32| {
            let i = (row * width + x) as usize;
            cycle.warmest[i] - cycle.mean[i]
        };
        let (land, sea) = (swing(width / 4), swing(3 * width / 4));
        assert!(land > 2.0 * sea, "land {land} sea {sea}");
    }

    #[test]
    fn the_northern_summer_comes_half_a_year_after_the_southern() {
        let (width, height) = (36, 18);
        let planet = Planet::new(width, height, |_, _| true);
        let cycle = planet.cycle();
        let north = (height / 5 * width) as usize;
        let south = ((height - 1 - height / 5) * width) as usize;
        let at = |i: usize, fraction: f64| {
            let phase = cycle.phase(fraction * cycle.year_seconds);
            cycle.anomaly(i, &phase)
        };
        let june = 0.47;
        assert!(at(north, june) > 0.0 && at(south, june) < 0.0);
        assert!(at(north, june + 0.5) < 0.0 && at(south, june + 0.5) > 0.0);
    }
}
