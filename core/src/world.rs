// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec;
use alloc::vec::Vec;

use crate::Seed;
use crate::atmosphere::{LAYERS, TROPOSPHERE};
use crate::climate::{Blocks, Climate, ClimateRun};
use crate::cycle::{self, Circulation, Flux, Grid, Pools, Sea, Settled, Substance, Weather};
use crate::energy::{
    EnergyBalance, Exchanged, SETTLED_KELVIN, SeasonalCycle, Surface, Transport, Wind,
    aether_warming_kelvin,
};
use crate::hash::Purpose;
use crate::hydrology::{self, Drainage};
use crate::life::{Cover, Flora, Site};
use crate::material::Materials;
use crate::material::Phase;
use crate::noise::Field;
use crate::ocean::{Currents, Ocean, SVERDRUP_CUBIC_METRES_PER_SECOND, equatorial_width_metres};
use crate::parallel::Parallel;
use crate::scale::{
    DAYS_PER_YEAR, FREEZING_DEPRESSION_CELSIUS_PER_GRAM_PER_KG, LAPSE_RATE_CELSIUS_PER_METRE,
    SECONDS_PER_DAY, SECONDS_PER_HOUR, Scale,
};
use crate::sky::{Heavens, Orbit, Sky};
use crate::wind::{self, Belts, Planet, Winds};
use core::f64::consts::TAU;

pub const SHALLOWEST_SEA_METRES: f64 = 1.0;

pub struct Recipe {
    pub equator_km: f64,
    pub sidereal_day_hours: f64,
    pub orbit: Orbit,
    pub land_fraction: f64,
    pub continent_frequency: f64,
    pub roughness: f64,
    pub peak_height_metres: f64,
    pub deepest_metres: f64,
    pub weather_systems_around: f64,
    pub weather_celsius: f64,
    pub climate_downscale: u32,
    pub salt_residence_years: f64,
    pub eddy_diffusivity_m2_per_s: f64,
    pub ipse_age_days: f64,
    pub ipse_mean_amount: f64,
}

const WEATHER_LAYERS: u64 = 2;

impl Recipe {
    pub fn surface_metres(&self, terrain: Terrain) -> f64 {
        match terrain {
            Terrain::Land { height } => height * self.peak_height_metres,
            Terrain::Ocean { depth } => -(depth * self.deepest_metres).max(SHALLOWEST_SEA_METRES),
        }
    }
}

pub struct WorldSpec {
    pub seed: Seed,
    pub width: u32,
    pub height: u32,
    pub recipe: Recipe,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Terrain {
    Ocean { depth: f64 },
    Land { height: f64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SubstanceCell {
    pub vapor: f64,
    pub condensed: f64,
    pub flux: Flux,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Water {
    pub column_mm: f64,
    pub boundary_humidity: f64,
    pub rain_mm_per_year: f64,
    pub evaporation_mm_per_year: f64,
    pub river_mm_per_year: f64,
    pub supply_mm_per_year: f64,
    pub basin_salt_per_year: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cell {
    pub terrain: Terrain,
    pub elevation: f64,
    pub temperature: f64,
    pub summer_temperature: f64,
    pub aether_warming_kelvin: f64,
    pub salinity_grams_per_kg: f64,
    pub gyre_sverdrups: f64,
    pub water: Water,
    pub ipse: SubstanceCell,
    pub cover: Cover,
}

pub struct GlobalMap {
    pub width: u32,
    pub height: u32,
    pub scale: Scale,
    pub climate_scale: Scale,
    pub climate_steps: u32,
    pub energy_years: u32,
    pub energy_settled: bool,
    pub coupling: Coupling,
    pub radiation: Radiation,
    pub radiative_contrast: f64,
    pub northward_heat: Vec<Transport>,
    pub equatorial_width_metres: f64,
    pub sea_rain_mm_per_year: f64,
    pub sea_evaporation_mm_per_year: f64,
    pub itcz_north_degrees: Span,
    pub overturning_sverdrups: Vec<f64>,
    pub belts: Belts,
    pub sea_salinity_grams_per_kg: f64,
    pub cells: Vec<Cell>,
    pub climate_moved_per_step: Vec<f64>,
    pub ipse_moved_per_step: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Coupling {
    pub residuals: Vec<Residual>,
    pub settled: bool,
    pub wind_metres_per_second: f64,
    pub atmosphere_work_w_per_m2: f64,
}

struct FineTemperatures {
    mean: Vec<f64>,
    summer: Vec<f64>,
    condensation: Vec<f64>,
}

impl FineTemperatures {
    fn grid<'a>(&'a self, spec: &WorldSpec, elevation: &'a [f64], ocean: &'a [bool]) -> Grid<'a> {
        Grid {
            width: spec.width,
            height: spec.height,
            elevation,
            celsius: &self.mean,
            condensation_celsius: &self.condensation,
            ocean,
            sun_latitude: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Radiation {
    pub surface: Surface,
    pub rain_heating_w_per_m2: f64,
    pub layer_humidity: [f64; LAYERS],
    pub layer_rain_share: [f64; LAYERS],
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Residual {
    pub energy_years: u32,
    pub wind_metres_per_second: f64,
    pub atmosphere_work_w_per_m2: f64,
    pub kelvin: f64,
    pub wind: f64,
    pub humidity: f64,
}

const ROUNDS_WITHOUT_PROGRESS: usize = 4;
const LOOSEST_ROUND_KELVIN: f64 = 10.0 * SETTLED_KELVIN;
const ROUND_TOLERANCE_PER_CHANGE: f64 = 0.1;

const MOISTURE_STEP_AT_MOST: f64 = 0.5;
const MOISTURE_STEP_AT_LEAST: f64 = 1.0 / 16.0;

fn relax(given: &mut [f64], answered: &[f64], step: f64) {
    for (value, &next) in given.iter_mut().zip(answered) {
        *value += step * (next - *value);
    }
}

fn row_means(field: &[[f64; LAYERS]], width: usize) -> Vec<[f64; LAYERS]> {
    field
        .chunks(width)
        .flat_map(|row| {
            let mean: [f64; LAYERS] = core::array::from_fn(|layer| {
                row.iter().map(|cell| cell[layer]).sum::<f64>() / row.len() as f64
            });
            core::iter::repeat_n(mean, row.len())
        })
        .collect()
}

fn largest_change(before: &[f64], after: &[f64]) -> f64 {
    before
        .iter()
        .zip(after)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Span {
    pub mean: f64,
    pub lowest: f64,
    pub highest: f64,
}

fn itcz_through_the_year(seasons: &SeasonalCycle) -> Span {
    let days = libm::ceil(seasons.year_seconds / SECONDS_PER_DAY) as usize;
    let north_degrees = |day: usize| {
        let phase = seasons.phase((day as f64 + 0.5) * seasons.year_seconds / days as f64);
        -seasons.flux_equator_at(&phase) * DEGREES_PER_GRID_LATITUDE
    };
    let (mut lowest, mut highest) = (f64::INFINITY, f64::NEG_INFINITY);
    for day in 0..days {
        lowest = lowest.min(north_degrees(day));
        highest = highest.max(north_degrees(day));
    }
    Span {
        mean: -seasons.flux_equator * DEGREES_PER_GRID_LATITUDE,
        lowest,
        highest,
    }
}

const DEGREES_PER_GRID_LATITUDE: f64 = 90.0;

fn annual_aether(sky: &Sky, grid: &Grid) -> Vec<f64> {
    if sky.companions.iter().all(|c| c.aether_watts <= 0.0) {
        return vec![0.0; grid.height as usize];
    }
    let days = libm::ceil(sky.year_seconds() / SECONDS_PER_DAY) as usize;
    let day_seconds = sky.year_seconds() / days as f64;
    (0..grid.height)
        .map(|row| {
            let north = -grid.latitude(row) * core::f64::consts::FRAC_PI_2;
            (0..days)
                .map(|day| sky.daily_aether((day as f64 + 0.5) * day_seconds, north))
                .sum::<f64>()
                / days as f64
        })
        .collect()
}

impl GlobalMap {
    pub fn cell(&self, x: u32, y: u32) -> &Cell {
        &self.cells[(y * self.width + x) as usize]
    }
}

fn substance_cell(settled: &Settled, i: usize) -> SubstanceCell {
    SubstanceCell {
        vapor: settled.pools.vapor[i],
        condensed: settled.pools.condensed[i],
        flux: settled.flux[i],
    }
}

impl WorldSpec {
    pub fn elevation(&self, x: i64, y: u32) -> f64 {
        Field {
            seed: self.seed,
            purpose: Purpose::Elevation,
            base_frequency: self.recipe.continent_frequency,
            roughness: self.recipe.roughness,
            max_layers: u64::MAX,
        }
        .sample_cylinder(x, y, self.width)
    }

    fn weather_pattern(&self) -> Vec<f64> {
        let field = Field {
            seed: self.seed,
            purpose: Purpose::Weather,
            base_frequency: self.recipe.weather_systems_around,
            roughness: self.recipe.roughness,
            max_layers: WEATHER_LAYERS,
        };
        self.positions()
            .map(|(x, y)| field.sample_cylinder(x, y, self.width))
            .collect()
    }

    fn positions(&self) -> impl Iterator<Item = (i64, u32)> + '_ {
        (0..self.height).flat_map(|y| (0..self.width).map(move |x| (i64::from(x), y)))
    }

    fn terrain(&self, elevations: &[f64]) -> Vec<Terrain> {
        let areas: Vec<f64> = (0..elevations.len())
            .map(|i| cycle::row_breadth(i as u32 / self.width, self.height))
            .collect();
        let sea_level = quantile(elevations, &areas, 1.0 - self.recipe.land_fraction);
        let lowest = elevations.iter().copied().fold(f64::INFINITY, f64::min);
        let highest = elevations.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        elevations
            .iter()
            .map(|&e| {
                if e > sea_level {
                    Terrain::Land {
                        height: (e - sea_level) / (highest - sea_level),
                    }
                } else {
                    Terrain::Ocean {
                        depth: (sea_level - e) / (sea_level - lowest),
                    }
                }
            })
            .collect()
    }

    fn altitude_cooling(&self, terrain: Terrain) -> f64 {
        match terrain {
            Terrain::Land { height } => {
                LAPSE_RATE_CELSIUS_PER_METRE * height * self.recipe.peak_height_metres
            }
            Terrain::Ocean { .. } => 0.0,
        }
    }

    pub fn scale(&self, wind_metres_per_second: f64) -> Scale {
        Scale::new(
            self.recipe.equator_km,
            self.width,
            self.recipe.eddy_diffusivity_m2_per_s,
            wind_metres_per_second,
        )
    }

    fn drained(
        &self,
        grid: &Grid,
        climate: &Climate,
        terrain: &[Terrain],
        materials: &Materials,
    ) -> (Drainage, f64) {
        let runoff: Vec<f64> = (0..terrain.len())
            .map(|i| {
                if grid.ocean[i] {
                    0.0
                } else {
                    (climate.rain_mm_per_year[i] - climate.evaporation_mm_per_year[i]).max(0.0)
                }
            })
            .collect();
        let drainage = hydrology::drain(grid, &runoff, materials.rock.solubility);
        let sea_depths: Vec<f64> = terrain.iter().map(|&t| self.sea_depth_metres(t)).collect();
        let salinity = hydrology::steady_salinity_grams_per_kg(
            drainage.sea_salt_per_year,
            self.recipe.salt_residence_years,
            grid.total(&sea_depths),
        );
        (drainage, salinity)
    }

    fn fine_temperatures(
        &self,
        blocks: &Blocks,
        seasons: &SeasonalCycle,
        terrain: &[Terrain],
    ) -> FineTemperatures {
        let at_altitude = |sea_level: Vec<f64>| -> Vec<f64> {
            sea_level
                .iter()
                .zip(terrain)
                .map(|(&celsius, &t)| celsius - self.altitude_cooling(t))
                .collect()
        };
        let mean = at_altitude(blocks.upsample(&seasons.mean));
        let boundary_cooling = LAPSE_RATE_CELSIUS_PER_METRE * TROPOSPHERE[0].middle_metres();
        FineTemperatures {
            condensation: mean.iter().map(|&t| t - boundary_cooling).collect(),
            summer: at_altitude(blocks.upsample(&seasons.warmest)),
            mean,
        }
    }

    fn sea_depth_metres(&self, terrain: Terrain) -> f64 {
        match terrain {
            Terrain::Ocean { depth } => {
                (depth * self.recipe.deepest_metres).max(SHALLOWEST_SEA_METRES)
            }
            Terrain::Land { .. } => 0.0,
        }
    }

    pub fn generate<P: Parallel>(
        &self,
        materials: &Materials,
        flora: &Flora,
        heavens: &Heavens,
        parallel: &P,
    ) -> GlobalMap {
        let mut elevations = vec![0.0; (self.width * self.height) as usize];
        let mut rows: Vec<&mut [f64]> = elevations.chunks_mut(self.width as usize).collect();
        parallel.for_each(&mut rows, |y, row| {
            for (x, elevation) in row.iter_mut().enumerate() {
                *elevation = self.elevation(x as i64, y as u32);
            }
        });
        let terrain = self.terrain(&elevations);
        let ocean: Vec<bool> = terrain
            .iter()
            .map(|t| matches!(t, Terrain::Ocean { .. }))
            .collect();
        let recipe = &self.recipe;
        let blocks = Blocks::new(self.width, self.height, recipe.climate_downscale);
        let coarse_ocean = blocks.majority(&ocean);
        let coarse_flat = vec![0.0; coarse_ocean.len()];
        let sky = heavens.sky(recipe.orbit);
        let coarse = Grid {
            width: blocks.width as u32,
            height: blocks.height as u32,
            elevation: &coarse_flat,
            celsius: &coarse_flat,
            condensation_celsius: &coarse_flat,
            ocean: &coarse_ocean,
            sun_latitude: 0.0,
        };
        let energy = EnergyBalance {
            grid: &coarse,
            cell_metres: recipe.equator_km * 1000.0 / blocks.width as f64,
            sky: &sky,
            materials,
            atmosphere_diffusivity_m2_per_s: recipe.eddy_diffusivity_m2_per_s,
        };
        let radiative = energy.radiative();
        let radiative_contrast = radiative.contrast;
        let radius_metres = recipe.equator_km * 1000.0 / TAU;
        let rotation_radians_per_second = TAU / (recipe.sidereal_day_hours * SECONDS_PER_HOUR);
        let planet = Planet {
            radius_metres,
            rotation_radians_per_second,
            tropopause_metres: TROPOSPHERE[LAYERS - 1].top_metres,
            radiative_contrast,
        };
        let belts = Belts::of(&planet);
        let air_width =
            wind::equatorial_width_metres(&planet, &materials.air, radiative.mean_kelvin);
        let equatorial_width = equatorial_width_metres(
            &materials.water,
            radiative.equator_to_pole_kelvin,
            radius_metres,
            rotation_radians_per_second,
        );
        let driven_by = |wind_metres_per_second: f64| {
            Ocean {
                grid: &coarse,
                belts,
                radius_metres,
                rotation_radians_per_second,
                water_kg_per_cubic_metre: materials.water.density,
                equatorial_width_metres: equatorial_width,
                wind_metres_per_second,
            }
            .currents()
        };
        let pattern = self.weather_pattern();
        let precision = SETTLED_KELVIN / radiative.mean_kelvin;
        let mut humidity = vec![[1.0; LAYERS]; coarse_ocean.len()];
        let mut wetness = vec![1.0; coarse_ocean.len()];
        let mut tolerance = LOOSEST_ROUND_KELVIN;
        let mut moisture_step = MOISTURE_STEP_AT_MOST;
        let mut previous = energy.run_within(
            &Exchanged {
                currents: &Currents::still(coarse_ocean.len()),
                layer_humidity: &humidity,
                wetness: &wetness,
                wind: Wind::FromHeatEngine,
            },
            None,
            parallel,
            tolerance,
        );
        let mut speed = previous.wind_metres_per_second;
        let mut residuals: Vec<Residual> = Vec::new();
        let (seasons, climate, fine, winds, currents, settled) = loop {
            let currents = driven_by(speed);
            let seasons = energy.run_within(
                &Exchanged {
                    currents: &currents,
                    layer_humidity: &humidity,
                    wetness: &wetness,
                    wind: Wind::Given(speed),
                },
                Some(&previous),
                parallel,
                tolerance,
            );
            let winds = Winds {
                planet,
                equatorial_width_metres: air_width,
                speed_metres_per_second: speed,
            };
            let fine = self.fine_temperatures(&blocks, &seasons, &terrain);
            let climate = ClimateRun {
                fine: &fine.grid(self, &elevations, &ocean),
                downscale: recipe.climate_downscale,
                seasons: &seasons,
                belts,
                equator_km: recipe.equator_km,
                eddy_diffusivity_m2_per_s: recipe.eddy_diffusivity_m2_per_s,
                weather_pattern: &pattern,
                weather_celsius: recipe.weather_celsius,
                water: &materials.water,
                air: &materials.air,
                layer_humidity: &humidity,
                soil_capacity_mm: materials.rock.holds_water_mm,
                winds,
            }
            .run();
            let answered = wind::surface_wind_metres_per_second(wind::ground_friction_w_per_m2(
                seasons.atmosphere_work_w_per_m2,
                &materials.air,
            ));
            let gap = answered - speed;
            let residual = Residual {
                wind_metres_per_second: speed,
                energy_years: seasons.years,
                atmosphere_work_w_per_m2: seasons.atmosphere_work_w_per_m2,
                kelvin: largest_change(&previous.mean, &seasons.mean)
                    .max(largest_change(&previous.ground_mean, &seasons.ground_mean))
                    / SETTLED_KELVIN,
                wind: gap.abs() / answered.max(speed).max(f64::MIN_POSITIVE) / precision,
                humidity: largest_change(
                    humidity.as_flattened(),
                    row_means(&climate.coarse_layer_humidity, blocks.width).as_flattened(),
                ) / precision,
            };
            residuals.push(residual);
            let best_round = (0..residuals.len())
                .min_by(|&a, &b| residuals[a].kelvin.total_cmp(&residuals[b].kelvin))
                .unwrap_or(0);
            let stalled = residuals.len() - 1 - best_round >= ROUNDS_WITHOUT_PROGRESS;
            let next = (speed + answered) / 2.0;
            let settled = residual.kelvin <= 1.0 && tolerance <= SETTLED_KELVIN;
            if settled || stalled {
                break (seasons, climate, fine, winds, currents, settled);
            }
            if let [.., before, now] = residuals.as_slice() {
                moisture_step = if now.humidity > before.humidity {
                    (moisture_step / 2.0).max(MOISTURE_STEP_AT_LEAST)
                } else {
                    (moisture_step * 1.25).min(MOISTURE_STEP_AT_MOST)
                };
            }
            relax(
                humidity.as_flattened_mut(),
                row_means(&climate.coarse_layer_humidity, blocks.width).as_flattened(),
                moisture_step,
            );
            relax(&mut wetness, &climate.coarse_wetness, moisture_step);
            tolerance = (residual.kelvin * SETTLED_KELVIN * ROUND_TOLERANCE_PER_CHANGE)
                .clamp(SETTLED_KELVIN, LOOSEST_ROUND_KELVIN);
            previous = seasons;
            speed = next;
        };
        let coupling = Coupling {
            residuals,
            settled,
            wind_metres_per_second: winds.speed_metres_per_second,
            atmosphere_work_w_per_m2: seasons.atmosphere_work_w_per_m2,
        };
        let rain_kg_per_m2_per_s = coarse.total(&climate.coarse_rain_mm_per_year)
            / coarse.total(&vec![1.0; coarse.len()])
            / (DAYS_PER_YEAR * SECONDS_PER_DAY);
        let coarse_area = coarse.total(&vec![1.0; coarse.len()]);
        let radiation = Radiation {
            surface: seasons.surface,
            rain_heating_w_per_m2: materials.water.latent_heat * rain_kg_per_m2_per_s,
            layer_humidity: core::array::from_fn(|layer| {
                let field: Vec<f64> = climate
                    .coarse_layer_humidity
                    .iter()
                    .map(|layers| layers[layer])
                    .collect();
                coarse.total(&field) / coarse_area
            }),
            layer_rain_share: {
                let totals: [f64; LAYERS] = core::array::from_fn(|layer| {
                    let field: Vec<f64> = climate
                        .coarse_layer_rain_mm
                        .iter()
                        .map(|layers| layers[layer])
                        .collect();
                    coarse.total(&field)
                });
                let all: f64 = totals.iter().sum::<f64>().max(f64::MIN_POSITIVE);
                core::array::from_fn(|layer| totals[layer] / all)
            },
        };
        let FineTemperatures {
            mean: temperature,
            summer: summer_temperature,
            condensation,
        } = &fine;
        let grid = fine.grid(self, &elevations, &ocean);
        let aether_by_row = annual_aether(&sky, &grid);
        let aether_warming: Vec<f64> = temperature
            .iter()
            .enumerate()
            .map(|(i, &celsius)| {
                aether_warming_kelvin(
                    aether_by_row[i / self.width as usize],
                    materials.ipse.aether_absorptance,
                    celsius,
                )
            })
            .collect();
        let felt_by_ipse = |field: &[f64]| -> Vec<f64> {
            field
                .iter()
                .zip(&aether_warming)
                .map(|(t, w)| t + w)
                .collect()
        };
        let (ipse_celsius, ipse_condensation) =
            (felt_by_ipse(temperature), felt_by_ipse(condensation));
        let ipse_grid = Grid {
            celsius: &ipse_celsius,
            condensation_celsius: &ipse_condensation,
            ..grid
        };
        let gyre_sverdrups: Vec<f64> = blocks
            .upsample(&currents.gyre)
            .iter()
            .zip(&ocean)
            .map(|(&psi, &sea)| {
                if sea {
                    psi / SVERDRUP_CUBIC_METRES_PER_SECOND
                } else {
                    0.0
                }
            })
            .collect();
        let scale = self.scale(winds.speed_metres_per_second);
        let circulation = Circulation {
            gas_spread: scale.mixing_share,
            wind_carry: scale.wind_carry,
            age_steps: scale.steps(recipe.ipse_age_days),
            belts,
            winds,
        };
        let deposit_weights: Vec<f64> = terrain
            .iter()
            .map(|t| match t {
                Terrain::Land { height } => *height,
                Terrain::Ocean { .. } => 0.0,
            })
            .collect();
        let ipse = cycle::settle(
            &ipse_grid,
            &Substance {
                material: &materials.ipse,
                sea: Sea::SinksThrough,
                soil_capacity: 0.0,
                evaporation_limit: None,
            },
            &circulation,
            &Weather {
                pattern: &pattern,
                amplitude_celsius: recipe.weather_celsius,
            },
            Pools::with_condensed(cycle::deposit(
                &grid,
                &deposit_weights,
                recipe.ipse_mean_amount,
            )),
        );

        let sea_only = |field: &[f64]| -> Vec<f64> {
            field
                .iter()
                .zip(&ocean)
                .map(|(&v, &sea)| if sea { v } else { 0.0 })
                .collect()
        };
        let sea_area = grid.total(&sea_only(&vec![1.0; ocean.len()]));
        let sea_mean =
            |field: &[f64]| grid.total(&sea_only(field)) / sea_area.max(f64::MIN_POSITIVE);
        let sea_rain = sea_mean(&climate.rain_mm_per_year);
        let sea_evaporation = sea_mean(&climate.evaporation_mm_per_year);
        let (drainage, sea_salinity) = self.drained(&grid, &climate, &terrain, materials);

        let cells = (0..terrain.len())
            .map(|i| {
                let land = !ocean[i];
                let rain = if land {
                    climate.rain_mm_per_year[i]
                } else {
                    0.0
                };
                let river = drainage.inflow_mm_per_year[i];
                let water = Water {
                    column_mm: climate.column_water_mm[i],
                    boundary_humidity: climate.boundary_humidity[i],
                    rain_mm_per_year: rain,
                    evaporation_mm_per_year: if land {
                        climate.evaporation_mm_per_year[i]
                    } else {
                        0.0
                    },
                    river_mm_per_year: river,
                    supply_mm_per_year: rain + river,
                    basin_salt_per_year: drainage.basin_salt_per_year[i],
                };
                let salinity_grams_per_kg = if land { 0.0 } else { sea_salinity };
                let ipse = substance_cell(&ipse, i);
                let summer_temperature = summer_temperature[i];
                let site = Surroundings {
                    terrain: terrain[i],
                    celsius: temperature[i],
                    summer_celsius: summer_temperature,
                    water_mm_per_year: water.supply_mm_per_year,
                    salinity_grams_per_kg,
                };
                Cell {
                    terrain: terrain[i],
                    elevation: elevations[i],
                    temperature: temperature[i],
                    summer_temperature,
                    aether_warming_kelvin: aether_warming[i],
                    salinity_grams_per_kg,
                    gyre_sverdrups: gyre_sverdrups[i],
                    cover: site.cover(&water, &ipse, materials, flora),
                    water,
                    ipse,
                }
            })
            .collect();

        GlobalMap {
            width: self.width,
            height: self.height,
            scale,
            climate_scale: climate.scale,
            climate_steps: climate.steps,
            energy_years: seasons.years,
            energy_settled: seasons.settled,
            coupling,
            radiative_contrast,
            northward_heat: seasons.northward.clone(),
            equatorial_width_metres: equatorial_width,
            sea_rain_mm_per_year: sea_rain,
            sea_evaporation_mm_per_year: sea_evaporation,
            itcz_north_degrees: itcz_through_the_year(&seasons),
            radiation,
            overturning_sverdrups: currents.overturning_sverdrups(blocks.width),
            belts,
            sea_salinity_grams_per_kg: sea_salinity,
            cells,
            climate_moved_per_step: climate.moved_per_step,
            ipse_moved_per_step: ipse.moved_per_step,
        }
    }
}

struct Surroundings {
    terrain: Terrain,
    celsius: f64,
    summer_celsius: f64,
    water_mm_per_year: f64,
    salinity_grams_per_kg: f64,
}

impl Surroundings {
    fn sea_freezes(&self, materials: &Materials) -> bool {
        let depression = FREEZING_DEPRESSION_CELSIUS_PER_GRAM_PER_KG * self.salinity_grams_per_kg;
        materials.water.pure_phase(self.celsius + depression) == Phase::Solid
    }

    fn cover(
        &self,
        water: &Water,
        ipse: &SubstanceCell,
        materials: &Materials,
        flora: &Flora,
    ) -> Cover {
        if matches!(self.terrain, Terrain::Ocean { .. }) {
            return if self.sea_freezes(materials) {
                Cover::SeaIce
            } else {
                Cover::Sea
            };
        }
        if water.rain_mm_per_year > 0.0
            && materials.water.pure_phase(self.summer_celsius) == Phase::Solid
        {
            return Cover::Ice;
        }
        let site = Site {
            celsius: self.celsius,
            water_mm_per_year: self.water_mm_per_year,
            ipse: ipse.vapor + ipse.condensed,
        };
        flora.dominant(&site).map_or(Cover::Barren, Cover::Plant)
    }
}

fn quantile(values: &[f64], weights: &[f64], fraction: f64) -> f64 {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
    let below = fraction * weights.iter().sum::<f64>();
    let mut covered = 0.0;
    for &i in &order {
        covered += weights[i];
        if covered > below {
            return values[i];
        }
    }
    order.last().map_or(0.0, |&i| values[i])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::earth;
    use crate::parallel::Sequential;
    use crate::scale::SECONDS_PER_DAY;

    fn generate(seed: u64) -> &'static GlobalMap {
        type Worlds = std::sync::Mutex<
            alloc::collections::BTreeMap<u64, &'static std::sync::OnceLock<GlobalMap>>,
        >;
        static WORLDS: std::sync::OnceLock<Worlds> = std::sync::OnceLock::new();
        let slot = *WORLDS
            .get_or_init(Default::default)
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .entry(seed)
            .or_insert_with(|| Box::leak(Box::default()));
        slot.get_or_init(|| {
            spec(seed).generate(
                &earth::materials(),
                &earth::flora(),
                &earth::heavens(),
                &Sequential,
            )
        })
    }

    fn spec(seed: u64) -> WorldSpec {
        WorldSpec {
            seed: Seed(seed),
            width: 64,
            height: 32,
            recipe: earth::recipe(),
        }
    }

    fn land(map: &GlobalMap) -> impl Iterator<Item = &Cell> {
        map.cells
            .iter()
            .filter(|c| matches!(c.terrain, Terrain::Land { .. }))
    }

    #[test]
    fn the_seed_determines_the_world() {
        let elevations = |seed| -> Vec<f64> {
            let spec = spec(seed);
            spec.positions()
                .map(|(x, y)| spec.elevation(x, y))
                .collect()
        };
        assert_eq!(elevations(1), elevations(1));
        assert_ne!(elevations(1), elevations(2));
    }

    #[test]
    fn the_world_wraps_east_west() {
        let spec = spec(7);
        for (x, y) in [(0, 0), (5, 10), (63, 31)] {
            assert_eq!(spec.elevation(x, y), spec.elevation(x + 64, y));
            assert_eq!(spec.elevation(x, y), spec.elevation(x - 64, y));
        }
    }

    #[test]
    fn land_fraction_sets_the_share_of_the_surface_that_is_land() {
        let map = generate(3);
        let area = |i: usize| cycle::row_breadth(i as u32 / map.width, map.height);
        let total: f64 = (0..map.cells.len()).map(area).sum();
        let land: f64 = (0..map.cells.len())
            .filter(|&i| matches!(map.cells[i].terrain, Terrain::Land { .. }))
            .map(area)
            .sum();
        let largest_cell = area((map.height / 2 * map.width) as usize);
        assert!(
            (land - 0.3 * total).abs() <= largest_cell,
            "{land} of {total}"
        );
    }

    #[test]
    fn an_earthlike_planet_puts_its_tropical_cell_edge_near_thirty_degrees() {
        let map = generate(3);
        let degrees = map.belts.hadley_edge * 90.0;
        assert!((26.0..36.0).contains(&degrees), "{degrees}");
        let polar_front = map.belts.polar_front * 90.0;
        assert!(polar_front > degrees && polar_front < 90.0);
    }

    #[test]
    fn a_slower_spinning_planet_has_wider_tropical_cells() {
        let mut slow = spec(3);
        slow.recipe.sidereal_day_hours *= 3.0;
        let slow_map = slow.generate(
            &earth::materials(),
            &earth::flora(),
            &earth::heavens(),
            &Sequential,
        );
        assert!(slow_map.belts.hadley_edge > generate(3).belts.hadley_edge);
    }

    #[test]
    fn the_wind_settles_where_ground_friction_spends_the_north_south_work() {
        let coupling = &generate(3).coupling;
        let answered = wind::surface_wind_metres_per_second(wind::ground_friction_w_per_m2(
            coupling.atmosphere_work_w_per_m2,
            &earth::materials().air,
        ));
        let gap = (answered - coupling.wind_metres_per_second).abs() / answered;
        assert!(
            gap < 0.01,
            "{answered} vs {}",
            coupling.wind_metres_per_second
        );
    }

    #[test]
    fn the_climate_runs_for_one_orbit() {
        let map = generate(3);
        let year = earth::heavens().sky(earth::ORBIT).year_seconds() / SECONDS_PER_DAY;
        let days = map.climate_scale.days(map.climate_steps);
        assert!(days >= year);
        assert!(map.climate_scale.days(map.climate_steps - 1) < year);
    }

    #[test]
    fn sea_ice_forms_exactly_where_sea_water_freezes() {
        let materials = earth::materials();
        let map = generate(5);
        let mut frozen = 0;
        for cell in map
            .cells
            .iter()
            .filter(|c| matches!(c.terrain, Terrain::Ocean { .. }))
        {
            let depression =
                FREEZING_DEPRESSION_CELSIUS_PER_GRAM_PER_KG * cell.salinity_grams_per_kg;
            let freezes = materials.water.pure_phase(cell.temperature + depression) == Phase::Solid;
            assert_eq!(cell.cover == Cover::SeaIce, freezes);
            frozen += usize::from(freezes);
        }
        assert!(frozen > 0, "salinity {}", map.sea_salinity_grams_per_kg);
    }

    #[test]
    fn salt_lowers_the_freezing_point_of_the_sea() {
        let materials = earth::materials();
        let sea = |salinity_grams_per_kg| Surroundings {
            terrain: Terrain::Ocean { depth: 0.5 },
            celsius: -1.0,
            summer_celsius: -1.0,
            water_mm_per_year: 0.0,
            salinity_grams_per_kg,
        };
        assert!(sea(0.0).sea_freezes(&materials));
        assert!(!sea(35.0).sea_freezes(&materials));
    }

    #[test]
    fn land_stays_icebound_only_where_even_summer_freezes() {
        let materials = earth::materials();
        let flora = earth::flora();
        let land = |summer_celsius| Surroundings {
            terrain: Terrain::Land { height: 0.5 },
            celsius: -5.0,
            summer_celsius,
            water_mm_per_year: 300.0,
            salinity_grams_per_kg: 0.0,
        };
        let snowy = Water {
            rain_mm_per_year: 300.0,
            ..Water::default()
        };
        let none = SubstanceCell::default();
        assert_eq!(
            land(-1.0).cover(&snowy, &none, &materials, &flora),
            Cover::Ice
        );
        assert_ne!(
            land(5.0).cover(&snowy, &none, &materials, &flora),
            Cover::Ice
        );
    }

    #[test]
    fn the_whole_sea_shares_one_salinity() {
        let map = generate(5);
        assert!(map.sea_salinity_grams_per_kg > 0.0);
        assert!(
            map.cells
                .iter()
                .filter(|c| matches!(c.terrain, Terrain::Ocean { .. }))
                .all(|c| c.salinity_grams_per_kg == map.sea_salinity_grams_per_kg)
        );
    }

    #[test]
    fn bare_land_carries_the_tallest_plant_that_survives_there() {
        let flora = earth::flora();
        let map = generate(11);
        let mut covered = 0;
        for cell in land(map).filter(|c| c.cover != Cover::Ice) {
            let site = Site {
                celsius: cell.temperature,
                water_mm_per_year: cell.water.supply_mm_per_year,
                ipse: cell.ipse.vapor + cell.ipse.condensed,
            };
            let expected = flora.dominant(&site).map_or(Cover::Barren, Cover::Plant);
            assert_eq!(cell.cover, expected);
            covered += usize::from(matches!(cell.cover, Cover::Plant(_)));
        }
        assert!(covered > 0);
    }

    #[test]
    fn plants_get_rain_and_the_rivers_that_pass_them() {
        let map = generate(9);
        let mut on_rivers = 0;
        for cell in land(map) {
            let w = cell.water;
            assert_eq!(
                w.supply_mm_per_year,
                w.rain_mm_per_year + w.river_mm_per_year
            );
            on_rivers += usize::from(w.river_mm_per_year > 0.0);
        }
        assert!(on_rivers > 0);
    }
}
