// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec;
use alloc::vec::Vec;

use crate::cycle::{Circulation, Grid, Weather, gather_mixing};
use crate::energy::SeasonalCycle;
use crate::material::{Material, Phase};
use crate::radiation;
use crate::scale::{
    GRAVITY_METRES_PER_SECOND_SQUARED, LAPSE_RATE_CELSIUS_PER_METRE,
    SEA_LEVEL_AIR_KG_PER_CUBIC_METRE, SURFACE_MOISTURE_TRANSFER_COEFFICIENT, Scale,
};
use crate::wind;

pub const LAYERS: usize = 3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Layer {
    pub base_metres: f64,
    pub top_metres: f64,
    pub base_hpa: f64,
    pub top_hpa: f64,
}

impl Layer {
    pub fn thickness_metres(&self) -> f64 {
        self.top_metres - self.base_metres
    }

    pub fn middle_metres(&self) -> f64 {
        (self.base_metres + self.top_metres) / 2.0
    }

    pub fn air_mass(&self) -> f64 {
        self.base_hpa - self.top_hpa
    }
}

pub const TROPOSPHERE: [Layer; LAYERS] = [
    Layer {
        base_metres: 0.0,
        top_metres: 1500.0,
        base_hpa: 1013.25,
        top_hpa: 845.6,
    },
    Layer {
        base_metres: 1500.0,
        top_metres: 5000.0,
        base_hpa: 845.6,
        top_hpa: 540.5,
    },
    Layer {
        base_metres: 5000.0,
        top_metres: 12_000.0,
        base_hpa: 540.5,
        top_hpa: 194.0,
    },
];

const PASCALS_PER_HECTOPASCAL: f64 = 100.0;

pub fn layer_air_kg_per_cubic_metre(layer: usize) -> f64 {
    TROPOSPHERE[layer].air_mass() * PASCALS_PER_HECTOPASCAL
        / GRAVITY_METRES_PER_SECOND_SQUARED
        / TROPOSPHERE[layer].thickness_metres()
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Slab {
    pub layer: usize,
    pub base_metres: f64,
    pub top_metres: f64,
    pub base_hpa: f64,
    pub top_hpa: f64,
}

impl Slab {
    pub fn thickness_metres(&self) -> f64 {
        self.top_metres - self.base_metres
    }

    pub fn middle_metres(&self) -> f64 {
        (self.base_metres + self.top_metres) / 2.0
    }

    pub fn air_kg_per_m2(&self) -> f64 {
        (self.base_hpa - self.top_hpa) * PASCALS_PER_HECTOPASCAL / GRAVITY_METRES_PER_SECOND_SQUARED
    }

    pub fn air_kg_per_cubic_metre(&self) -> f64 {
        self.air_kg_per_m2() / self.thickness_metres()
    }
}

pub const SLABS: usize = 28;

pub const RADIATION_SLABS: [Slab; SLABS] = [
    Slab {
        layer: 0,
        base_metres: 0.0,
        top_metres: 88.1,
        base_hpa: 1013.25,
        top_hpa: 1002.72,
    },
    Slab {
        layer: 0,
        base_metres: 88.1,
        top_metres: 179.4,
        base_hpa: 1002.72,
        top_hpa: 991.89,
    },
    Slab {
        layer: 0,
        base_metres: 179.4,
        top_metres: 274.2,
        base_hpa: 991.89,
        top_hpa: 980.75,
    },
    Slab {
        layer: 0,
        base_metres: 274.2,
        top_metres: 372.8,
        base_hpa: 980.75,
        top_hpa: 969.27,
    },
    Slab {
        layer: 0,
        base_metres: 372.8,
        top_metres: 475.5,
        base_hpa: 969.27,
        top_hpa: 957.43,
    },
    Slab {
        layer: 0,
        base_metres: 475.5,
        top_metres: 582.7,
        base_hpa: 957.43,
        top_hpa: 945.20,
    },
    Slab {
        layer: 0,
        base_metres: 582.7,
        top_metres: 694.7,
        base_hpa: 945.20,
        top_hpa: 932.55,
    },
    Slab {
        layer: 0,
        base_metres: 694.7,
        top_metres: 812.0,
        base_hpa: 932.55,
        top_hpa: 919.45,
    },
    Slab {
        layer: 0,
        base_metres: 812.0,
        top_metres: 935.2,
        base_hpa: 919.45,
        top_hpa: 905.86,
    },
    Slab {
        layer: 0,
        base_metres: 935.2,
        top_metres: 1064.8,
        base_hpa: 905.86,
        top_hpa: 891.73,
    },
    Slab {
        layer: 0,
        base_metres: 1064.8,
        top_metres: 1201.5,
        base_hpa: 891.73,
        top_hpa: 877.02,
    },
    Slab {
        layer: 0,
        base_metres: 1201.5,
        top_metres: 1346.3,
        base_hpa: 877.02,
        top_hpa: 861.67,
    },
    Slab {
        layer: 0,
        base_metres: 1346.3,
        top_metres: 1500.0,
        base_hpa: 861.67,
        top_hpa: 845.60,
    },
    Slab {
        layer: 1,
        base_metres: 1500.0,
        top_metres: 1677.4,
        base_hpa: 845.60,
        top_hpa: 827.37,
    },
    Slab {
        layer: 1,
        base_metres: 1677.4,
        top_metres: 1868.5,
        base_hpa: 827.37,
        top_hpa: 808.10,
    },
    Slab {
        layer: 1,
        base_metres: 1868.5,
        top_metres: 2075.5,
        base_hpa: 808.10,
        top_hpa: 787.63,
    },
    Slab {
        layer: 1,
        base_metres: 2075.5,
        top_metres: 2301.2,
        base_hpa: 787.63,
        top_hpa: 765.79,
    },
    Slab {
        layer: 1,
        base_metres: 2301.2,
        top_metres: 2549.4,
        base_hpa: 765.79,
        top_hpa: 742.34,
    },
    Slab {
        layer: 1,
        base_metres: 2549.4,
        top_metres: 2825.0,
        base_hpa: 742.34,
        top_hpa: 716.99,
    },
    Slab {
        layer: 1,
        base_metres: 2825.0,
        top_metres: 3134.7,
        base_hpa: 716.99,
        top_hpa: 689.35,
    },
    Slab {
        layer: 1,
        base_metres: 3134.7,
        top_metres: 3488.1,
        base_hpa: 689.35,
        top_hpa: 658.87,
    },
    Slab {
        layer: 1,
        base_metres: 3488.1,
        top_metres: 3899.2,
        base_hpa: 658.87,
        top_hpa: 624.78,
    },
    Slab {
        layer: 1,
        base_metres: 3899.2,
        top_metres: 4390.5,
        base_hpa: 624.78,
        top_hpa: 585.94,
    },
    Slab {
        layer: 1,
        base_metres: 4390.5,
        top_metres: 5000.0,
        base_hpa: 585.94,
        top_hpa: 540.50,
    },
    Slab {
        layer: 2,
        base_metres: 5000.0,
        top_metres: 5703.2,
        base_hpa: 540.50,
        top_hpa: 491.83,
    },
    Slab {
        layer: 2,
        base_metres: 5703.2,
        top_metres: 6671.7,
        base_hpa: 491.83,
        top_hpa: 430.67,
    },
    Slab {
        layer: 2,
        base_metres: 6671.7,
        top_metres: 8217.2,
        base_hpa: 430.67,
        top_hpa: 345.92,
    },
    Slab {
        layer: 2,
        base_metres: 8217.2,
        top_metres: 12000.0,
        base_hpa: 345.92,
        top_hpa: 194.00,
    },
];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Profile {
    pub celsius: [f64; LAYERS],
    pub convecting: bool,
}

fn stirred_celsius(air_celsius: f64, metres: f64) -> f64 {
    air_celsius - LAPSE_RATE_CELSIUS_PER_METRE * metres
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlabAir {
    pub celsius: [f64; SLABS],
    pub saturation: [f64; SLABS],
}

pub fn slab_air(
    guesses: Option<&[f64; SLABS]>,
    air_celsius: f64,
    boundary_relative_humidity: f64,
    air: &Material,
    water: &Material,
) -> SlabAir {
    let carried = carried_heat(air_celsius, boundary_relative_humidity, air, water);
    let mut celsius = [0.0; SLABS];
    let mut saturation = [0.0; SLABS];
    let mut below = None;
    for (slab, at) in RADIATION_SLABS.iter().enumerate() {
        let (metres, density) = (at.middle_metres(), at.air_kg_per_cubic_metre());
        let stirred = stirred_celsius(air_celsius, metres);
        let (stirred_heat, stirred_saturation) =
            saturated_heat(stirred, metres, density, air, water);
        if stirred_heat < carried {
            let guess = guesses.map(|known| known[slab].max(stirred)).or(below);
            let (parcel, parcel_saturation) = solve_rising(
                guess.unwrap_or(stirred),
                carried,
                metres,
                density,
                air,
                water,
            );
            celsius[slab] = parcel;
            saturation[slab] = parcel_saturation;
            below = Some(parcel);
        } else {
            celsius[slab] = stirred;
            saturation[slab] = stirred_saturation;
            below = None;
        }
    }
    SlabAir {
        celsius,
        saturation,
    }
}

fn solve_rising(
    start: f64,
    carried: f64,
    metres: f64,
    density: f64,
    air: &Material,
    water: &Material,
) -> (f64, f64) {
    let mut celsius = start;
    let mut last_step = f64::INFINITY;
    loop {
        let (saturation, per_kelvin) = water.saturation_with_slope(celsius);
        let excess = air.specific_heat * celsius
            + GRAVITY_METRES_PER_SECOND_SQUARED * metres
            + water.latent_heat * saturation / density
            - carried;
        let step = excess / (air.specific_heat + water.latent_heat * per_kelvin / density);
        if step == 0.0 || step.abs() >= last_step {
            return (celsius, saturation);
        }
        celsius -= step;
        last_step = step.abs();
    }
}

pub fn moist_parcel_celsius(
    air_celsius: f64,
    boundary_relative_humidity: f64,
    layer: usize,
    air: &Material,
    water: &Material,
) -> f64 {
    let metres = TROPOSPHERE[layer].middle_metres();
    let dry = air_celsius - GRAVITY_METRES_PER_SECOND_SQUARED / air.specific_heat * metres;
    let carried = carried_heat(air_celsius, boundary_relative_humidity, air, water);
    let (root, _) = solve_rising(
        dry,
        carried,
        metres,
        layer_air_kg_per_cubic_metre(layer),
        air,
        water,
    );
    root.max(dry)
}

pub fn profile(
    air_celsius: f64,
    boundary_relative_humidity: f64,
    air: &Material,
    water: &Material,
) -> Profile {
    profile_from(None, air_celsius, boundary_relative_humidity, air, water).0
}

pub fn profile_from(
    guesses: Option<&[f64; LAYERS]>,
    air_celsius: f64,
    boundary_relative_humidity: f64,
    air: &Material,
    water: &Material,
) -> (Profile, [f64; LAYERS]) {
    let carried = carried_heat(air_celsius, boundary_relative_humidity, air, water);
    let mut celsius = [0.0; LAYERS];
    let mut rising = [false; LAYERS];
    let mut below = None;
    for layer in 0..LAYERS {
        let (metres, density) = (
            TROPOSPHERE[layer].middle_metres(),
            layer_air_kg_per_cubic_metre(layer),
        );
        let stirred = stirred_celsius(air_celsius, metres);
        if saturated_heat(stirred, metres, density, air, water).0 < carried {
            let guess = guesses.map(|known| known[layer].max(stirred)).or(below);
            celsius[layer] = solve_rising(
                guess.unwrap_or(stirred),
                carried,
                metres,
                density,
                air,
                water,
            )
            .0;
            rising[layer] = true;
            below = Some(celsius[layer]);
        } else {
            celsius[layer] = stirred;
            below = None;
        }
    }
    (
        Profile {
            celsius,
            convecting: rising[1],
        },
        celsius,
    )
}

fn carried_heat(
    air_celsius: f64,
    boundary_relative_humidity: f64,
    air: &Material,
    water: &Material,
) -> f64 {
    air.specific_heat * air_celsius
        + water.latent_heat * boundary_relative_humidity * water.saturation(air_celsius)
            / air.density
}

fn saturated_heat(
    celsius: f64,
    metres: f64,
    density: f64,
    air: &Material,
    water: &Material,
) -> (f64, f64) {
    let saturation = water.saturation(celsius);
    (
        air.specific_heat * celsius
            + GRAVITY_METRES_PER_SECOND_SQUARED * metres
            + water.latent_heat * saturation / density,
        saturation,
    )
}

pub fn evaporation_share(scale: &Scale) -> f64 {
    let per_second = SEA_LEVEL_AIR_KG_PER_CUBIC_METRE
        * SURFACE_MOISTURE_TRANSFER_COEFFICIENT
        * scale.wind_metres_per_second
        / boundary_air_kg_per_square_metre();
    1.0 - libm::exp(-per_second * scale.step_seconds)
}

pub struct Ground {
    pub soil_capacity_mm: f64,
    pub evaporation_share: f64,
}

pub struct Air {
    pub water_mm: [Vec<f64>; LAYERS],
    pub soil_mm: Vec<f64>,
    pub snow_mm: Vec<f64>,
}

impl Air {
    pub fn total(&self, grid: &Grid) -> f64 {
        self.water_mm
            .iter()
            .map(|layer| grid.total(layer))
            .sum::<f64>()
            + grid.total(&self.soil_mm)
            + grid.total(&self.snow_mm)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Exchange {
    pub precipitation_mm: f64,
    pub evaporation_mm: f64,
    pub runoff_mm: f64,
    pub boundary_humidity: f64,
    pub layer_humidity: [f64; LAYERS],
    pub layer_rain_mm: [f64; LAYERS],
    pub wetness: f64,
}

pub struct Outcome {
    pub air: Air,
    pub per_step: Vec<Exchange>,
    pub moved_per_step: Vec<f64>,
}

struct Temperatures {
    surface: Vec<f64>,
    ground: Vec<f64>,
    layers: [Vec<f64>; LAYERS],
    convecting: Vec<bool>,
    parcels: Vec<Option<[f64; LAYERS]>>,
}

impl Temperatures {
    fn new(cells: usize) -> Self {
        Self {
            surface: vec![0.0; cells],
            ground: vec![0.0; cells],
            layers: core::array::from_fn(|_| vec![0.0; cells]),
            convecting: vec![false; cells],
            parcels: vec![None; cells],
        }
    }
}

struct RowDrift {
    east: f64,
    east_step: i8,
    south: f64,
    south_step: i8,
}

impl RowDrift {
    fn new(east: f64, south_fraction: f64, wind_carry: f64, row: usize, height: usize) -> Self {
        let step = |v: f64| match v.partial_cmp(&0.0) {
            Some(core::cmp::Ordering::Greater) => 1,
            Some(core::cmp::Ordering::Less) => -1,
            _ => 0,
        };
        let south_step = step(south_fraction);
        let target_exists = match south_step {
            1 => row + 1 < height,
            -1 => row > 0,
            _ => false,
        };
        let east_step = step(east);
        Self {
            east: if east_step == 0 {
                0.0
            } else {
                wind_carry * (1.0 - south_fraction.abs())
            },
            east_step,
            south: if target_exists {
                wind_carry * south_fraction.abs()
            } else {
                0.0
            },
            south_step: if target_exists { south_step } else { 0 },
        }
    }
}

pub(crate) fn boundary_air_kg_per_square_metre() -> f64 {
    TROPOSPHERE[0].air_mass() * PASCALS_PER_HECTOPASCAL / GRAVITY_METRES_PER_SECOND_SQUARED
}

#[derive(Clone, Copy)]
pub struct Seasons<'a> {
    pub cycle: &'a SeasonalCycle,
    pub step_seconds: f64,
    pub air: &'a Material,
    pub layer_humidity: &'a [[f64; LAYERS]],
}

pub struct Atmosphere<'a> {
    pub grid: &'a Grid<'a>,
    pub water: &'a Material,
    pub circulation: Circulation,
    pub weather: &'a Weather<'a>,
    pub seasons: Option<Seasons<'a>>,
    pub ground: Ground,
}

impl Atmosphere<'_> {
    fn south_fraction(&self, layer: usize, row: u32, belt_latitude: f64) -> f64 {
        let boundary = wind::boundary_south_fraction(
            self.grid.latitude(row),
            belt_latitude,
            &self.circulation.belts,
            &self.circulation.winds,
        );
        match layer {
            0 => boundary,
            1 => 0.0,
            _ => {
                let returning =
                    boundary * TROPOSPHERE[0].air_mass() / TROPOSPHERE[layer].air_mass();
                -returning.clamp(-1.0, 1.0)
            }
        }
    }
    fn capacity_mm(&self, layer: usize, celsius: f64) -> f64 {
        self.water.saturation(celsius) * TROPOSPHERE[layer].thickness_metres()
    }

    fn column(&self, i: usize, surface_celsius: f64) -> Profile {
        self.column_from(None, i, surface_celsius).0
    }

    fn column_from(
        &self,
        guesses: Option<&[f64; LAYERS]>,
        i: usize,
        surface_celsius: f64,
    ) -> (Profile, Option<[f64; LAYERS]>) {
        match self.seasons {
            Some(seasons) => {
                let (column, parcels) = profile_from(
                    guesses,
                    surface_celsius,
                    seasons.layer_humidity[i][0],
                    seasons.air,
                    self.water,
                );
                (column, Some(parcels))
            }
            None => (
                Profile {
                    celsius: core::array::from_fn(|layer| {
                        surface_celsius
                            - LAPSE_RATE_CELSIUS_PER_METRE * TROPOSPHERE[layer].middle_metres()
                    }),
                    convecting: false,
                },
                None,
            ),
        }
    }

    fn season(&self, step: u32, scratch: &mut Temperatures) {
        let grid = self.grid;
        scratch.surface.copy_from_slice(grid.celsius);
        scratch.ground.copy_from_slice(grid.celsius);
        let Some(seasons) = self.seasons else {
            return;
        };
        let cycle = seasons.cycle;
        let phase = cycle.phase(f64::from(step) * seasons.step_seconds);
        for (i, (surface, ground)) in scratch
            .surface
            .iter_mut()
            .zip(scratch.ground.iter_mut())
            .enumerate()
        {
            *surface += cycle.anomaly(i, &phase);
            *ground += cycle.ground_mean[i] - cycle.mean[i] + cycle.ground_anomaly(i, &phase);
        }
    }

    fn season_and_belt(&self, step: u32, scratch: &mut Temperatures) -> f64 {
        self.season(step, scratch);
        self.seasons.as_ref().map_or(0.0, |seasons| {
            let phase = seasons.cycle.phase(f64::from(step) * seasons.step_seconds);
            seasons.cycle.flux_equator_at(&phase)
        })
    }

    fn temperatures(&self, step: u32, scratch: &mut Temperatures) {
        let width = self.grid.width as usize;
        self.weather.apply(
            &scratch.surface,
            width,
            f64::from(step) * self.circulation.wind_carry,
            &mut scratch.layers[0],
        );
        for i in 0..scratch.surface.len() {
            let anomaly = scratch.surface[i] - scratch.layers[0][i];
            let (column, parcels) =
                self.column_from(scratch.parcels[i].as_ref(), i, scratch.surface[i]);
            scratch.parcels[i] = parcels;
            for layer in 1..LAYERS {
                scratch.layers[layer][i] = column.celsius[layer] - anomaly;
            }
            scratch.convecting[i] = column.convecting;
        }
    }

    pub fn initial(&self) -> Air {
        let cells = self.grid.len();
        let water_mm = core::array::from_fn(|layer| {
            (0..cells)
                .map(|i| {
                    let surface = self.grid.celsius[i];
                    let celsius = if layer == 0 {
                        surface
                    } else {
                        self.column(i, surface).celsius[layer]
                    };
                    self.capacity_mm(layer, celsius)
                })
                .collect()
        });
        Air {
            water_mm,
            soil_mm: (0..cells)
                .map(|i| {
                    if self.grid.ocean[i] {
                        0.0
                    } else {
                        self.ground.soil_capacity_mm
                    }
                })
                .collect(),
            snow_mm: vec![0.0; cells],
        }
    }

    fn surface(&self, i: usize, temps: &Temperatures, air: &mut Air, exchange: &mut Exchange) {
        let boundary_room = (self.capacity_mm(0, temps.ground[i]) - air.water_mm[0][i]).max(0.0);
        let evaporable = boundary_room * self.ground.evaporation_share;
        if self.grid.ocean[i] {
            air.water_mm[0][i] += evaporable;
            exchange.evaporation_mm += evaporable;
            return;
        }
        let thawed = self.water.pure_phase(temps.ground[i]) != Phase::Solid;
        let melt = if thawed {
            core::mem::take(&mut air.snow_mm[i])
        } else {
            0.0
        };
        self.soak(i, melt, air, exchange);
        let capacity = self.ground.soil_capacity_mm;
        if capacity > 0.0 {
            let wetness = air.soil_mm[i] / capacity;
            let breathed = air.soil_mm[i].min(evaporable * wetness);
            air.soil_mm[i] -= breathed;
            air.water_mm[0][i] += breathed;
            exchange.evaporation_mm += breathed;
        }
    }

    fn convect(
        &self,
        i: usize,
        temps: &Temperatures,
        convergence: f64,
        air: &mut Air,
        exchange: &mut Exchange,
    ) {
        let Some(seasons) = self.seasons else {
            return;
        };
        let air_mass = |layer: usize| {
            TROPOSPHERE[layer].air_mass() * PASCALS_PER_HECTOPASCAL
                / GRAVITY_METRES_PER_SECOND_SQUARED
        };
        let static_energy = |layer: usize| {
            seasons.air.specific_heat * temps.layers[layer][i]
                + GRAVITY_METRES_PER_SECOND_SQUARED * TROPOSPHERE[layer].middle_metres()
        };
        let rising_kg_per_m2_per_s = convergence * air_mass(0) / seasons.step_seconds;
        let lifted_heating = -rising_kg_per_m2_per_s * (static_energy(2) - static_energy(0));
        let needed = (seasons.cycle.free_air_cooling_w_per_m2[i] - lifted_heating).max(0.0);
        let top = LAYERS - 1;
        if !temps.convecting[i] {
            let warming_per_kg = static_energy(top) - static_energy(0);
            if needed <= 0.0 || warming_per_kg <= 0.0 {
                return;
            }
            let sinking_kg_per_m2 = needed / warming_per_kg * seasons.step_seconds;
            let leaving: [f64; LAYERS] = core::array::from_fn(|layer| {
                if layer == 0 {
                    0.0
                } else {
                    (sinking_kg_per_m2 / air_mass(layer)).min(1.0) * air.water_mm[layer][i]
                }
            });
            for layer in 0..top {
                air.water_mm[layer][i] += leaving[layer + 1] - leaving[layer];
            }
            air.water_mm[top][i] -= leaving[top];
            return;
        }
        let middle = top - 1;
        let lapse = (temps.layers[middle][i] - temps.layers[top][i])
            / (TROPOSPHERE[top].middle_metres() - TROPOSPHERE[middle].middle_metres());
        let outflow_celsius = temps.layers[top][i]
            - lapse * (TROPOSPHERE[top].top_metres - TROPOSPHERE[top].middle_metres());
        let holds_per_kg = |layer: usize| {
            let celsius = if layer == top {
                outflow_celsius
            } else {
                temps.layers[layer][i]
            };
            self.capacity_mm(layer, celsius) / air_mass(layer)
        };
        let condensing = air.water_mm[0][i] / air_mass(0) - holds_per_kg(top);
        if needed <= 0.0 || condensing <= 0.0 {
            return;
        }
        let flux_kg_per_m2_per_s = needed / (self.water.latent_heat * condensing);
        let moved_kg_per_m2 = (flux_kg_per_m2_per_s * seasons.step_seconds).min(air_mass(0));
        let leaving: [f64; LAYERS] = core::array::from_fn(|layer| {
            (moved_kg_per_m2 / air_mass(layer)).min(1.0) * air.water_mm[layer][i]
        });
        let mut carried = leaving[0];
        let mut fallen = 0.0;
        for layer in 1..LAYERS {
            let holds = moved_kg_per_m2 * holds_per_kg(layer);
            let condensed = (carried - holds).max(0.0);
            carried -= condensed;
            fallen += condensed;
            exchange.layer_rain_mm[layer] += condensed;
        }
        air.water_mm[0][i] += leaving[1] - leaving[0];
        for layer in 1..top {
            air.water_mm[layer][i] += leaving[layer + 1] - leaving[layer];
        }
        air.water_mm[top][i] += carried - leaving[top];
        self.land(i, fallen, temps, air, exchange);
    }

    fn soak(&self, i: usize, water_mm: f64, air: &mut Air, exchange: &mut Exchange) {
        let room = (self.ground.soil_capacity_mm - air.soil_mm[i]).max(0.0);
        let soaked = water_mm.min(room);
        air.soil_mm[i] += soaked;
        exchange.runoff_mm += water_mm - soaked;
    }

    fn precipitate(&self, i: usize, temps: &Temperatures, air: &mut Air, exchange: &mut Exchange) {
        let mut fallen = 0.0;
        for layer in 0..LAYERS {
            let capacity = self.capacity_mm(layer, temps.layers[layer][i]);
            let excess = (air.water_mm[layer][i] - capacity).max(0.0);
            air.water_mm[layer][i] -= excess;
            fallen += excess;
            exchange.layer_rain_mm[layer] += excess;
        }
        self.land(i, fallen, temps, air, exchange);
    }

    fn land(
        &self,
        i: usize,
        fallen: f64,
        temps: &Temperatures,
        air: &mut Air,
        exchange: &mut Exchange,
    ) {
        exchange.precipitation_mm += fallen;
        if self.grid.ocean[i] {
            return;
        }
        if self.water.pure_phase(temps.surface[i]) == Phase::Solid {
            air.snow_mm[i] += fallen;
        } else {
            self.soak(i, fallen, air, exchange);
        }
    }

    fn convergence(&self, belt_latitude: f64) -> Vec<f64> {
        let grid = self.grid;
        let carry = self.circulation.wind_carry;
        let height = grid.height as usize;
        let meridional: Vec<(f64, f64)> = (0..grid.height)
            .map(|row| {
                let fraction = self.south_fraction(0, row, belt_latitude);
                (fraction, carry * fraction.abs())
            })
            .collect();
        (0..height)
            .map(|y| {
                let (south, out) = meridional[y];
                let leaves = if (south > 0.0 && y + 1 < height) || (south < 0.0 && y > 0) {
                    out
                } else {
                    0.0
                };
                let per_area = |n: usize| grid.breadth(n as u32) / grid.breadth(y as u32);
                let from_north = if y > 0 && meridional[y - 1].0 > 0.0 {
                    meridional[y - 1].1 * per_area(y - 1)
                } else {
                    0.0
                };
                let from_south = if y + 1 < height && meridional[y + 1].0 < 0.0 {
                    meridional[y + 1].1 * per_area(y + 1)
                } else {
                    0.0
                };
                from_north + from_south - leaves
            })
            .collect()
    }

    fn rise_and_sink(&self, i: usize, convergence: f64, air: &mut Air) {
        let mass = |layer: usize| TROPOSPHERE[layer].air_mass();
        if convergence > 0.0 {
            let up_from_boundary = (convergence).min(1.0) * air.water_mm[0][i];
            let up_from_middle = (convergence * mass(0) / mass(1)).min(1.0) * air.water_mm[1][i];
            air.water_mm[0][i] -= up_from_boundary;
            air.water_mm[1][i] += up_from_boundary - up_from_middle;
            air.water_mm[2][i] += up_from_middle;
        } else if convergence < 0.0 {
            let sinking = -convergence;
            let down_from_middle = (sinking * mass(0) / mass(1)).min(1.0) * air.water_mm[1][i];
            let down_from_top = (sinking * mass(0) / mass(2)).min(1.0) * air.water_mm[2][i];
            air.water_mm[2][i] -= down_from_top;
            air.water_mm[1][i] += down_from_top - down_from_middle;
            air.water_mm[0][i] += down_from_middle;
        }
    }

    fn carry(&self, belt_latitude: f64, air: &mut Air) -> f64 {
        let grid = self.grid;
        let width = grid.width as usize;
        let height = grid.height as usize;
        let share = self.circulation.gas_spread / 8.0;
        let mut moved = 0.0;
        let mut gathered = vec![0.0; width];
        for layer in 0..LAYERS {
            let rows: Vec<RowDrift> = (0..height)
                .map(|y| {
                    let flow = wind::layer_flow(
                        layer,
                        grid.latitude(y as u32) - belt_latitude,
                        &self.circulation.belts,
                    );
                    RowDrift::new(
                        flow.east,
                        self.south_fraction(layer, y as u32, belt_latitude),
                        self.circulation.wind_carry,
                        y,
                        height,
                    )
                })
                .collect();
            let water: Vec<f64> = air.water_mm[layer].iter().map(|&mm| mm.max(0.0)).collect();
            let next = &mut air.water_mm[layer];
            for y in 0..height {
                let neighbors = gather_mixing(grid, &water, y, &mut gathered);
                let row = &water[y * width..(y + 1) * width];
                let leaving = share * neighbors + rows[y].east + rows[y].south;
                let arriving_from =
                    |n: usize| rows[n].south * grid.breadth(n as u32) / grid.breadth(y as u32);
                let from_north = y
                    .checked_sub(1)
                    .filter(|&n| rows[n].south_step == 1)
                    .map(|n| (&water[n * width..(n + 1) * width], arriving_from(n)));
                let from_south = (y + 1 < height && rows[y + 1].south_step == -1).then(|| {
                    (
                        &water[(y + 1) * width..(y + 2) * width],
                        arriving_from(y + 1),
                    )
                });
                for x in 0..width {
                    let own = row[x];
                    let out = own * leaving;
                    let mut arrived = share * (gathered[x] - own);
                    let upwind = match rows[y].east_step {
                        1 => Some(if x == 0 { width - 1 } else { x - 1 }),
                        -1 => Some(if x + 1 == width { 0 } else { x + 1 }),
                        _ => None,
                    };
                    if let Some(source) = upwind {
                        arrived += row[source] * rows[y].east;
                    }
                    if let Some((north, fraction)) = from_north {
                        arrived += north[x] * fraction;
                    }
                    if let Some((south, fraction)) = from_south {
                        arrived += south[x] * fraction;
                    }
                    next[y * width + x] = own - out + arrived;
                    moved += out;
                }
            }
        }
        moved
    }

    fn step(
        &self,
        step: u32,
        air: &mut Air,
        exchange: &mut [Exchange],
        scratch: &mut Temperatures,
    ) -> f64 {
        let belt_latitude = self.season_and_belt(step, scratch);
        let convergence = self.convergence(belt_latitude);
        self.temperatures(step, scratch);
        let width = self.grid.width as usize;
        for (i, exchange) in exchange.iter_mut().enumerate() {
            self.surface(i, scratch, air, exchange);
            self.convect(i, scratch, convergence[i / width], air, exchange);
            self.precipitate(i, scratch, air, exchange);
            self.rise_and_sink(i, convergence[i / width], air);
            exchange.boundary_humidity +=
                air.water_mm[0][i] / self.capacity_mm(0, scratch.surface[i]);
            let capacity_celsius = radiation::capacity_celsius(
                scratch.surface[i],
                &core::array::from_fn(|layer| scratch.layers[layer][i]),
            );
            for (layer, held) in exchange.layer_humidity.iter_mut().enumerate() {
                *held += air.water_mm[layer][i] / self.capacity_mm(layer, capacity_celsius[layer]);
            }
            let capacity = self.ground.soil_capacity_mm;
            exchange.wetness += if self.grid.ocean[i] {
                1.0
            } else if capacity > 0.0 {
                air.soil_mm[i] / capacity
            } else {
                0.0
            };
        }
        self.carry(belt_latitude, air)
    }

    pub fn run(&self, steps: u32) -> Outcome {
        self.run_from(self.initial(), steps)
    }

    pub fn run_from(&self, mut air: Air, steps: u32) -> Outcome {
        let cells = self.grid.len();
        let mut totals = vec![Exchange::default(); cells];
        let mut scratch = Temperatures::new(cells);
        let mut moved_per_step = Vec::with_capacity(steps as usize);
        for s in 0..steps {
            moved_per_step.push(self.step(s, &mut air, &mut totals, &mut scratch));
        }
        let per_step = f64::from(steps.max(1));
        for t in &mut totals {
            t.precipitation_mm /= per_step;
            t.evaporation_mm /= per_step;
            t.runoff_mm /= per_step;
            t.boundary_humidity /= per_step;
            t.wetness /= per_step;
            for held in t.layer_humidity.iter_mut().chain(&mut t.layer_rain_mm) {
                *held /= per_step;
            }
        }
        Outcome {
            air,
            per_step: totals,
            moved_per_step,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::earth::water;
    use crate::scale::SECONDS_PER_DAY;

    const CIRCULATION: Circulation = Circulation {
        gas_spread: 0.08,
        wind_carry: 0.25,
        age_steps: 0,
        belts: crate::wind::Belts::EARTHLIKE,
        winds: wind::EARTH_WINDS,
    };

    struct World {
        width: u32,
        height: u32,
        elevation: Vec<f64>,
        celsius: Vec<f64>,
        ocean: Vec<bool>,
    }

    impl World {
        fn new(width: u32, height: u32, ocean_columns: u32, celsius: f64) -> Self {
            let cells = (width * height) as usize;
            let ocean: Vec<bool> = (0..cells)
                .map(|i| (i as u32 % width) < ocean_columns)
                .collect();
            Self {
                width,
                height,
                elevation: ocean.iter().map(|&o| if o { -1.0 } else { 0.0 }).collect(),
                celsius: vec![celsius; cells],
                ocean,
            }
        }

        fn grid(&self) -> Grid<'_> {
            Grid {
                width: self.width,
                height: self.height,
                elevation: &self.elevation,
                celsius: &self.celsius,
                condensation_celsius: &self.celsius,
                ocean: &self.ocean,
                sun_latitude: 0.0,
            }
        }
    }

    fn atmosphere<'a>(grid: &'a Grid<'a>, water: &'a Material, ground: Ground) -> Atmosphere<'a> {
        Atmosphere {
            grid,
            water,
            circulation: CIRCULATION,
            weather: &Weather::CALM,
            seasons: None,
            ground,
        }
    }

    const EARTH_MEAN_KELVIN: f64 = 288.0;

    #[test]
    fn the_atmospheric_equatorial_width_matches_earths_observed_scale() {
        let width_km = wind::equatorial_width_metres(
            &wind::EARTH_WINDS.planet,
            &crate::earth::air(),
            EARTH_MEAN_KELVIN,
        ) / 1000.0;
        assert!((1000.0..1600.0).contains(&width_km), "{width_km} km");
    }

    #[test]
    fn the_trades_cross_isobars_at_earths_observed_speed() {
        let world = World::new(8, 36, 8, 25.0);
        let grid = world.grid();
        let material = water();
        let atmosphere = atmosphere(
            &grid,
            &material,
            Ground {
                soil_capacity_mm: 0.0,
                evaporation_share: 0.0,
            },
        );
        let middle_of_the_trades = 14;
        let speed = atmosphere.south_fraction(0, middle_of_the_trades, 0.0)
            * wind::EARTH_WINDS.speed_metres_per_second;
        assert!((0.3..3.0).contains(&speed), "{speed} m/s");
    }

    #[test]
    fn layers_follow_the_standard_atmosphere() {
        let total: f64 = TROPOSPHERE.iter().map(Layer::air_mass).sum();
        assert_eq!(total, 1013.25 - 194.0);
    }

    #[test]
    fn moist_warm_air_warms_the_upper_layers_and_cold_air_keeps_the_stirred_lapse() {
        let (air, water) = (crate::earth::air(), water());
        let stirred =
            |surface: f64| surface - LAPSE_RATE_CELSIUS_PER_METRE * TROPOSPHERE[1].middle_metres();
        let tropical = profile(27.0, 0.8, &air, &water);
        assert!(tropical.convecting);
        assert!(tropical.celsius[1] > stirred(27.0), "{tropical:?}");
        let polar = profile(-20.0, 0.8, &air, &water);
        assert!(!polar.convecting);
        assert_eq!(polar.celsius[1], stirred(-20.0));
    }

    #[test]
    fn a_lifted_parcel_keeps_its_moist_static_energy() {
        let (air, water) = (crate::earth::air(), water());
        let (surface, relative) = (25.0, 0.8);
        let parcel = moist_parcel_celsius(surface, relative, 2, &air, &water);
        let carried = air.specific_heat * surface
            + water.latent_heat * relative * water.saturation(surface) / air.density;
        let aloft = air.specific_heat * parcel
            + GRAVITY_METRES_PER_SECOND_SQUARED * TROPOSPHERE[2].middle_metres()
            + water.latent_heat * water.saturation(parcel) / layer_air_kg_per_cubic_metre(2);
        assert!(
            (aloft - carried).abs() < 1e-6 * carried,
            "{aloft} vs {carried}"
        );
    }

    #[test]
    fn colder_layers_hold_less_water() {
        let world = World::new(4, 4, 0, 20.0);
        let grid = world.grid();
        let material = water();
        let air = atmosphere(
            &grid,
            &material,
            Ground {
                soil_capacity_mm: 0.0,
                evaporation_share: 0.0,
            },
        )
        .initial();
        let per_metre =
            |layer: usize| air.water_mm[layer][0] / TROPOSPHERE[layer].thickness_metres();
        assert!(per_metre(0) > per_metre(1));
        assert!(per_metre(1) > per_metre(2));
    }

    #[test]
    fn the_recorded_exchanges_account_for_every_change_in_the_air() {
        let world = World::new(8, 12, 4, 25.0);
        let grid = world.grid();
        let material = water();
        let atmosphere = atmosphere(
            &grid,
            &material,
            Ground {
                soil_capacity_mm: 150.0,
                evaporation_share: 0.5,
            },
        );
        let start = atmosphere.initial();
        let before = start.total(&grid);
        let steps = 60;
        let outcome = atmosphere.run_from(start, steps);
        let summed = |pick: fn(&Exchange) -> f64, sea: bool| {
            let field: Vec<f64> = outcome
                .per_step
                .iter()
                .zip(&world.ocean)
                .map(|(e, &o)| if o == sea { pick(e) } else { 0.0 })
                .collect();
            grid.total(&field) * f64::from(steps)
        };
        let sea_evaporation = summed(|e| e.evaporation_mm, true);
        let gained =
            sea_evaporation - summed(|e| e.precipitation_mm, true) - summed(|e| e.runoff_mm, false);
        let after = outcome.air.total(&grid);
        assert!(sea_evaporation > 0.0);
        assert!(
            (after - before - gained).abs() <= before * 1e-9,
            "{before} + {gained} != {after}"
        );
    }

    #[test]
    fn water_is_conserved_when_nothing_leaves() {
        let world = World::new(8, 12, 0, 15.0);
        let grid = world.grid();
        let material = water();
        let atmosphere = atmosphere(
            &grid,
            &material,
            Ground {
                soil_capacity_mm: 1.0e9,
                evaporation_share: 0.5,
            },
        );
        let mut start = atmosphere.initial();
        start.soil_mm.fill(0.0);
        let before = start.total(&grid);
        let after = atmosphere.run_from(start, 60).air.total(&grid);
        assert!(
            (after - before).abs() <= before * 1e-12,
            "{before} -> {after}"
        );
    }

    #[test]
    fn converging_air_lifts_water_into_higher_layers() {
        let world = World::new(4, 4, 0, 20.0);
        let grid = world.grid();
        let material = water();
        let atmosphere = atmosphere(
            &grid,
            &material,
            Ground {
                soil_capacity_mm: 0.0,
                evaporation_share: 0.0,
            },
        );
        let mut air = atmosphere.initial();
        let before = [air.water_mm[0][0], air.water_mm[1][0], air.water_mm[2][0]];
        atmosphere.rise_and_sink(0, 0.1, &mut air);
        assert!(air.water_mm[0][0] < before[0]);
        assert!(air.water_mm[2][0] > before[2]);
        let total = |a: &Air| a.water_mm[0][0] + a.water_mm[1][0] + a.water_mm[2][0];
        assert!((total(&air) - before.iter().sum::<f64>()).abs() < 1e-12);
    }

    #[test]
    fn sinking_air_arrives_dry() {
        let world = World::new(4, 4, 0, 20.0);
        let grid = world.grid();
        let material = water();
        let atmosphere = atmosphere(
            &grid,
            &material,
            Ground {
                soil_capacity_mm: 0.0,
                evaporation_share: 0.0,
            },
        );
        let mut air = atmosphere.initial();
        let capacity = air.water_mm[0][0];
        air.water_mm[0][0] = 0.0;
        atmosphere.rise_and_sink(0, -0.1, &mut air);
        let gained = air.water_mm[0][0];
        assert!(gained > 0.0);
        assert!(
            gained < 0.1 * capacity,
            "air from aloft brings little water: {gained} of {capacity}"
        );
    }

    #[test]
    fn the_sea_evaporates_in_proportion_to_how_dry_the_air_is() {
        let world = World::new(1, 1, 1, 25.0);
        let grid = world.grid();
        let material = water();
        let share = evaporation_share(&Scale {
            step_seconds: SECONDS_PER_DAY,
            ..Scale::new(
                40_000.0,
                400,
                0.0,
                wind::EARTH_WINDS.speed_metres_per_second,
            )
        });
        let atmosphere = atmosphere(
            &grid,
            &material,
            Ground {
                soil_capacity_mm: 0.0,
                evaporation_share: share,
            },
        );
        let warm = one_cell(25.0, 25.0);
        let saturated = atmosphere.initial();
        let capacity = saturated.water_mm[0][0];
        let evaporated_into = |humidity: f64| {
            let mut air = atmosphere.initial();
            air.water_mm[0][0] = capacity * humidity;
            atmosphere.surface(0, &warm, &mut air, &mut Exchange::default());
            air.water_mm[0][0] - capacity * humidity
        };
        let humid = evaporated_into(0.9);
        let dry = evaporated_into(0.6);
        assert!((dry / humid - 4.0).abs() < 1e-9, "humid {humid} dry {dry}");
        let typical = evaporated_into(0.8);
        assert!(
            (3.0..8.0).contains(&typical),
            "a 20% deficit at 25 °C evaporates a few mm a day, as on Earth: {typical}"
        );
    }

    #[test]
    fn snow_waits_for_the_thaw() {
        let mut world = World::new(1, 1, 0, -5.0);
        let material = water();
        let ground = || Ground {
            soil_capacity_mm: 100.0,
            evaporation_share: 0.0,
        };
        let mut air = {
            let grid = world.grid();
            let mut air = atmosphere(&grid, &material, ground()).initial();
            air.soil_mm[0] = 0.0;
            air
        };
        let frozen = one_cell(-5.0, -60.0);
        let mut exchange = Exchange::default();
        {
            let grid = world.grid();
            atmosphere(&grid, &material, ground()).precipitate(0, &frozen, &mut air, &mut exchange);
        }
        let snow = air.snow_mm[0];
        assert!(snow > 0.0);
        assert_eq!(air.soil_mm[0], 0.0);

        world.celsius[0] = 5.0;
        let thawed = one_cell(5.0, 5.0);
        let grid = world.grid();
        atmosphere(&grid, &material, ground()).surface(0, &thawed, &mut air, &mut exchange);
        assert_eq!(air.snow_mm[0], 0.0);
        assert_eq!(air.soil_mm[0], snow);
    }

    fn one_cell(surface: f64, layers: f64) -> Temperatures {
        Temperatures {
            surface: vec![surface],
            ground: vec![surface],
            layers: core::array::from_fn(|_| vec![layers]),
            convecting: vec![false],
            parcels: vec![None],
        }
    }

    fn belt_at_step(step: u32, flux_equator: f64, swing: f64) -> f64 {
        let (width, height) = (8, 32);
        let world = World::new(width, height, width, 20.0);
        let cells = (width * height) as usize;
        let cycle = SeasonalCycle {
            mean: world.celsius.clone(),
            warmest: world.celsius.clone(),
            harmonics: vec![[0.0; 4]; cells],
            ground_mean: world.celsius.clone(),
            ground_harmonics: vec![[0.0; 4]; cells],
            year_seconds: 4.0,
            years: 1,
            settled: true,
            northward: Vec::new(),
            flux_equator,
            flux_equator_harmonics: [swing, 0.0, 0.0, 0.0],
            atmosphere_work_w_per_m2: 0.0,
            wind_metres_per_second: wind::EARTH_WINDS.speed_metres_per_second,
            surface: crate::energy::Surface::default(),
            free_air_cooling_w_per_m2: vec![0.0; cells],
            end_of_year: world.celsius.clone(),
            ground_end_of_year: world.celsius.clone(),
        };
        let grid = world.grid();
        let material = water();
        let air = crate::earth::air();
        let humidity = vec![[0.8; LAYERS]; cells];
        let atmosphere = Atmosphere {
            seasons: Some(Seasons {
                cycle: &cycle,
                step_seconds: 1.0,
                air: &air,
                layer_humidity: &humidity,
            }),
            ..atmosphere(
                &grid,
                &material,
                Ground {
                    soil_capacity_mm: 0.0,
                    evaporation_share: 0.0,
                },
            )
        };
        let mut scratch = Temperatures::new(grid.len());
        atmosphere.season_and_belt(step, &mut scratch)
    }

    #[test]
    fn wind_belts_follow_the_energy_flux_equator_through_the_year() {
        let (mean, swing) = (-0.02, -0.1);
        assert!((belt_at_step(0, mean, swing) - (mean + swing)).abs() < 1e-12);
        assert!((belt_at_step(2, mean, swing) - (mean - swing)).abs() < 1e-12);
    }

    #[test]
    fn the_hadley_cell_wets_the_equator_and_dries_thirty_degrees() {
        let world = World::new(24, 24, 6, 25.0);
        let grid = world.grid();
        let material = water();
        let atmosphere = atmosphere(
            &grid,
            &material,
            Ground {
                soil_capacity_mm: 150.0,
                evaporation_share: 0.2,
            },
        );
        let outcome = atmosphere.run(400);
        let land_rain = |rows: &[u32]| -> f64 {
            rows.iter()
                .flat_map(|&row| (6..24).map(move |x| (row * 24 + x) as usize))
                .map(|i| outcome.per_step[i].precipitation_mm)
                .sum::<f64>()
                / rows.len() as f64
        };
        let degrees = |row: u32| (grid.latitude(row) * 90.0).abs();
        assert!(degrees(11) < 5.0 && degrees(7) > 25.0 && degrees(7) < 35.0);
        let equator = land_rain(&[11, 12]);
        let subtropics = land_rain(&[7, 16]);
        assert!(
            equator > 2.0 * subtropics,
            "equator {equator} subtropics {subtropics}"
        );
    }
}
