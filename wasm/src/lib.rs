// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use orbis::Seed;
use orbis::atlas::Atlas;
use orbis::earth;
use orbis::life::Cover;
use orbis::parallel::Sequential;
use orbis::relief::{self, Field, Level, Rock};
use orbis::wind;
use orbis::world::{GlobalMap, Recipe, Terrain, WorldSpec};
use wasm_bindgen::prelude::*;

const REGION_REFINEMENT: usize = 8;
const NARROWEST_RIVER_METRES: f64 = 1.0;

fn spec(seed: u64, width: u32, height: u32) -> WorldSpec {
    WorldSpec {
        seed: Seed(seed),
        width,
        height,
        recipe: earth::recipe(),
    }
}

fn floats(values: impl Iterator<Item = f64>) -> Vec<f32> {
    values.map(|value| value as f32).collect()
}

/// The ground alone: available at once, without waiting for the climate.
#[wasm_bindgen]
pub struct Shape {
    seed: Seed,
    recipe: Recipe,
    terrains: Vec<Terrain>,
    atlas: Atlas,
    rock: Rock,
    relief: Option<Field>,
}

#[wasm_bindgen]
impl Shape {
    #[wasm_bindgen(constructor)]
    pub fn new(seed: u64, width: u32, height: u32) -> Shape {
        let spec = spec(seed, width, height);
        let terrains = spec.terrains();
        let atlas = Atlas::of_terrain(width, height, &terrains, &spec.recipe);
        Shape {
            seed: spec.seed,
            rock: Rock::of(&earth::materials().rock),
            recipe: spec.recipe,
            terrains,
            atlas,
            relief: None,
        }
    }

    #[wasm_bindgen(getter)]
    pub fn width(&self) -> u32 {
        self.atlas.width
    }

    #[wasm_bindgen(getter)]
    pub fn height(&self) -> u32 {
        self.atlas.height
    }

    /// Surface height of each cell in metres, negative under the sea.
    pub fn elevation(&self) -> Vec<f32> {
        floats(
            self.terrains
                .iter()
                .map(|&terrain| self.recipe.surface_metres(terrain)),
        )
    }

    /// The eroded planet-wide height field, `reliefWidth × reliefHeight`, row 0 north.
    pub fn relief(&mut self) -> Vec<f32> {
        floats(self.planet().land.height.iter().copied())
    }

    /// Area draining through each relief cell, in square metres.
    #[wasm_bindgen(js_name = drainageArea)]
    pub fn drainage_area(&mut self) -> Vec<f32> {
        let land = &self.planet().land;
        floats(
            land.drainage(&vec![1.0; land.height.len()])
                .area_m2
                .into_iter(),
        )
    }

    #[wasm_bindgen(getter, js_name = reliefWidth)]
    pub fn relief_width(&self) -> u32 {
        self.atlas.width * relief::PLANET_REFINEMENT as u32
    }

    #[wasm_bindgen(getter, js_name = reliefHeight)]
    pub fn relief_height(&self) -> u32 {
        self.atlas.height * relief::PLANET_REFINEMENT as u32
    }

    #[wasm_bindgen(getter, js_name = circumferenceMetres)]
    pub fn circumference_metres(&self) -> f64 {
        self.atlas.circumference_metres()
    }

    #[wasm_bindgen(getter, js_name = regionCells)]
    pub fn region_cells(&self) -> u32 {
        relief::REGION_CELLS as u32
    }

    #[wasm_bindgen(getter, js_name = regionSpacingMetres)]
    pub fn region_spacing_metres(&self) -> f64 {
        self.level().spacing_metres()
    }

    /// Heights of the detailed region tile `(i, j)`, `regionCells` square.
    pub fn region(&mut self, i: i64, j: i64, rain_m_per_year: f64) -> Vec<f32> {
        let (seed, rock, level) = (self.seed, self.rock, self.level());
        let planet = self.planet();
        let region = relief::region(
            seed,
            level,
            (i, j),
            |x, y| planet.metres(x, y),
            rain_m_per_year,
            NARROWEST_RIVER_METRES,
            rock,
        );
        floats(region.field.land.height.iter().copied())
    }

    #[wasm_bindgen(getter, js_name = yearDays)]
    pub fn year_days(&self) -> f64 {
        earth::heavens().sky(self.recipe.orbit).year_seconds() / 86_400.0
    }

    /// Latitude where the sun stands overhead, in degrees, at a fraction of the year.
    #[wasm_bindgen(js_name = sunDeclination)]
    pub fn sun_declination(&self, year_fraction: f64) -> f64 {
        let heavens = earth::heavens();
        let sky = heavens.sky(self.recipe.orbit);
        sky.at(year_fraction * sky.year_seconds())
            .declination
            .to_degrees()
    }
}

impl Shape {
    fn planet(&mut self) -> &Field {
        self.relief
            .get_or_insert_with(|| relief::planet(self.seed, &self.atlas, self.rock))
    }

    fn level(&self) -> Level {
        Level {
            layer: 1,
            parent_cell_metres: self.atlas.cell_metres / relief::PLANET_REFINEMENT as f64,
            refinement: REGION_REFINEMENT,
        }
    }
}

#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layer {
    Elevation,
    Temperature,
    SummerTemperature,
    Rain,
    Ipse,
    Cover,
    River,
    Gyre,
    Evaporation,
    AetherWarming,
}

/// The settled climate of a world.
#[wasm_bindgen]
pub struct World {
    recipe: Recipe,
    map: GlobalMap,
}

#[wasm_bindgen]
impl World {
    #[wasm_bindgen(constructor)]
    pub fn new(seed: u64, width: u32, height: u32) -> World {
        let spec = spec(seed, width, height);
        let map = spec.generate(
            &earth::materials(),
            &earth::flora(),
            &earth::heavens(),
            &Sequential,
        );
        World {
            recipe: spec.recipe,
            map,
        }
    }

    #[wasm_bindgen(getter)]
    pub fn width(&self) -> u32 {
        self.map.width
    }

    #[wasm_bindgen(getter)]
    pub fn height(&self) -> u32 {
        self.map.height
    }

    /// Whether the climate and its winds reached a steady state.
    #[wasm_bindgen(getter)]
    pub fn settled(&self) -> bool {
        self.map.coupling.settled && self.map.energy_settled
    }

    #[wasm_bindgen(getter, js_name = seaSalinity)]
    pub fn sea_salinity(&self) -> f64 {
        self.map.sea_salinity_grams_per_kg
    }

    #[wasm_bindgen(getter, js_name = windMetresPerSecond)]
    pub fn wind_metres_per_second(&self) -> f64 {
        self.map.coupling.wind_metres_per_second
    }

    /// Cover codes: 0 sea, 1 sea ice, 2 land ice, 3 barren, 4 + n the n-th plant.
    pub fn layer(&self, layer: Layer) -> Vec<f32> {
        floats(self.map.cells.iter().map(|cell| match layer {
            Layer::Elevation => self.recipe.surface_metres(cell.terrain),
            Layer::Temperature => cell.temperature,
            Layer::SummerTemperature => cell.summer_temperature,
            Layer::Rain => cell.water.rain_mm_per_year,
            Layer::Ipse => cell.ipse.vapor + cell.ipse.condensed,
            Layer::Cover => match cell.cover {
                Cover::Sea => 0.0,
                Cover::SeaIce => 1.0,
                Cover::Ice => 2.0,
                Cover::Barren => 3.0,
                Cover::Plant(id) => (4 + id) as f64,
            },
            Layer::River => cell.water.river_mm_per_year,
            Layer::Gyre => cell.gyre_sverdrups,
            Layer::Evaporation => cell.water.evaporation_mm_per_year,
            Layer::AetherWarming => cell.aether_warming_kelvin,
        }))
    }

    /// Direction of the surface wind at a latitude, as `[east, north]`.
    #[wasm_bindgen(js_name = surfaceWind)]
    pub fn surface_wind(&self, north_degrees: f64) -> Vec<f64> {
        let flow = wind::layer_flow(0, -north_degrees / 90.0, &self.map.belts);
        vec![flow.east, -flow.south]
    }
}
