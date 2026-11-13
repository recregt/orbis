// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use orbis::Seed;
use orbis::atlas::Atlas;
use orbis::earth;
use orbis::life::Cover;
use orbis::parallel::Sequential;
use orbis::relief::{self, Field, Level, Rock};
use orbis::world::{GlobalMap, Recipe, WorldSpec};
use wasm_bindgen::prelude::*;

const REGION_REFINEMENT: usize = 8;
const NARROWEST_RIVER_METRES: f64 = 1.0;

#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layer {
    Elevation,
    Temperature,
    SummerTemperature,
    Rain,
    Ipse,
    Cover,
}

#[wasm_bindgen]
pub struct World {
    seed: Seed,
    recipe: Recipe,
    rock: Rock,
    map: GlobalMap,
    atlas: Atlas,
    relief: Option<Field>,
}

#[wasm_bindgen]
impl World {
    #[wasm_bindgen(constructor)]
    pub fn new(seed: u64, width: u32, height: u32) -> World {
        let (materials, flora, heavens) = (earth::materials(), earth::flora(), earth::heavens());
        let spec = WorldSpec {
            seed: Seed(seed),
            width,
            height,
            recipe: earth::recipe(),
        };
        let map = spec.generate(&materials, &flora, &heavens, &Sequential);
        let atlas = Atlas::of(&map, &flora, &spec.recipe);
        World {
            seed: spec.seed,
            recipe: spec.recipe,
            rock: Rock::of(&materials.rock),
            map,
            atlas,
            relief: None,
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

    pub fn layer(&self, layer: Layer) -> Vec<f32> {
        self.map
            .cells
            .iter()
            .map(|cell| match layer {
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
            } as f32)
            .collect()
    }

    pub fn relief(&mut self) -> Vec<f32> {
        heights(self.planet())
    }

    #[wasm_bindgen(getter, js_name = reliefWidth)]
    pub fn relief_width(&self) -> u32 {
        self.map.width * relief::PLANET_REFINEMENT as u32
    }

    #[wasm_bindgen(getter, js_name = reliefHeight)]
    pub fn relief_height(&self) -> u32 {
        self.map.height * relief::PLANET_REFINEMENT as u32
    }

    #[wasm_bindgen(getter, js_name = regionCells)]
    pub fn region_cells(&self) -> u32 {
        relief::REGION_CELLS as u32
    }

    #[wasm_bindgen(getter, js_name = regionSpacingMetres)]
    pub fn region_spacing_metres(&self) -> f64 {
        self.level().spacing_metres()
    }

    #[wasm_bindgen(getter, js_name = circumferenceMetres)]
    pub fn circumference_metres(&self) -> f64 {
        self.atlas.circumference_metres()
    }

    pub fn region(&mut self, i: i64, j: i64) -> Vec<f32> {
        let level = self.level();
        let spacing = level.spacing_metres();
        let (seed, rock) = (self.seed, self.rock);
        let rain = self
            .atlas
            .blend(i as f64 * spacing, j as f64 * spacing)
            .rain_m_per_year;
        let planet = self.planet();
        let region = relief::region(
            seed,
            level,
            (i, j),
            |x, y| planet.metres(x, y),
            rain,
            NARROWEST_RIVER_METRES,
            rock,
        );
        heights(&region.field)
    }
}

impl World {
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

fn heights(field: &Field) -> Vec<f32> {
    field.land.height.iter().map(|&h| h as f32).collect()
}
