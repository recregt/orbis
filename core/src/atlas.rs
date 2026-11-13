// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec::Vec;
use serde::{Deserialize, Serialize};

use crate::life::{Cover, Flora};
use crate::world::{GlobalMap, Recipe, Terrain};

pub const PERSON_HEIGHT_METRES: f64 = 1.7;
const MILLIMETRES_PER_METRE: f64 = 1000.0;
const METRES_PER_KILOMETRE: f64 = 1000.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Ground {
    Water,
    Sand,
    Grass,
    Snow,
}

impl Ground {
    pub const ALL: [Ground; 4] = [Ground::Water, Ground::Sand, Ground::Grass, Ground::Snow];

    pub fn level(self) -> usize {
        self as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Place {
    pub ground: Ground,
    pub forest: bool,
    pub celsius: f64,
    pub summer_celsius: f64,
    pub metres: f64,
    pub rain_m_per_year: f64,
    pub ipse: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atlas {
    pub width: u32,
    pub height: u32,
    pub cell_metres: f64,
    pub row_metres: f64,
    pub places: Vec<Place>,
}

pub struct Blend {
    pub weights: [f64; 4],
    pub forest: f64,
    pub metres: f64,
    pub rain_m_per_year: f64,
    pub celsius: f64,
    pub summer_celsius: f64,
    pub ipse: f64,
}

impl Atlas {
    pub fn of(map: &GlobalMap, flora: &Flora, recipe: &Recipe) -> Self {
        let plants = flora.plants();
        let places = map
            .cells
            .iter()
            .map(|cell| {
                let (ground, forest) = match cell.cover {
                    Cover::Sea => (Ground::Water, false),
                    Cover::SeaIce | Cover::Ice => (Ground::Snow, false),
                    Cover::Barren => (Ground::Sand, false),
                    Cover::Plant(id) => (
                        Ground::Grass,
                        plants[id].height_metres > PERSON_HEIGHT_METRES,
                    ),
                };
                Place {
                    ground,
                    forest,
                    celsius: cell.temperature,
                    summer_celsius: cell.summer_temperature,
                    metres: recipe.surface_metres(cell.terrain),
                    rain_m_per_year: cell.water.rain_mm_per_year / MILLIMETRES_PER_METRE,
                    ipse: cell.ipse.vapor + cell.ipse.condensed,
                }
            })
            .collect();
        Self::laid_out(map.width, map.height, recipe, places)
    }

    pub fn of_terrain(width: u32, height: u32, terrains: &[Terrain], recipe: &Recipe) -> Self {
        let places = terrains
            .iter()
            .map(|&terrain| Place {
                ground: match terrain {
                    Terrain::Ocean { .. } => Ground::Water,
                    Terrain::Land { .. } => Ground::Sand,
                },
                forest: false,
                celsius: 0.0,
                summer_celsius: 0.0,
                metres: recipe.surface_metres(terrain),
                rain_m_per_year: 0.0,
                ipse: 0.0,
            })
            .collect();
        Self::laid_out(width, height, recipe, places)
    }

    fn laid_out(width: u32, height: u32, recipe: &Recipe, places: Vec<Place>) -> Self {
        let equator_metres = recipe.equator_km * METRES_PER_KILOMETRE;
        Self {
            width,
            height,
            cell_metres: equator_metres / f64::from(width),
            row_metres: equator_metres / 2.0 / f64::from(height),
            places,
        }
    }

    pub fn circumference_metres(&self) -> f64 {
        self.cell_metres * f64::from(self.width)
    }

    pub fn meridian_metres(&self) -> f64 {
        self.row_metres * f64::from(self.height)
    }

    pub fn place(&self, column: i64, row: i64) -> &Place {
        let column = column.rem_euclid(i64::from(self.width));
        let row = row.clamp(0, i64::from(self.height) - 1);
        &self.places[(row * i64::from(self.width) + column) as usize]
    }

    pub fn centre_metres(&self, column: i64, row: i64) -> (f64, f64) {
        (
            (column as f64 + 0.5) * self.cell_metres,
            (row as f64 + 0.5) * self.row_metres,
        )
    }

    pub fn blend(&self, x_metres: f64, y_metres: f64) -> Blend {
        let fx = x_metres / self.cell_metres - 0.5;
        let fy = y_metres / self.row_metres - 0.5;
        let (column, row) = (libm::floor(fx), libm::floor(fy));
        let (tx, ty) = (fx - column, fy - row);
        let (column, row) = (column as i64, row as i64);
        let mut weights = [0.0; 4];
        let mut forest = 0.0;
        let mut metres = 0.0;
        let mut rain_m_per_year = 0.0;
        let mut celsius = 0.0;
        let mut summer_celsius = 0.0;
        let mut ipse = 0.0;
        for (dx, dy, share) in [
            (0, 0, (1.0 - tx) * (1.0 - ty)),
            (1, 0, tx * (1.0 - ty)),
            (0, 1, (1.0 - tx) * ty),
            (1, 1, tx * ty),
        ] {
            let place = self.place(column + dx, row + dy);
            weights[place.ground.level()] += share;
            if place.forest {
                forest += share;
            }
            metres += share * place.metres;
            rain_m_per_year += share * place.rain_m_per_year;
            celsius += share * place.celsius;
            summer_celsius += share * place.summer_celsius;
            ipse += share * place.ipse;
        }
        Blend {
            weights,
            forest,
            metres,
            rain_m_per_year,
            celsius,
            summer_celsius,
            ipse,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn place(ground: Ground) -> Place {
        Place {
            ground,
            forest: false,
            celsius: 10.0,
            summer_celsius: 16.0,
            metres: 0.0,
            rain_m_per_year: 1.0,
            ipse: 0.0,
        }
    }

    fn strip() -> Atlas {
        Atlas {
            width: 2,
            height: 1,
            cell_metres: 1000.0,
            row_metres: 1000.0,
            places: alloc::vec![place(Ground::Water), place(Ground::Grass)],
        }
    }

    #[test]
    fn a_cell_centre_takes_that_cells_ground_alone() {
        let atlas = strip();
        let (x, y) = atlas.centre_metres(1, 0);
        let blend = atlas.blend(x, y);
        assert_eq!(blend.weights, [0.0, 0.0, 1.0, 0.0]);
    }

    #[test]
    fn halfway_between_two_cells_each_counts_half_and_the_planet_wraps() {
        let atlas = strip();
        let halfway = atlas.blend(1000.0, 500.0);
        assert_eq!(halfway.weights, [0.5, 0.0, 0.5, 0.0]);
        let across_the_seam = atlas.blend(0.0, 500.0);
        assert_eq!(across_the_seam.weights, [0.5, 0.0, 0.5, 0.0]);
        let shifted = atlas.blend(250.0 + atlas.circumference_metres(), 500.0);
        assert_eq!(shifted.weights, atlas.blend(250.0, 500.0).weights);
    }
}
