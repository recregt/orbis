// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::string::String;
use alloc::vec::Vec;
use core::fmt;

use serde::Deserialize;

#[derive(Debug, Clone, Copy, PartialEq, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Tolerance {
    pub min: Option<f64>,
    pub max: Option<f64>,
}

impl Tolerance {
    pub const fn above(min: f64) -> Self {
        Self {
            min: Some(min),
            max: None,
        }
    }

    pub fn admits(&self, value: f64) -> bool {
        self.min.is_none_or(|min| value >= min) && self.max.is_none_or(|max| value <= max)
    }

    fn is_ordered(&self) -> bool {
        match (self.min, self.max) {
            (Some(min), Some(max)) => min <= max,
            _ => true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Thrives {
    #[serde(default)]
    pub celsius: Tolerance,
    #[serde(default)]
    pub water_mm_per_year: Tolerance,
    #[serde(default)]
    pub ipse: Tolerance,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Plant {
    pub name: String,
    pub height_metres: f64,
    pub thrives: Thrives,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Site {
    pub celsius: f64,
    pub water_mm_per_year: f64,
    pub ipse: f64,
}

impl Plant {
    pub fn survives(&self, site: &Site) -> bool {
        self.thrives.celsius.admits(site.celsius)
            && self
                .thrives
                .water_mm_per_year
                .admits(site.water_mm_per_year)
            && self.thrives.ipse.admits(site.ipse)
    }
}

pub type PlantId = usize;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cover {
    Sea,
    SeaIce,
    Ice,
    Plant(PlantId),
    Barren,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PlantError {
    DuplicatePlant(String),
    NonPositiveHeight(String),
    UnorderedTolerance(String),
}

impl fmt::Display for PlantError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicatePlant(name) => write!(f, "plant `{name}` is defined twice"),
            Self::NonPositiveHeight(name) => write!(f, "plant `{name}`: height must be positive"),
            Self::UnorderedTolerance(name) => {
                write!(
                    f,
                    "plant `{name}`: every tolerance `min` must not exceed its `max`"
                )
            }
        }
    }
}

pub struct Flora {
    plants: Vec<Plant>,
}

impl Flora {
    pub fn new(plants: Vec<Plant>) -> Result<Self, Vec<PlantError>> {
        let mut errors = Vec::new();
        for (i, plant) in plants.iter().enumerate() {
            let name = || plant.name.clone();
            if plants[..i].iter().any(|p| p.name == plant.name) {
                errors.push(PlantError::DuplicatePlant(name()));
            }
            if plant.height_metres <= 0.0 {
                errors.push(PlantError::NonPositiveHeight(name()));
            }
            let t = &plant.thrives;
            if ![t.celsius, t.water_mm_per_year, t.ipse]
                .iter()
                .all(Tolerance::is_ordered)
            {
                errors.push(PlantError::UnorderedTolerance(name()));
            }
        }
        if errors.is_empty() {
            Ok(Self { plants })
        } else {
            Err(errors)
        }
    }

    pub fn dominant(&self, site: &Site) -> Option<PlantId> {
        self.plants
            .iter()
            .enumerate()
            .filter(|(_, plant)| plant.survives(site))
            .max_by(|(a, pa), (b, pb)| pa.height_metres.total_cmp(&pb.height_metres).then(b.cmp(a)))
            .map(|(id, _)| id)
    }

    pub fn plants(&self) -> &[Plant] {
        &self.plants
    }
}

#[cfg(test)]
pub(crate) fn plant(name: &str, height_metres: f64, thrives: Thrives) -> Plant {
    Plant {
        name: name.into(),
        height_metres,
        thrives,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::earth;

    fn site(celsius: f64, water_mm_per_year: f64) -> Site {
        Site {
            celsius,
            water_mm_per_year,
            ipse: 0.0,
        }
    }

    #[test]
    fn the_tallest_plant_that_survives_wins() {
        let flora = earth::flora();
        assert_eq!(flora.dominant(&site(15.0, 800.0)), Some(1));
        assert_eq!(flora.dominant(&site(15.0, 200.0)), Some(0));
        assert_eq!(flora.dominant(&site(-10.0, 800.0)), None);
    }

    #[test]
    fn equally_tall_plants_resolve_by_order() {
        let flora = Flora::new(alloc::vec![
            plant("first", 1.0, Thrives::default()),
            plant("second", 1.0, Thrives::default()),
        ])
        .unwrap();
        assert_eq!(flora.dominant(&site(0.0, 0.0)), Some(0));
    }

    #[test]
    fn plant_sheets_are_checked() {
        let backwards = Thrives {
            celsius: Tolerance {
                min: Some(10.0),
                max: Some(0.0),
            },
            ..Thrives::default()
        };
        let errors = Flora::new(alloc::vec![
            plant("odd", 0.0, backwards),
            plant("odd", 1.0, Thrives::default()),
        ])
        .err()
        .unwrap();
        assert_eq!(
            errors,
            [
                PlantError::NonPositiveHeight("odd".into()),
                PlantError::UnorderedTolerance("odd".into()),
                PlantError::DuplicatePlant("odd".into()),
            ]
        );
    }
}
