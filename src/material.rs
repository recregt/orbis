use alloc::string::String;
use alloc::vec::Vec;
use core::fmt;

use serde::Deserialize;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    Solid,
    Liquid,
    Gas,
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Threshold {
    pub celsius: f64,
    pub at: Option<f64>,
    #[serde(default)]
    pub per_doubling: f64,
}

impl Threshold {
    pub const fn fixed(celsius: f64) -> Self {
        Self {
            celsius,
            at: None,
            per_doubling: 0.0,
        }
    }

    pub const fn shifting(celsius: f64, at: f64, per_doubling: f64) -> Self {
        Self {
            celsius,
            at: Some(at),
            per_doubling,
        }
    }

    fn reference(&self) -> Option<f64> {
        self.at.filter(|_| self.per_doubling != 0.0)
    }

    pub fn celsius_at(&self, concentration: f64) -> f64 {
        match self.reference() {
            Some(at) => self.celsius + self.per_doubling * libm::log2(concentration / at),
            None => self.celsius,
        }
    }

    pub fn concentration_at(&self, celsius: f64) -> f64 {
        match self.reference() {
            Some(at) => at * libm::exp2((celsius - self.celsius) / self.per_doubling),
            None if celsius >= self.celsius => f64::INFINITY,
            None => 0.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    pub name: String,
    pub density: f64,
    pub melts_at: Threshold,
    pub boils_at: Threshold,
    #[serde(default)]
    pub solubility: f64,
    #[serde(default)]
    pub holds_water_mm: f64,
    #[serde(default)]
    pub specific_heat: f64,
    #[serde(default)]
    pub conductivity: f64,
    #[serde(default)]
    pub latent_heat: f64,
    #[serde(default)]
    pub expands_per_kelvin: f64,
    #[serde(default)]
    pub aether_absorptance: f64,
    #[serde(default)]
    pub longwave_absorption_m2_per_kg: f64,
    #[serde(default)]
    pub sunlight_absorptance: f64,
    #[serde(default)]
    pub longwave_window_share: f64,
    #[serde(default)]
    pub free_air_dissipation_share: f64,
    #[serde(default)]
    pub erodibility_per_year: f64,
    #[serde(default)]
    pub creep_m2_per_year: f64,
    #[serde(default)]
    pub relief_erosion_per_year: f64,
}

impl Material {
    pub fn thermal_diffusivity(&self) -> f64 {
        self.conductivity / (self.density * self.specific_heat)
    }

    pub fn saturation(&self, celsius: f64) -> f64 {
        self.boils_at.concentration_at(celsius)
    }

    pub fn saturation_per_kelvin(&self, celsius: f64) -> f64 {
        self.saturation_with_slope(celsius).1
    }

    pub fn saturation_with_slope(&self, celsius: f64) -> (f64, f64) {
        let saturation = self.saturation(celsius);
        match self.boils_at.reference() {
            Some(_) => (
                saturation,
                saturation * core::f64::consts::LN_2 / self.boils_at.per_doubling,
            ),
            None => (saturation, 0.0),
        }
    }

    pub fn condensed_phase(&self, celsius: f64, amount: f64) -> Phase {
        if celsius < self.melts_at.celsius_at(amount) {
            Phase::Solid
        } else {
            Phase::Liquid
        }
    }

    pub fn pure_phase(&self, celsius: f64) -> Phase {
        if self.saturation(celsius) >= self.density {
            Phase::Gas
        } else {
            self.condensed_phase(celsius, self.density)
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContentError {
    DuplicateMaterial(String),
    NonPositiveDensity(String),
    FallingThreshold(String),
    MissingReference(String),
    ShareOutOfRange {
        material: String,
        property: &'static str,
    },
    NegativeWaterHolding(String),
    NegativeHeatProperty(String),
    MissingMaterial(&'static str),
    MissingHeatProperty(&'static str),
    MissingProperty {
        material: &'static str,
        property: &'static str,
        needed_by: &'static str,
    },
}

impl fmt::Display for ContentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateMaterial(name) => write!(f, "material `{name}` is defined twice"),
            Self::NonPositiveDensity(name) => {
                write!(f, "material `{name}`: density must be positive")
            }
            Self::FallingThreshold(name) => write!(
                f,
                "material `{name}`: thresholds may only rise with concentration (per_doubling >= 0)"
            ),
            Self::MissingReference(name) => write!(
                f,
                "material `{name}`: a threshold that shifts needs a positive reference concentration `at`"
            ),
            Self::ShareOutOfRange { material, property } => write!(
                f,
                "material `{material}`: {property} is a share and must lie between 0 and 1"
            ),
            Self::NegativeWaterHolding(name) => {
                write!(f, "material `{name}`: holds_water_mm cannot be negative")
            }
            Self::NegativeHeatProperty(name) => write!(
                f,
                "material `{name}`: specific_heat and conductivity cannot be negative"
            ),
            Self::MissingMaterial(name) => {
                write!(f, "material `{name}` is required by the world's laws")
            }
            Self::MissingProperty {
                material,
                property,
                needed_by,
            } => write!(
                f,
                "material `{material}` needs a positive {property} for {needed_by}"
            ),
            Self::MissingHeatProperty(name) => write!(
                f,
                "material `{name}` needs a positive specific_heat and conductivity for the energy balance"
            ),
        }
    }
}

pub struct Materials {
    pub air: Material,
    pub water: Material,
    pub ipse: Material,
    pub rock: Material,
    all: Vec<Material>,
}

impl Materials {
    pub fn new(all: Vec<Material>) -> Result<Self, Vec<ContentError>> {
        let mut errors = Vec::new();
        for (i, material) in all.iter().enumerate() {
            let name = || material.name.clone();
            if all[..i].iter().any(|m| m.name == material.name) {
                errors.push(ContentError::DuplicateMaterial(name()));
            }
            if material.density <= 0.0 {
                errors.push(ContentError::NonPositiveDensity(name()));
            }
            let thresholds = [material.melts_at, material.boils_at];
            if thresholds.iter().any(|t| t.per_doubling < 0.0) {
                errors.push(ContentError::FallingThreshold(name()));
            }
            if thresholds
                .iter()
                .any(|t| t.per_doubling != 0.0 && t.at.is_none_or(|at| at <= 0.0))
            {
                errors.push(ContentError::MissingReference(name()));
            }
            for (share, property) in [
                (material.solubility, "solubility"),
                (material.aether_absorptance, "aether_absorptance"),
                (material.sunlight_absorptance, "sunlight_absorptance"),
                (material.longwave_window_share, "longwave_window_share"),
                (
                    material.free_air_dissipation_share,
                    "free_air_dissipation_share",
                ),
            ] {
                if !(0.0..=1.0).contains(&share) {
                    errors.push(ContentError::ShareOutOfRange {
                        material: name(),
                        property,
                    });
                }
            }
            if material.holds_water_mm < 0.0 {
                errors.push(ContentError::NegativeWaterHolding(name()));
            }
            if material.specific_heat < 0.0 || material.conductivity < 0.0 {
                errors.push(ContentError::NegativeHeatProperty(name()));
            }
        }
        let find = |name: &'static str, errors: &mut Vec<ContentError>| {
            let found = all.iter().find(|m| m.name == name).cloned();
            if found.is_none() {
                errors.push(ContentError::MissingMaterial(name));
            }
            found
        };
        let air = find("air", &mut errors);
        let water = find("water", &mut errors);
        let ipse = find("ipse", &mut errors);
        let rock = find("rock", &mut errors);
        for (name, material) in [("air", &air), ("water", &water), ("rock", &rock)] {
            if material
                .as_ref()
                .is_some_and(|m| m.specific_heat <= 0.0 || m.conductivity <= 0.0)
            {
                errors.push(ContentError::MissingHeatProperty(name));
            }
        }
        let mut required = Vec::new();
        if let Some(water) = &water {
            required.extend([
                (
                    "water",
                    water.latent_heat,
                    "latent_heat",
                    "the moist energy balance",
                ),
                (
                    "water",
                    water.expands_per_kelvin,
                    "expands_per_kelvin",
                    "the equatorial ocean width",
                ),
                (
                    "water",
                    water.longwave_absorption_m2_per_kg,
                    "longwave_absorption_m2_per_kg",
                    "the radiation through the air",
                ),
            ]);
        }
        if let Some(air) = &air {
            required.extend([
                (
                    "air",
                    air.longwave_absorption_m2_per_kg,
                    "longwave_absorption_m2_per_kg",
                    "the radiation through the air",
                ),
                (
                    "air",
                    air.sunlight_absorptance,
                    "sunlight_absorptance",
                    "the radiation through the air",
                ),
                (
                    "air",
                    air.longwave_window_share,
                    "longwave_window_share",
                    "the radiation through the air",
                ),
                (
                    "air",
                    air.free_air_dissipation_share,
                    "free_air_dissipation_share",
                    "the wind's energy budget",
                ),
            ]);
        }
        if let Some(rock) = &rock {
            required.extend([
                (
                    "rock",
                    rock.erodibility_per_year,
                    "erodibility_per_year",
                    "the erosion of the land",
                ),
                (
                    "rock",
                    rock.creep_m2_per_year,
                    "creep_m2_per_year",
                    "the erosion of the land",
                ),
                (
                    "rock",
                    rock.relief_erosion_per_year,
                    "relief_erosion_per_year",
                    "the uplift of the land",
                ),
            ]);
        }
        for (material, value, property, needed_by) in required {
            if value <= 0.0 {
                errors.push(ContentError::MissingProperty {
                    material,
                    property,
                    needed_by,
                });
            }
        }
        match (air, water, ipse, rock) {
            (Some(air), Some(water), Some(ipse), Some(rock)) if errors.is_empty() => Ok(Self {
                air,
                water,
                ipse,
                rock,
                all,
            }),
            _ => Err(errors),
        }
    }

    pub fn all(&self) -> &[Material] {
        &self.all
    }
}

#[cfg(test)]
pub(crate) fn air() -> Material {
    Material {
        name: "air".into(),
        density: 1.2,
        melts_at: Threshold::fixed(-213.0),
        boils_at: Threshold::fixed(-194.0),
        solubility: 0.0,
        holds_water_mm: 0.0,
        specific_heat: 1005.0,
        conductivity: 0.026,
        latent_heat: 0.0,
        expands_per_kelvin: 0.0,
        aether_absorptance: 0.0,
        longwave_absorption_m2_per_kg: 1.24e-4,
        sunlight_absorptance: 0.23,
        longwave_window_share: 0.1,
        free_air_dissipation_share: 0.5,
        erodibility_per_year: 0.0,
        creep_m2_per_year: 0.0,
        relief_erosion_per_year: 0.0,
    }
}

#[cfg(test)]
pub(crate) fn water() -> Material {
    Material {
        name: "water".into(),
        density: 1000.0,
        melts_at: Threshold::fixed(0.0),
        boils_at: Threshold::shifting(20.0, 0.0173, 11.0),
        solubility: 0.0,
        holds_water_mm: 0.0,
        specific_heat: 4186.0,
        conductivity: 0.6,
        latent_heat: 2.5e6,
        expands_per_kelvin: 2.07e-4,
        aether_absorptance: 0.0,
        longwave_absorption_m2_per_kg: 0.324,
        sunlight_absorptance: 0.0,
        longwave_window_share: 0.0,
        free_air_dissipation_share: 0.0,
        erodibility_per_year: 0.0,
        creep_m2_per_year: 0.0,
        relief_erosion_per_year: 0.0,
    }
}

#[cfg(test)]
pub(crate) fn ipse() -> Material {
    Material {
        name: "ipse".into(),
        density: 4200.0,
        melts_at: Threshold::shifting(20.0, 120.0, 40.0),
        boils_at: Threshold::shifting(20.0, 40.0, 40.0),
        solubility: 0.0,
        holds_water_mm: 0.0,
        specific_heat: 0.0,
        conductivity: 0.0,
        latent_heat: 0.0,
        expands_per_kelvin: 0.0,
        aether_absorptance: 1.0,
        longwave_absorption_m2_per_kg: 0.0,
        sunlight_absorptance: 0.0,
        longwave_window_share: 0.0,
        free_air_dissipation_share: 0.0,
        erodibility_per_year: 0.0,
        creep_m2_per_year: 0.0,
        relief_erosion_per_year: 0.0,
    }
}

#[cfg(test)]
pub(crate) fn rock() -> Material {
    Material {
        name: "rock".into(),
        density: 2700.0,
        melts_at: Threshold::fixed(1215.0),
        boils_at: Threshold::fixed(2230.0),
        solubility: 0.0001,
        holds_water_mm: 150.0,
        specific_heat: 790.0,
        conductivity: 2.8,
        latent_heat: 0.0,
        expands_per_kelvin: 0.0,
        aether_absorptance: 0.0,
        longwave_absorption_m2_per_kg: 0.0,
        sunlight_absorptance: 0.0,
        longwave_window_share: 0.0,
        free_air_dissipation_share: 0.0,
        erodibility_per_year: 1.0e-5,
        creep_m2_per_year: 0.01,
        relief_erosion_per_year: 1.535e-7,
    }
}

#[cfg(test)]
pub(crate) fn test_materials() -> Materials {
    Materials::new(alloc::vec![air(), water(), ipse(), rock()])
        .unwrap_or_else(|errors| panic!("test materials are valid: {errors:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn warm_air_holds_more_vapor() {
        let water = water();
        assert!(water.saturation(30.0) > water.saturation(20.0));
        assert_eq!(water.saturation(20.0), 0.0173);
        assert_eq!(water.saturation(31.0), 0.0173 * 2.0);
    }

    #[test]
    fn pure_substances_keep_their_everyday_phase() {
        let water = water();
        assert_eq!(water.pure_phase(-5.0), Phase::Solid);
        assert_eq!(water.pure_phase(20.0), Phase::Liquid);
        assert_eq!(ipse().pure_phase(20.0), Phase::Solid);
    }

    #[test]
    fn collected_ipse_crystallizes_and_cold_lowers_the_amount_needed() {
        let ipse = ipse();
        assert_eq!(ipse.condensed_phase(20.0, 60.0), Phase::Liquid);
        assert_eq!(ipse.condensed_phase(20.0, 200.0), Phase::Solid);
        assert_eq!(ipse.condensed_phase(-20.0, 80.0), Phase::Solid);
    }

    #[test]
    fn content_rules_are_enforced() {
        let mut weightless = water();
        weightless.density = 0.0;
        let mut falling = ipse();
        falling.boils_at.per_doubling = -1.0;
        let mut unanchored = ipse();
        unanchored.melts_at.at = None;
        let mut overly_soluble = rock();
        overly_soluble.solubility = 2.0;
        let errors = Materials::new(alloc::vec![weightless, falling, unanchored, overly_soluble])
            .err()
            .unwrap();
        assert_eq!(
            errors,
            [
                ContentError::NonPositiveDensity("water".into()),
                ContentError::FallingThreshold("ipse".into()),
                ContentError::DuplicateMaterial("ipse".into()),
                ContentError::MissingReference("ipse".into()),
                ContentError::ShareOutOfRange {
                    material: "rock".into(),
                    property: "solubility",
                },
                ContentError::MissingMaterial("air"),
            ]
        );
    }

    #[test]
    fn laws_require_their_materials() {
        let errors = Materials::new(alloc::vec![water()]).err().unwrap();
        assert_eq!(
            errors,
            [
                ContentError::MissingMaterial("air"),
                ContentError::MissingMaterial("ipse"),
                ContentError::MissingMaterial("rock"),
            ]
        );
    }

    #[test]
    fn radiation_needs_the_air_and_vapor_to_absorb() {
        let mut clear_air = air();
        clear_air.longwave_absorption_m2_per_kg = 0.0;
        clear_air.sunlight_absorptance = 0.0;
        clear_air.longwave_window_share = 0.0;
        let mut clear_water = water();
        clear_water.longwave_absorption_m2_per_kg = 0.0;
        let errors = Materials::new(alloc::vec![clear_air, clear_water, ipse(), rock()])
            .err()
            .unwrap();
        let missing = |material, property| ContentError::MissingProperty {
            material,
            property,
            needed_by: "the radiation through the air",
        };
        assert_eq!(
            errors,
            [
                missing("water", "longwave_absorption_m2_per_kg"),
                missing("air", "longwave_absorption_m2_per_kg"),
                missing("air", "sunlight_absorptance"),
                missing("air", "longwave_window_share"),
            ]
        );
    }

    #[test]
    fn the_energy_balance_needs_heat_properties() {
        let mut insulating_rock = rock();
        insulating_rock.conductivity = 0.0;
        let mut strange_air = air();
        strange_air.specific_heat = -1.0;
        let errors = Materials::new(alloc::vec![strange_air, water(), ipse(), insulating_rock])
            .err()
            .unwrap();
        assert_eq!(
            errors,
            [
                ContentError::NegativeHeatProperty("air".into()),
                ContentError::MissingHeatProperty("air"),
                ContentError::MissingHeatProperty("rock"),
            ]
        );
    }
}
