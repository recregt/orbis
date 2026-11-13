// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec;
use alloc::vec::Vec;

use crate::life::{Flora, Plant, Thrives, Tolerance};
use crate::material::{Material, Materials, Threshold};
use crate::sky::{Heavens, Orbit, Star};
use crate::world::Recipe;

pub fn air() -> Material {
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

pub fn water() -> Material {
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

pub fn ipse() -> Material {
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

pub fn rock() -> Material {
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

pub fn materials() -> Materials {
    Materials::new(vec![air(), water(), ipse(), rock()])
        .unwrap_or_else(|errors| panic!("the earth materials are valid: {errors:?}"))
}

pub fn flora() -> Flora {
    Flora::new(vec![
        Plant {
            name: "grass".into(),
            height_metres: 0.5,
            thrives: Thrives {
                celsius: Tolerance::above(-5.0),
                ..Thrives::default()
            },
        },
        Plant {
            name: "tree".into(),
            height_metres: 20.0,
            thrives: Thrives {
                celsius: Tolerance::above(0.0),
                water_mm_per_year: Tolerance::above(500.0),
                ..Thrives::default()
            },
        },
    ])
    .unwrap_or_else(|errors| panic!("the earth flora is valid: {errors:?}"))
}

pub fn sun() -> Star {
    Star {
        name: "sun".into(),
        luminosity_watts: 3.828e26,
        mass_kg: 1.988_47e30,
    }
}

pub fn heavens() -> Heavens {
    Heavens {
        star: sun(),
        companions: Vec::new(),
    }
}

pub const ORBIT: Orbit = Orbit {
    semi_major_axis_metres: 1.495_978_707e11,
    eccentricity: 0.0167,
    axial_tilt_degrees: 23.44,
    sun_longitude_at_perihelion_degrees: 282.9,
};

pub fn recipe() -> Recipe {
    Recipe {
        equator_km: 40_075.0,
        sidereal_day_hours: 23.934,
        orbit: ORBIT,
        land_fraction: 0.3,
        continent_frequency: 2.0,
        roughness: 0.5,
        peak_height_metres: 6000.0,
        deepest_metres: 11_000.0,
        weather_systems_around: 16.0,
        weather_celsius: 6.0,
        climate_downscale: 4,
        salt_residence_years: 13.0e6,
        eddy_diffusivity_m2_per_s: 1.0e6,
        ipse_age_days: 3.6,
        ipse_mean_amount: 10.0,
    }
}
